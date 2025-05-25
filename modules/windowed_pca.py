'''
Window parser and GT PCA functions.
'''


## IMPORT CONFIG
from . import config

## IMPORT PACKAGES & SET NUMBER OF THREADS FOR PCANGSD
import os
if config.VAR_FMT == 'GT':
    import allel
else:
    from pcangsd.shared import emMAF
    from pcangsd.reader_cy import filterArrays                                  # pylint: disable=E0611
    from pcangsd.covariance import emPCA
    os.environ["OMP_NUM_THREADS"] = str(config.N_THREADS)
    os.environ["OPENBLAS_NUM_THREADS"] = str(config.N_THREADS)
    os.environ["MKL_NUM_THREADS"] = str(config.N_THREADS)
import sys
import gzip
import numpy as np

## MODULES
from modules.log import Log

## INSTANTIATE LOGGER
log = Log()


## CLASSES

class WPCA:

    '''
    Parse hard-called genotypes (GT) from VCF or TSV and apply a function.
    '''

    def __init__(self,
                 variant_file_path,
                 file_fmt,
                 var_fmt,
                 sample_lst,
                 n_pcs, pc_a, pc_b,
                 chrom, start, stop,
                 w_size, w_step,
                 gt_min_var_per_w,
                 gl_pl_min_var_per_w,
                 skip_monomorphic,
                 gt_mean_impute,
                 vcf_pass_filter,
                 proc_trail_trunc_w,
                 min_maf,
                 x_mode,
                 n_threads,
                 ):

        # variant file/region
        self.variant_file_path = variant_file_path
        self.file_fmt = file_fmt
        self.var_fmt = var_fmt
        self.sample_lst = sample_lst
        self.n_pcs = int(n_pcs)
        self.pc_a = int(pc_a)
        self.pc_b = int(pc_b)
        self.chrom = chrom
        self.start = start
        self.stop = stop

        # parameters
        self.w_size = w_size
        self.w_step = w_step
        self.gt_min_var_per_w = gt_min_var_per_w
        self.gl_pl_min_var_per_w = gl_pl_min_var_per_w
        self.skip_monomorphic = skip_monomorphic
        self.gt_mean_impute = gt_mean_impute
        self.vcf_pass_filter = vcf_pass_filter
        self.proc_trail_trunc_w = proc_trail_trunc_w
        self.min_maf = min_maf
        self.x_mode = x_mode
        self.n_threads = n_threads

        # transient variables
        self.n_windows = None
        self.w_start = None
        self.w_stop = None
        self.w_pos = None
        self.pos = None
        self.win = None
        self.w_gt_arr = None
        self.w_gl_arr = None
        self.w_pl_arr = None
        self.gl_min_maf_arr = None
        self.w_idx = None
        self.pct_complete = 0

        # instantiate results dict
        self.out_dct = {
        }

        # genotype encoding
        self.gt_code_dct = {
            '0/0':  0,
            '0|0':  0,
            '0/1':  1,
            '0|1':  1,
            '1/0':  1,
            '1|0':  1,
            '1/1':  2,
            '1|1':  2,
            './.': np.nan,
            '.|.': np.nan,
            '0/.': np.nan,
            '0|.': np.nan,
            './0': np.nan,
            '.|0': np.nan,
            '1/.': np.nan,
            '1|.': np.nan,
            './1': np.nan,
            '.|1': np.nan,
            '.':   np.nan,
        }


    def parse_variant(self, process_func, values):
        '''
        Parse GT/GL/PL values of a variant (fixed window size)
        '''
        while self.pos > self.w_stop:
            process_func()
            self.init_win()
        if self.pos >= self.stop:
            if (self.proc_trail_trunc_w
                and any([x[0] > self.w_stop for x in self.win])
            ):
                self.w_stop = self.stop
                process_func()
                self.init_win()
            return True
        if self.pos >= self.w_start: 
            self.win.append([self.pos] + values)
        return False


    def parse_variant_x(self, process_func, values):
        '''
        Parse GT/GL/PL values of a variant (variable window size)
        '''
        if len(self.win) < self.w_size:
            self.win.append([self.pos] + values)
        if len(self.win) == self.w_size:
            self.w_start = self.win[0][0]
            self.w_stop = self.win[-1][0]
            process_func()
            self.init_win_x()
        if self.pos >= self.stop:
            if self.proc_trail_trunc_w:
                if (len(self.win) > (self.w_size-self.w_step)):
                    self.w_start = self.win[0][0]
                    self.w_stop = self.stop
                    process_func()
            return True
        return False


    def init_win(self):
        '''
        Initialize new window by shifting one w_step and dropping obsolete
        variants from previous window.
        '''

        self.w_start = self.w_start + self.w_step
        self.w_stop = self.w_start + self.w_size-1
        self.w_idx += 1
        self.win = [x for x in self.win if x[0] >= self.w_start]

        log.info(f'Processed {self.w_idx}/{self.n_windows} windows')


    def gt_mean_imputation(self):
        '''
        Mean impute missing GT calls.
        '''

        # calculate row/site means (ignoring NaN)
        row_mean_arr = np.nanmean(self.w_gt_arr, axis=1)

        # fetch indices of rows/sites WITH missing call(s)
        miss_mask_arr = np.isnan(self.w_gt_arr)

        # replace NaN with row/site means
        self.w_gt_arr[miss_mask_arr] = np.take(
            row_mean_arr, np.where(miss_mask_arr)[0]
        )


    def init_win_x(self):
        '''
        Initialize new window by shifting one w_step and dropping obsolete
        variants from previous window.
        '''

        self.w_idx += 1
        self.win = self.win[self.w_step:]

        pct_complete = round((self.pos/self.stop)*100)
        if pct_complete > self.pct_complete and pct_complete % 10 == 0:
            self.pct_complete = pct_complete
            log.info(f'Processed {self.pct_complete}%')


    def gt_drop_missing_sites(self):
        '''
        remove any site with at least one missing GT call.
        '''

        # fetch indices of rows/sites WITHOUT missing call(s)
        mask = ~np.isnan(self.w_gt_arr).any(axis=1)

        # drop those rows/sites
        self.w_gt_arr = self.w_gt_arr[mask]


    def gt_min_maf_filter(self):
        '''
        Drop SNPs with minor allele frequency below specified value.
        '''

        # count # GTs per row/site (np.sum counts True), *2 because diploid
        allele_counts_arr = np.sum(~np.isnan(self.w_gt_arr), axis=1) * 2

        # count # alleles per row/site
        allele_sums_arr = np.nansum(self.w_gt_arr, axis=1)

        # calculate allel frequencies and multiple with -1 if AF > 0.5 (because
        # input data may not be polarized by major/minor allel)
        af_arr = allele_sums_arr / allele_counts_arr
        af_arr[af_arr > 0.5] = 1 - af_arr[af_arr > 0.5]

        # keep only sites where AF >= min_maf
        self.w_gt_arr = self.w_gt_arr[af_arr >= self.min_maf]


    def gt_process_win(self):
        '''
        Remove POS info, convert to numpy array and apply target function or
        return empty array if there are no variants.
        '''
        # non-empty: trim off pos info, convert to numpy arr, drop missing
        # sites or mean impute, apply min_maf filter
        if self.win:
            self.w_gt_arr = np.array(
                [x[1:] for x in self.win], dtype=np.float32
                )

        # empty: convert to empty numpy arr
        else:
            self.w_gt_arr = np.empty((0,0))

        # call target function
        self.pca()


    def gl_min_maf_filter(self):
        '''
        Drop variants with estimated minor allele frequency below specified
        value using PCAngsd code.
        '''

        # PCAngsd
        # emMAF expects different numbers of args depending on the version
        try:
            self.gl_min_maf_arr = emMAF(self.w_gl_arr, 200, 1e-4)
        except:
            self.gl_min_maf_arr = emMAF(
                self.w_gl_arr, 200, 1e-4, self.n_threads
            )
            log.info('PCAngsd version is outdated, updating to the latest '
                     'version is recommended')

        # filter by minor allel frequency
        if self.min_maf > 0.0:
            min_maf_mask = (self.gl_min_maf_arr >= self.min_maf) \
                & (self.gl_min_maf_arr <= (1 - self.min_maf))

            # Update arrays
            m = np.sum(min_maf_mask)
            tmp_mask = min_maf_mask.astype(np.uint8)
            filterArrays(self.w_gl_arr, self.gl_min_maf_arr, tmp_mask)
            self.w_gl_arr = self.w_gl_arr[:m,:]


    def gl_process_win(self):
        '''
        Remove POS info and convert to numpy array, return empty array if there
        are no variants.
        '''

        # mute STDOUT by redirecting STDOUT tp /dev/null
        old_stdout = sys.stdout
        sys.stdout = open(os.devnull, 'w')

        # non-empty: trim off pos info, convert to numpy arr, apply min_maf
        # filter, convert PL to GL
        if self.win:
            self.w_gl_arr = np.array([x[1:] for x in self.win],
                                     dtype=np.float32)
            self.gl_min_maf_filter()

        # empty: convert to empty numpy arr & apply target function
        else:
            self.w_gl_arr = np.empty((0,0))

        # call target function
        self.pcangsd()

        # restore sys.stdout
        sys.stdout = old_stdout


    def pl_convert_to_gl(self):
        '''
        Convert a 2D array of phred scaled genotype likelihoods (PL) to
        normalized genotype likelihoods (GL) and return a 2D array, omitting
        the third GL, which is the expected input for PCAngsd.
        '''

        # fetch dimensions
        n_rows, n_cols = self.w_pl_arr.shape

        # reshape pl_arr to have separation by sample as first dimension
        # (to vectorize normalization) --> pl_arr_3d dimensions: (samples,'
        # variants, 3 pl_values)
        w_pl_arr_3d = self.w_pl_arr.reshape(n_rows, -1, 3).transpose(1, 0, 2)   # pylint: disable=E1121

        # unphred
        w_pl_arr_3d = np.power(10, -w_pl_arr_3d/10)

        # preallocate output w_gl_arr (PCAngsd expects only first two GL values
        # for memory efficiency)
        self.w_gl_arr = np.zeros((n_rows, n_cols*2//3), dtype=np.float32)

        # normalize all GLs partitioned by sample
        for idx in range(0, w_pl_arr_3d.shape[0]):
            ind_arr = w_pl_arr_3d[idx]
            ind_arr = 1-(ind_arr / ind_arr.sum(axis=1).reshape(n_rows, 1))
            self.w_gl_arr[:,idx*2], self.w_gl_arr[:,idx*2+1] = \
                ind_arr[:,0], ind_arr[:,1]


    def pl_process_win(self):
        '''
        Remove POS info and convert to numpy array, return empty array if there
        are no variants.
        '''

        # mute STDOUT by redirecting STDOUT tp /dev/null
        old_stdout = sys.stdout
        sys.stdout = open(os.devnull, 'w')

        # non-empty: trim off pos info, convert to numpy arr, apply min_maf
        # filter, convert PL to GL
        if self.win:
            self.w_pl_arr = np.array([x[1:] for x in self.win], dtype=np.int16)
            self.pl_convert_to_gl()
            self.gl_min_maf_filter()

        # empty: convert to empty numpy arr & apply target function
        else:
            self.w_gl_arr = np.empty((0,0))

        # call target function
        self.pcangsd()

        # restore sys.stdout
        sys.stdout = old_stdout


    def pca(self):
        '''
        Apply min_maf filter, drop rows/sites with missing call(s) OR mean
        impute, conduct PCA, but if (n_var < min_var_per_w) generate
        empty/dummy output instead.
        '''

        # get window mid for X value
        self.w_pos = int(self.w_start + self.w_size/2-1)

        # min_maf filter
        if self.min_maf:
            self.gt_min_maf_filter()

        # count variants
        n_var = self.w_gt_arr.shape[0]

        # count missing sites per sample
        n_miss_arr = np.isnan(self.w_gt_arr).sum(axis=0)

        # calculate non-missing sites per sample
        n_var_arr = n_var - n_miss_arr

        # calculate % het sites relative to n_var
        hetp_lst = list(np.sum(self.w_gt_arr == 1, axis=0)/n_var_arr)

        # mean impute
        if self.gt_mean_impute:
            self.gt_mean_imputation()
        
        # drop missing sites (& re-count)
        else:

            # drop missing sites
            self.gt_drop_missing_sites()

            # re-count count remaining variants
            n_var = self.w_gt_arr.shape[0]

        # if # variants passes specified threshold
        if n_var >= self.gt_min_var_per_w:
            # pca
            pca = allel.pca(
                self.w_gt_arr,
                n_components=self.n_pcs,
                copy=True,
                scaler='patterson',
                ploidy=2,
            )

            # compile to output
            out = {
                'pos': self.w_pos,
                'pc_a': pca[0][:, self.pc_a-1],
                'pc_b': pca[0][:, self.pc_b-1],
                f'pc_{self.pc_a}_ve': round(
                    pca[1].explained_variance_ratio_[self.pc_a-1]*100, 2
                ),
                f'pc_{self.pc_b}_ve': round(
                    pca[1].explained_variance_ratio_[self.pc_b-1]*100, 2
                ),
                'hetp': hetp_lst,
                'n_miss': n_miss_arr,
                'n_var': n_var
            }

        # else create empty output & print info
        else:
            empty_lst = [None] * self.w_gt_arr.shape[1]
            out = {
                'pos': self.w_pos,
                'pc_a': empty_lst,
                'pc_b': empty_lst,
                f'pc_{self.pc_a}_ve': None,
                f'pc_{self.pc_b}_ve': None,
                'hetp': empty_lst,
                'n_miss': empty_lst,
                'n_var': n_var
            }
            log.info('Skipped window'
                     f' {self.w_start}-{self.w_start + self.w_size-1} with'
                     f' {n_var} variants (theshold: {self.gt_min_var_per_w}'
                      ' variants)')

        # append output
        self.out_dct[self.w_idx] = out


    def pcangsd(self):
        '''
        Conduct PCAngsd PCA, but if (n_variants < min_var_per_w) generate
        empty/dummy output instead.
        '''

        # get window mid for X value
        self.w_pos = int(self.w_start + self.w_size/2-1)

        # count variants
        n_var = self.w_gl_arr.shape[0]

        # if # variants passes specified threshold
        if n_var >= self.gl_pl_min_var_per_w:

            # compute covariance matrix with PCAngsd
            # emPCA expects different numbers of args depending on the version
            try:
                cov_arr, _, _, _, _ = emPCA(
                    self.w_gl_arr,
                    self.gl_min_maf_arr,
                    self.n_pcs, 100, 1e-5,
                )
            except:
                cov_arr, _, _, _, _ = emPCA(
                    self.w_gl_arr,
                    self.gl_min_maf_arr,
                    self.n_pcs, 100, 1e-5,
                    self.n_threads
                )

            # eigendecomposition
            eigenval_arr, eigenvec_arr = np.linalg.eigh(cov_arr)

            # sort by eigenvalue
            idx = eigenval_arr.argsort()[::-1]
            eigenval_arr = eigenval_arr[idx]
            eigenvec_arr = eigenvec_arr[:,idx]

            # calculate % variance explained
            pct_exp_arr = [x/sum(eigenval_arr)*100 for x in eigenval_arr]

            out = {
                'pos': self.w_pos,
                'pc_a': eigenvec_arr[:, self.pc_a-1],
                'pc_b': eigenvec_arr[:, self.pc_b-1],
                f'pc_{self.pc_a}_ve': round(pct_exp_arr[self.pc_a-1], 2),
                f'pc_{self.pc_b}_ve': round(pct_exp_arr[self.pc_b-1], 2),
                'hetp': [None for x in eigenvec_arr[:, 0]],
                'n_miss': [None for x in eigenvec_arr[:, 0]],
                'n_var': n_var,
            }

        # else create empty output
        else:
            empty_lst = [None] * (self.w_gl_arr.shape[1]//2)
            out = {
                'pos': self.w_pos,
                'pc_a': empty_lst,
                'pc_b': empty_lst,
                f'pc_{self.pc_a}_ve': None,
                f'pc_{self.pc_b}_ve': None,
                'hetp': empty_lst,
                'n_miss': empty_lst,
                'n_var': n_var,
            }
            log.info('Skipped window'
                     f' {self.w_start}-{self.w_start + self.w_size-1} with'
                     f' {n_var} variants (theshold: {self.gl_pl_min_var_per_w}'
                      ' variants)')

        # append output
        self.out_dct[self.w_idx] = out


    def window_parser(self):
        '''
        Apply a target function to windows of variants (GT field) in an
        (optionally gzipped) VCF file.
        '''

        # calculate total number of windows
        self.n_windows = len(
            list(range(self.start, self.stop-self.w_size+2, self.w_step))
        )

        # determin compression
        read_func = (
            gzip.open if self.variant_file_path.endswith('.gz') else open
        )

        # for VCF: read first 10000 lines to check if multiallellic
        if self.file_fmt == 'VCF':
            with read_func(self.variant_file_path, 'rt') as check_gts:
                rows = [next(check_gts) for _ in range(10000)]
                rows = [x for x in rows if not x.startswith('#')]
                alts = [x.split('\t')[4] for x in rows]
                if any([len(x.split(',')) > 1 for x in alts]):
                    log.error_nl(
                        f'VARIANT_FILE: {self.variant_file_path}'
                            ' contains sites with more than 2 alleles'
                    )

        # open iput file
        with read_func(self.variant_file_path, 'rt') as variant_file:

            # fetch sample ids from different header types
            if self.file_fmt == 'VCF':
                for line in variant_file:
                    if line.startswith('#CHROM'):
                        var_file_sample_lst = line.strip().split('\t')[9:]
                        break
            if self.file_fmt == 'TSV':
                var_file_sample_lst = \
                    variant_file.readline().strip().split('\t')[2:]
            if self.file_fmt == 'BEAGLE':
                var_file_sample_lst = \
                    variant_file.readline().strip().split('\t')[3:]

            # use var_file_sample_lst (drop duplicates) if no samples specified
            if self.sample_lst is None:
                self.sample_lst = list(dict.fromkeys(var_file_sample_lst))

            # obtain index positions (returns first hit)
            sample_idx_lst = [
                var_file_sample_lst.index(x) for x in self.sample_lst
            ]

            # for PL: modify to account for 3 columns per sample
            if self.var_fmt in ['GL', 'PL'] and \
                self.file_fmt in ['TSV', 'BEAGLE']:

                # for GL, parse first two columns per sample, for PL all three
                if self.var_fmt == 'GL':
                    sample_idx_lst = [[i, i+1] for i in sample_idx_lst]
                if self.var_fmt == 'PL':
                    sample_idx_lst = [[i, i+1, i+2] for i in sample_idx_lst]
                sample_idx_lst = [x for i in sample_idx_lst for x in i]

        # initiate first window
        self.w_start = self.start
        self.w_stop = self.w_start + self.w_size-1
        self.w_idx = 0
        self.win = []

        # iterrate over windows
        with read_func(self.variant_file_path, 'rt') as variant_file:

            # GT
            if self.var_fmt == 'GT':

                if self.file_fmt == 'VCF':

                    for line in variant_file:
                        line = line.strip().split('\t')
                        q_chrom = line[0]
                        if q_chrom != self.chrom: continue
                        filter_field = line[6]
                        if filter_field != 'PASS' and self.vcf_pass_filter:
                            continue
                        self.pos = int(line[1])
                        gts = [
                            line[9:][idx].split(':')[0] \
                                for idx in sample_idx_lst
                        ]
                        gts = [self.gt_code_dct[x] for x in gts]
                        if self.skip_monomorphic \
                            and len(set(gts)-{np.nan}) == 1:
                            continue
                        if self.x_mode:
                            if self.parse_variant_x(
                                self.gt_process_win,
                                gts,
                            ): break
                        else:
                            if self.parse_variant(
                                self.gt_process_win,
                                gts,
                            ): break

                if self.file_fmt == 'TSV':
                    for line in variant_file:
                        line = line.strip().split('\t')
                        q_chrom = line[0]
                        if q_chrom != self.chrom: continue
                        self.pos = int(line[1])
                        gts = [line[2:][idx] for idx in sample_idx_lst]
                        if self.skip_monomorphic \
                            and len(set(gts)-{np.nan}) == 1:
                            continue
                        if self.x_mode:
                            if self.parse_variant_x(
                                self.gt_process_win,
                                gts,
                            ): break
                        else:
                            if self.parse_variant(
                                self.gt_process_win,
                                gts,
                            ): break

            # GL
            if self.var_fmt == 'GL':

                if self.file_fmt == 'VCF':
                    for line in variant_file:
                        line = line.strip().split('\t')
                        q_chrom = line[0]
                        if q_chrom != self.chrom: continue
                        format_field = line[8].split(':')
                        if not self.var_fmt in format_field: continue
                        filter_field = line[6]
                        if filter_field != 'PASS' and self.vcf_pass_filter:
                            continue
                        self.pos = int(line[1])
                        gt_fields = line[9:]
                        # get each sample's GL field as a list
                        gls = [
                            gt_fields[idx].split(':')[
                                format_field.index(self.var_fmt)
                            ].split(',') \
                                for idx in sample_idx_lst
                        ]
                        # only keep biallellic GLs (i.e. AA, AB, BB)
                        gls = [x for gl in gls if len(gl) == 3 for x in gl]
                        # GATK encodes missing data as '.' --> drop lines
                        # where length of GLs != 3* n_samples
                        gls = \
                            [] if (len(gls)) != len(sample_idx_lst)* 3 else gls
                        if gls == []: continue
                        # delete 3rd field for each GL (expected by PCAngsd)
                        gls = np.delete(gls, np.s_[2::3], axis=1)
                        if self.x_mode:
                            if self.parse_variant_x(
                                self.gl_process_win,
                                gts,
                            ): break
                        else:
                            if self.parse_variant(
                                self.gl_process_win,
                                gts,
                            ): break

                if self.file_fmt == 'TSV':
                    for line in variant_file:
                        line = line.strip().split('\t')
                        q_chrom = line[0]
                        if q_chrom != self.chrom: continue
                        self.pos = int(line[1])
                        gls = [line[2:][idx] for idx in sample_idx_lst]
                        if self.x_mode:
                            if self.parse_variant_x(
                                self.gl_process_win,
                                gts,
                            ): break
                        else:
                            if self.parse_variant(
                                self.gl_process_win,
                                gts,
                            ): break

                if self.file_fmt == 'BEAGLE':
                    for line in variant_file:
                        line = line.strip().split('\t')
                        q_chrom = line[0].rsplit('_', 1)[0]
                        if q_chrom != self.chrom: continue
                        self.pos = int(line[0].rsplit('_', 1)[1])
                        gls = [line[3:][idx] for idx in sample_idx_lst]
                        if self.x_mode:
                            if self.parse_variant_x(
                                self.gl_process_win,
                                gts,
                            ): break
                        else:
                            if self.parse_variant(
                                self.gl_process_win,
                                gts,
                            ): break


            # PL
            if self.var_fmt == 'PL':

                if self.file_fmt == 'VCF':
                    for line in variant_file:
                        line = line.strip().split('\t')
                        q_chrom = line[0]
                        if q_chrom != self.chrom: continue
                        format_field = line[8].split(':')
                        if not self.var_fmt in format_field: continue
                        self.pos = int(line[1])
                        filter_field = line[6]
                        if filter_field != 'PASS' and self.vcf_pass_filter:
                            continue
                        gt_fields = line[9:]
                        # get each sample's PL field as a list
                        pls = [
                            gt_fields[idx].split(':')[
                                format_field.index(self.var_fmt)
                            ].split(',') \
                                for idx in sample_idx_lst
                        ]
                        # only keep biallellic PLs (i.e. AA, AB, BB)
                        pls = [x for pl in pls if len(pl) == 3 for x in pl]
                        # GATK encodes missing PL data as '.' --> drop lines
                        # where length of PLs != 3* n_samples
                        pls = \
                            [] if (len(pls)) != len(sample_idx_lst)*3 else pls
                        if pls == []: continue
                        if self.x_mode:
                            if self.parse_variant_x(
                                self.pl_process_win,
                                gts,
                            ): break
                        else:
                            if self.parse_variant(
                                self.pl_process_win,
                                gts,
                            ): break

                if self.file_fmt == 'TSV':
                    for line in variant_file:
                        line = line.strip().split('\t')
                        q_chrom = line[0]
                        if q_chrom != self.chrom: continue
                        self.pos = int(line[1])
                        pls = [line[2:][idx] for idx in sample_idx_lst]
                        if self.x_mode:
                            if self.parse_variant_x(
                                self.pl_process_win,
                                gts,
                            ): break
                        else:
                            if self.parse_variant(
                                self.pl_process_win,
                                gts,
                            ): break

        # check if any windows were processed
        if len(self.out_dct) == 0:
            log.error('No windows found. Please check if chromosome name was'
            ' specified correctly')
            log.newline()

        # print exit message
        log.newline()
        log.info('Processed all windows')
