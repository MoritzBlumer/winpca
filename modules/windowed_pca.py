'''
Window parser and GT PCA functions.
'''

from modules import config                                                      # DELETE

## SET NUMBER OF THREADS

# load PCAngsd if working on GL/PL                                              # IMPROVE IMPORT SECTION
import os
#if config.VAR_FMT == 'GT':
import allel
#else:
os.environ["OMP_NUM_THREADS"] = str(config.N_THREADS)
os.environ["OPENBLAS_NUM_THREADS"] = str(config.N_THREADS)
os.environ["MKL_NUM_THREADS"] = str(config.N_THREADS)


## IMPORT PACKAGES
import sys, os
import gzip
import numpy as np
from numba import njit, prange, set_num_threads                                 # ???


## IMPORT MODULES
from pcangsd.shared import emMAF
from pcangsd.reader_cy import filterArrays
from pcangsd.covariance import emPCA


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
                 chrom, start, stop,
                 w_size, w_step, 
                 min_var_per_w,
                 skip_monomorphic,
                 min_maf,
                 n_threads,
                 ):

        # variant file/region
        self.variant_file_path = variant_file_path
        self.file_fmt = file_fmt
        self.var_fmt = var_fmt
        self.sample_lst = sample_lst
        self.chrom = chrom
        self.start = start
        self.stop = stop

        # parameters (defaults from config)
        self.w_size = w_size
        self.w_step = w_step
        self.min_var_per_w = min_var_per_w
        self.skip_monomorphic = skip_monomorphic
        self.min_maf = min_maf
        self.n_threads = n_threads

        # transient variables
        self.n_windows = None
        self.w_start = None
        self.w_stop = None
        self.pos = None
        self.w_idx = None
        self.win = None
        self.w_gt_arr = None
        self.w_gl_arr = None
        self.w_pl_arr = None
        self.gl_min_maf_arr = None

        # instantiate results dict
        self.out_dct = {
        }

        # genotype encoding
        self.gt_code_dct = {
            '0/0': 0,
            '0|0': 0,
            '0/1': 1,
            '0|1': 1,
            '1/0': 1,
            '1|0': 1,
            '1/1': 2,
            '1|1': 2,
            './.': -1,
            '.|.': -1,
            '0/.': -1,
            '0|.': -1,
            './0': -1,
            '.|0': -1,
            '1/.': -1,
            '1|.': -1,
            './1': -1,
            '.|1': -1
        }


    def init_win(self):
        '''
        Initialize new window by shifting one w_step and dropping obsolete
        variants from previous window.
        '''

        self.w_start = self.w_start + self.w_step
        self.w_stop = self.w_start + self.w_size-1
        self.w_idx += 1
        self.win = [x for x in self.win if x[0] >= self.w_start]

        print(
            '[INFO] Processed ' + '' + str(self.w_idx) + ' of ' +
            str(self.n_windows) + ' windows',
                file=sys.stderr, flush=True,
            )


#    @njit()                                                                    ### BRING BACK
    def gt_min_maf_filter(self):
        '''
        Drop SNPs with minor allele frequency below specified value.
        '''

        # allele count
        n_alleles = 2 * self.w_gt_arr.shape[1]

        # calculate allel frequencies and multiple with -1 if AF > 0.5 (because
        # input data may not be polarized by major/minor allel)
        afs = np.sum(self.w_gt_arr, axis=1) / n_alleles
        afs[afs > 0.5] = 1 - afs[afs > 0.5]

        # keep only sites where AF >= min_maf
        self.w_gt_arr = self.w_gt_arr[afs >= self.min_maf]


    def gt_process_win(self):
        '''
        Remove POS info, convert to numpy array, apply min_maf filter, call 
        target function or return empty array if there are no variants.
        '''

        # non-empty: trim off pos info, convert to numpy arr, apply min_maf
        # filter
        if self.win:
            self.w_gt_arr = np.array([x[1:] for x in self.win], dtype=np.int8)
            if self.min_maf:
                self.gt_min_maf_filter()

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
        self.gl_min_maf_arr = emMAF(self.w_gl_arr, 200, 1e-4, self.n_threads)

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


    # @njit(parallel=True)                                                      # BRING BACK
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
        w_pl_arr_3d = self.w_pl_arr.reshape(n_rows, -1, 3).transpose(1, 0, 2)

        # unphred
        w_pl_arr_3d = np.power(10, -w_pl_arr_3d/10)

        # preallocate output w_gl_arr (PCAngsd expects only first two GL values
        # for memory efficiency)
        self.w_gl_arr = np.zeros((n_rows, n_cols*2//3), dtype=np.float32)

        # normalize all GLs partitioned by sample in parallel using numba prange
        for idx in prange(0, w_pl_arr_3d.shape[0]):
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


    def pca(self):                                                              # incorporate skip_monomorphic here?
        '''
        Conduct PCA, but if (n_var < min_var_per_w) generate empty/dummy
        output instead.
        '''

        # get window mid for X value
        pos = int(self.w_start + self.w_size/2-1)

        # count variants
        n_var = self.w_gt_arr.shape[0]

        # count missing sites per sample
        n_miss_lst = list(np.sum(self.w_gt_arr == -1, axis=0))

        # if # variants passes specified threshold
        if n_var >= self.min_var_per_w:

            # pca
            pca = allel.pca(
                self.w_gt_arr,
                n_components=2,
                copy=True,
                scaler='patterson',
                ploidy=2,
            )

            # calculate % het sites relative to n_var
            hetp_lst = list(np.sum(self.w_gt_arr == 1, axis=0)/n_var)

            # compile to output
            out = {
                'pos': pos,
                'pc_1': pca[0][:, 0],
                'pc_2': pca[0][:, 1],
                'pc_1_ve': round(pca[1].explained_variance_ratio_[0]*100, 2),
                'pc_2_ve': round(pca[1].explained_variance_ratio_[1]*100, 2),
                'hetp': hetp_lst,
                'n_miss': n_miss_lst,
                'n_var': n_var
            }

        # else create empty output & print info
        else:
            print(
                '[INFO] Skipped window ' + str(self.w_start) + '-'
                + str(self.w_start + self.w_size-1) + ' with ' + str(n_var) +
                ' variants (threshold is ' + str(self.min_var_per_w) +
                ' variants per window)',
                file=sys.stderr, flush=True,
            )
            empty_lst = [None] * self.w_gt_arr.shape[1]
            out = {
                'pos': pos,
                'pc_1': empty_lst,
                'pc_2': empty_lst,
                'pc_1_ve': None,
                'pc_2_ve': None,
                'hetp': empty_lst,
                'n_miss': empty_lst,
                'n_var': n_var
            }

        # append output
        self.out_dct[self.w_idx] = out


    def pcangsd(self):
        '''
        Conduct PCAngsd PCA, but if (n_variants < min_var_per_w) generate
        empty/dummy output instead.
        '''

        # get window mid for X value
        pos = int(self.w_start + self.w_size/2-1)

        # count variants
        n_var = self.w_gl_arr.shape[0]

        # if # variants passes specified threshold  
        if n_var >= self.min_var_per_w:

            # compute covariance matrix with PCAngsd
            cov_arr, _, _, _, _ = emPCA(
                self.w_gl_arr,
                self.gl_min_maf_arr, 
                0, 100, 1e-5,
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
                'pos': pos,
                'pc_1': eigenvec_arr[:, 0],
                'pc_2': eigenvec_arr[:, 1],
                'pc_1_ve': round(pct_exp_arr[0], 2),
                'pc_2_ve': round(pct_exp_arr[1], 2),
                'hetp': [None for x in eigenvec_arr[:, 0]],
                'n_miss': [None for x in eigenvec_arr[:, 0]],
                'n_var': n_var,
            }

        # else create empty output
        else:
            print(
                '[INFO] Skipped window ' + str(self.w_start) + '-' + str(self.w_start + self.w_size-1) + ' with ' +
                str(n_var) + ' variants (threshold is ' + str(self.min_var_per_w) +
                ' variants per window)',
                file=sys.stderr, flush=True,
            )

            empty_lst = [None] * (self.w_gl_arr.shape[1]//2)
            out = {
                'pos': pos,
                'pc_1': empty_lst,
                'pc_2': empty_lst,
                'pc_1_ve': None,
                'pc_2_ve': None,
                'hetp': empty_lst,
                'n_miss': empty_lst,
                'n_var': n_var,
            }

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

        # open iput file
        read_func = gzip.open if self.variant_file_path.endswith('.gz') else open
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
                
            #print(f'sample_lst: {self.sample_lst}', flush=True)                # DELETE
            #print(f'var_file_sample_lst: {var_file_sample_lst}', flush=True)   # DELETE
            # obtain index positions
            sample_idx_lst = [
                var_file_sample_lst.index(x) for x in self.sample_lst
            ]
            #print(f'sample_idx_lst: {sample_idx_lst}', flush=True)             # DELETE
            # for PL: modify to account for 3 columns per sample
            if self.var_fmt in ['GL', 'PL'] and \
                self.file_fmt in ['TSV', 'BEAGLE']:

                # for GL, parse first two columns per sample, for PL all three
                if self.var_fmt == 'GL':
                    sample_idx_lst = [[i, i+1] for i in sample_idx_lst]
                if self.var_fmt == 'PL':
                    sample_idx_lst = [[i, i+1, i+2] for i in sample_idx_lst]
                sample_idx_lst = [x for i in sample_idx_lst for x in i]
                #print(f'sample_idx_lst: {sample_idx_lst}', flush=True)         # DELETE

                # remove duplicates from var_file_sample_lst (preserve order)
                var_file_sample_lst = list(dict.fromkeys(var_file_sample_lst))
                #print(f'var_file_sample_lst: {var_file_sample_lst}', flush=True)# DELETE

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
                        if filter_field != 'PASS': continue
                        pos = int(line[1])
                        gts = [line[9:][idx].split(':')[0] for idx in sample_idx_lst]
                        gts = [self.gt_code_dct[x] for x in gts]
                        if self.skip_monomorphic and len(set(gts)) == 1: continue
                        while self.w_stop < pos:
                            if self.win: self.gt_process_win()
                            if self.stop < self.w_stop: break
                            self.init_win()
                        if pos > self.w_start: self.win.append([pos] + gts)
                        if self.w_stop <= pos:
                            self.gt_process_win() 
                            if self.stop < self.w_stop: break
                            self.init_win()

                if self.file_fmt == 'TSV':
                    for line in variant_file:
                        line = line.strip().split('\t')
                        q_chrom = line[0]
                        if q_chrom != self.chrom: continue
                        pos = int(line[1])
                        gts = [line[2:][idx] for idx in sample_idx_lst]
                        if self.skip_monomorphic and len(set(gts)) == 1: continue
                        while self.w_stop < pos:
                            if self.win: self.gt_process_win()
                            if self.stop < self.w_stop: break
                            self.init_win()
                        if pos > self.w_start: self.win.append([pos] + gts)
                        if self.w_stop <= pos:
                            self.gt_process_win() 
                            if self.stop < self.w_stop: break
                            self.init_win()

            # GL
            if self.var_fmt == 'GL':                                            ### ADD
                if self.file_fmt == 'VCF':
                    pass

                if self.file_fmt == 'VCF':
                    for line in variant_file:
                        line = line.strip().split('\t')
                        q_chrom = line[0]
                        if q_chrom != self.chrom: continue
                        format = line[8].split(':') 
                        if not self.var_fmt in format: continue 
                        pos = int(line[1])
                        gt_fields = line[9:]
                        # get each sample's GL field as a list #### REFINE for 2 cols
                        gls = [
                            gt_fields[idx].split(':')[
                                format.index(self.var_fmt)
                            ].split(',') \
                                for idx in sample_idx_lst
                        ]
                        # only keep biallellic GLs (i.e. AA, AB, BB)
                        gls = [x for gl in gls if len(gl) == 3 for x in gl]
                        # delete 3rd field for each GL (expected by PCAngsd)
                        gls = np.delete(gls, np.s_[2::3], axis=1)
                        # GATK encodes missing data as '.' --> drop lines
                        # where length of PLs != 3* n_samples
                        gls = [] if (len(gls)) != len(sample_idx_lst)*3 else gls            
                        if pls == []: continue

                        while self.w_stop < pos: 
                            if self.win: self.pl_process_win()
                            if self.stop < self.w_stop: break
                            self.init_win() 
                        if pos > self.w_start: self.win.append([pos] + pls) 
                        if self.w_stop <= pos:
                            self.pl_process_win() 
                            if self.stop < self.w_stop: break
                            self.init_win()

                if self.file_fmt == 'TSV':
                    for line in variant_file:
                        line = line.strip().split('\t')
                        q_chrom = line[0]
                        if q_chrom != self.chrom: continue
                        pos = int(line[1])
                        gls = [line[2:][idx] for idx in sample_idx_lst]
                        while self.w_stop < pos:
                            if self.win: self.gl_process_win()
                            if self.stop < self.w_stop: break 
                            self.init_win()
                        if pos > self.w_start: self.win.append([pos] + gls) 
                        if self.w_stop <= pos: 
                            self.gl_process_win() 
                            if self.stop < self.w_stop: break

                if self.file_fmt == 'BEAGLE':
                    for line in variant_file:
                        line = line.strip().split('\t')
                        q_chrom = line[0].rsplit('_', 1)[0]
                        if q_chrom != self.chrom: continue
                        pos = int(line[0].rsplit('_', 1)[1])
                        gls = [line[2:][idx] for idx in sample_idx_lst]
                        while self.w_stop < pos:
                            if self.win: self.gl_process_win()
                            if self.stop < self.w_stop: break 
                            self.init_win()
                        if pos > self.w_start: self.win.append([pos] + gls) 
                        if self.w_stop <= pos: 
                            self.gl_process_win() 
                            if self.stop < self.w_stop: break

            # PL
            if self.var_fmt == 'PL':

                if self.file_fmt == 'VCF':
                    for line in variant_file:
                        line = line.strip().split('\t')
                        q_chrom = line[0]
                        if q_chrom != self.chrom: continue
                        format = line[8].split(':') 
                        if not self.var_fmt in format: continue 
                        pos = int(line[1])
                        gt_fields = line[9:]
                        # get each sample's PL field as a list
                        pls = [
                            gt_fields[idx].split(':')[
                                format.index(self.var_fmt)
                            ].split(',') \
                                for idx in sample_idx_lst
                        ]
                        # only keep biallellic PLs (i.e. AA, AB, BB)
                        pls = [x for pl in pls if len(pl) == 3 for x in pl]
                        # GATK encodes missing PL data as '.' --> drop lines
                        # where length of PLs != 3* n_samples
                        pls = [] if (len(pls)) != len(sample_idx_lst)*3 else pls               
                        if pls == []: continue
                        while self.w_stop < pos: 
                            if self.win: self.pl_process_win()
                            if self.stop < self.w_stop: break
                            self.init_win() 
                        if pos > self.w_start: self.win.append([pos] + pls) 
                        if self.w_stop <= pos:
                            self.pl_process_win() 
                            if self.stop < self.w_stop: break
                            self.init_win()

                if self.file_fmt == 'TSV':
                    for line in variant_file:
                        line = line.strip().split('\t')
                        q_chrom = line[0]
                        if q_chrom != self.chrom: continue
                        pos = int(line[1])
                        pls = [line[2:][idx] for idx in sample_idx_lst]
                        # GATK encodes missing PL data as '.' --> drop lines
                        # where length of PLs != 3* n_samples
                        pls = [] if (len(pls)) != len(sample_idx_lst) else pls  
                        if pls == []: continue
                        while self.w_stop < pos: 
                            if self.win: self.pl_process_win()
                            if self.stop < self.w_stop: break
                            self.init_win() 
                        if pos > self.w_start: self.win.append([pos] + pls) 
                        if self.w_stop <= pos:
                            self.pl_process_win() 
                            if self.stop < self.w_stop: break
                            self.init_win()

        # print exit message
        print(
            '\n[INFO] Processed all windows',
            file=sys.stderr, flush=True,
        )

#            # traverse input file
#            for line in variant_file:                                           # ALL
#
#                # skip header
#                if line.startswith('#'): continue                               # DELETE             (only applies to VCF-GT, VCF-PL) ### ACTUALLY DELETE, SINCE POINTER ISN'T RESET AFTER PARSING THE HEADER
#
#
#                # FETCH LINE
#                line = line.strip().split('\t')                                 # ALL
#
#
#                # PROCESS LINE
#
#                q_chrom = line[0]                                               # ALL_BUT_BEAGLE
#                # q_chrom = line[0].rsplit('_', 1)[0]                           #     BEAGLE
#
#                # skip other than the specified chromosome
#                if q_chrom != self.chrom: continue                              # ALL (VCF-PL has this after q_chrom = line[0] to enhance performance --> makes sense for all)       
#
#                # keep only PASS sites
#                filter_field = line[6]                                          # VCF-GT
#                if filter_field != 'PASS': continue                             # VCF-GT
#
#                pos = int(line[1])                                              # ALL_BUT_BEAGLE
#                # pos = int(line[0].rsplit('_', 1)[1])                          #     BEAGLE
#
#                # format = line[8].split(':')                                   # VCF-PL (makes sense for VCF-GL as well)
#                # if not 'PL' in format: continue                               # VCF-PL (makes sense for VCF-GL as well)
#
#
#                # FETCH VARIANTS
#                gts = [line[9:][idx].split(':')[0] for idx in sample_idx_lst]   # VCF-GT
#                gts = [self.gt_code_dct[x] for x in gts]                        # VCF-GT
#                # gts = [line[2:][idx] for idx in sample_idx_lst]               # TSV-GT
#                # pls = [line[2:][idx] for idx in sample_idx_lst]               # TSV-PL
#                # pls = pl_fetch_target_pls(                                    # VCF-PL
#                #       line[9:], format.index('PL'), sample_idx_lst)           # VCF-PL
#                # gls = [line[2:][idx] for idx in sample_idx_lst]               # TSV-GL
#                # gls = [line[3:][idx] for idx in sample_idx_lst]               # BEAGLE
#
#
#
#                # PROCESS VARIANTS
#
#                # skip monomorphic sites if specified                           
#                if self.skip_monomorphic and len(set(gts)) == 1: continue       # VCF-GT TSV-GT                                                   ### MODIFY (?)
#
#                # # skip non-biallellic sites                                   # VCF-PL
#                # if pls == []: continue                                        # VCF-PL 
#
#
#                # PROCESS WINDOWS
#
#
#
#                # case: pos exceeds current window
#                while self.w_stop < pos:                                        # ALL
#
#                    # apply min_maf filter if specified and if window contains
#                    # variants: apply function
#                    if self.win: self.gt_process_win()                          # VCF-GT TSV-GT
#                  # if win: pl_process_win(win, w_start, w_size, min_maf, func) #               VCF-PL TSV-PL
#                  # if win: gl_process_win(win, w_start, w_size, min_maf, func) #                             TSV-GL BEAGLE    
#                    if self.stop < self.w_stop: break                           # ALL
#                    self.init_win()                                             # ALL
#
#                # append pos (and genotypes) to current window
#                if pos > self.w_start: self.win.append([pos] + gts)             # VCF-GT TSV-GT              
#                # if pos > w_start: win.append([pos] + pls)                     #               VCF-PL TSV-PL
#                # if pos > w_start: win.append([pos] + gls)                     #                             TSV-GL BEAGLE
#
#                # if end of window is reached: apply min_maf filter, function \
#                # & initialize
#                if self.w_stop <= pos:                                          # ALL
#                    self.gt_process_win()                                       # VCF-GT TSV-GT 
#                    # pl_process_win(win, w_start, w_size, min_maf, func)                       VCF-PL TSV-PL
#                    # gl_process_win(win, w_start, w_size, min_maf, func)       #                             TSV-GL BEAGLE
#                    if self.stop < self.w_stop: break                           # ALL
#                    self.init_win()                                             # ALL


