'''
Central window parsers and PCA functions.
'''

## IMPORT PACKAGES
import sys
import gzip
import numpy as np
import allel


## IMPORT MODULES
from modules import config


## CLASSES

class GT_WPCA:

    '''
    Parse hard-called genotypes (GT) from TSV or VCF and apply a function.
    '''

    def __init__(self, variant_file_path, sample_lst,
                 chrom, start, stop,
                 w_size=config.W_SIZE,
                 w_step=config.W_STEP,
                 min_var_per_w=config.MIN_VAR_PER_W,
                 skip_monomorphic=config.SKIP_MONOMORPHIC,
                 min_maf=config.MIN_MAF,):

        # variant file/region
        self.variant_file_path = variant_file_path
        self.sample_lst = sample_lst
        self.chrom = chrom
        self.start = start
        self.stop = stop

        # parameters (defaults from config)
        self.w_size = w_size
        self.w_step = w_step
        self.min_var_per_w = min_var_per_w
        self.min_maf = min_maf
        self.skip_monomorphic = skip_monomorphic

        # transient variables
        self.n_windows = None
        self.w_start = None
        self.w_stop = None
        self.pos = None
        self.w_idx = None
        self.win = None
        self.w_gt_arr = None

        # results
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


#    @njit()                                                                        ### BRING BACK
    def gt_min_maf_filter(self):
        '''
        Drop SNPs with minor allele frequency below specified value.
        '''

        # allele count
        n_alleles = 2 * self.w_gt_arr.shape[1]

        # calculate allel frequencies and multiple with -1 if AF > 0.05 (because
        # input data may not be polarized by major/minor allel)
        afs = np.sum(self.w_gt_arr, axis=1) / n_alleles
        afs[afs > 0.5] = 1 - afs[afs > 0.5]

        # keep only sites where AF >= min_maf
        self.w_gt_arr = self.w_gt_arr[afs >= self.min_maf]


    def gt_process_win(self):
        '''
        Remove POS info, convert to numpy array, apply minMAF filter and return
        empty array if there are no variants. Call target function:
        func(w_gt_arr, w_start, w_size).
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
        self.func()


    def func(self): # incorporate skip_monomorphic here
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


    def win_gt_file(self):
        '''
        Apply a target function to windows of variants in an (optionally
        gzipped) genotype TSV file.
        '''

        # calculate total number of windows
        self.n_windows = len(
            list(range(self.start, self.stop-self.w_size+2, self.w_step))
        )

        # open uncompressed or gzip-compressed input file
        read_func = gzip.open if self.variant_file_path.endswith('.gz') else open
        with read_func(self.variant_file_path, 'rt') as variant_file:

            # fetch sample ids from header and obtain index positions
            variant_file_sample_lst = \
                variant_file.readline().strip().split('\t')[2:]
            sample_idx_lst = [
                variant_file_sample_lst.index(x) for x in self.sample_lst
            ]

            # initiate first window
            self.w_start = self.start
            self.w_stop = self.w_start + self.w_size-1
            self.w_idx = 0
            self.win = []

            # traverse input file
            for line in variant_file:
                line = line.strip().split('\t')
                q_chrom = line[0]
                pos = int(line[1])

                # skip other than the specified chromosome
                if q_chrom != self.chrom:
                    continue

                # fetch genotypes
                gts = [line[2:][idx] for idx in sample_idx_lst]

                # skip monomorphic sites if specified
                if self.skip_monomorphic and len(set(gts)) == 1:
                    continue                                                    ### MODIFY (?)

                # case: pos exceeds current window
                while self.w_stop < pos:

                    # apply min_maf filter if specified and if window contains
                    # variants: apply function
                    if self.win:
                        self.gt_process_win()
                    if self.stop < self.w_stop:
                        break
                    self.init_win()

                # append pos (and genotypes) to current window if larger than
                # window start
                if pos > self.w_start:
                    self.win.append([pos] + gts)

                # if end of window is reached: apply min_maf filter, function
                # & initialize
                if self.w_stop <= pos:
                    self.gt_process_win()
                    if self.stop < self.w_stop:
                        break
                    self.init_win()

        # print info
        print(
            '\n[INFO] Processed all windows',
            file=sys.stderr, flush=True,
        )


    def win_vcf_gt(self):
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

            for line in variant_file:
                if line.startswith('#CHROM'):

                    # fetch sample ids from header and obtain index positions
                    variant_file_sample_lst = line.strip().split('\t')[9:]
                    sample_idx_lst = [
                        variant_file_sample_lst.index(x) for x in self.sample_lst
                    ]
                    break

        # initiate first window
        self.w_start = self.start
        self.w_stop = self.w_start + self.w_size-1
        self.w_idx = 0
        self.win = []

        with read_func(self.variant_file_path, 'rt') as variant_file:

            # traverse input file
            for line in variant_file:
                if line.startswith('#'):
                    continue

                line = line.strip().split('\t')
                q_chrom = line[0]
                pos = int(line[1])
                filter_field = line[6]

                # skip other than the specified chromosome
                if q_chrom != self.chrom:
                    continue

                # keep only PASS sites
                if filter_field != 'PASS':
                    continue

                # fetch genotypes
                gts = [line[9:][idx].split(':')[0] for idx in sample_idx_lst]
                gts = [self.gt_code_dct[x] for x in gts]

                # skip monomorphic sites if specified
                if self.skip_monomorphic and len(set(gts)) == 1:
                    continue                                                    ### MODIFY (?)

                # case: pos exceeds current window
                while self.w_stop < pos:

                    # apply min_maf filter if specified and if window contains
                    # variants: apply function
                    if self.win:
                        self.gt_process_win()
                    if self.stop < self.w_stop:
                        break
                    self.init_win()

                # append pos (and genotypes) to current window
                if pos > self.w_start:
                    self.win.append([pos] + gts)

                # if end of window is reached: apply min_maf filter, function \
                # & initialize
                if self.w_stop <= pos:
                    self.gt_process_win()
                    if self.stop < self.w_stop:
                        break
                    self.init_win()

        # print exit message
        print(
            '\n[INFO] Processed all windows',
            file=sys.stderr, flush=True,
        )
