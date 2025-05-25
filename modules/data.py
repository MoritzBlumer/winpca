'''
Central class to store, read, write and modify PCA and associated data.
'''

## IMPORT PACKAGES
import os
import pandas as pd

## MODULES
from modules.log import Log

## INSTANTIATE LOGGER
log = Log()


## CLASSES

class WPCAData:
    '''
    Central container for PCA and associated data. Compiles object from
    existing output files (via prefix) or directly from w_pca instance. Data
    can be modified through indluded modifier methods.
    '''

    def __init__(self, prefix, pc_a, pc_b, w_pca_obj=None):

        # input variables
        self.prefix = prefix
        self.w_pca_obj = w_pca_obj
        self.pc_a = pc_a
        self.pc_b = pc_b

        # instance variables
        self.pc_a_df = None
        self.pc_b_df = None
        self.hetp_df = None
        self.miss_df = None
        self.stat_df = None

        # set file suffixes
        self.suffix_dct = {
            'pc_a': f'pc_{self.pc_a}.tsv.gz',
            'pc_b': f'pc_{self.pc_b}.tsv.gz',
            'hetp': 'hetp.tsv.gz',
            'miss': 'miss.tsv.gz',
            'stat': 'stat.tsv.gz',
        }

        # compile object from w_pca instance
        if w_pca_obj:
            self.from_w_pca()

        # compile object from files
        else:
            self.from_files()


    def from_w_pca(self):
        '''
        Compile from out_dct of w_pca instance.
        '''

        # convert to df
        out_df = pd.DataFrame.from_dict(self.w_pca_obj.out_dct, orient='index')

        # pc_a
        self.pc_a_df = pd.DataFrame(
            data=out_df['pc_a'].to_list(),
            index=out_df['pos'],
            columns=self.w_pca_obj.sample_lst,
        )

        # pc_b
        self.pc_b_df = pd.DataFrame(
            data=out_df['pc_b'].to_list(),
            index=out_df['pos'],
            columns=self.w_pca_obj.sample_lst,
        )

        # hetp
        self.hetp_df = pd.DataFrame(
            data=out_df['hetp'].to_list(),
            index=out_df['pos'],
            columns=self.w_pca_obj.sample_lst,
        )

        # miss
        self.miss_df = pd.DataFrame(
            data=out_df['n_miss'].to_list(),
            index=out_df['pos'],
            columns=self.w_pca_obj.sample_lst,
            dtype='Int64',
        )

        # stats: window index, pc_a/pc_b variance explained, number of variants
        self.stat_df = out_df[
            ['pos', 'n_var', f'pc_{self.pc_a}_ve', f'pc_{self.pc_b}_ve', ]
        ].copy()
        self.stat_df = self.stat_df.reset_index()
        self.stat_df = self.stat_df.rename(columns={'index': 'w_idx'})
        self.stat_df = self.stat_df.set_index('pos')


    def from_files(self):
        '''
        Compile from out_dct of WPCA instance.
        '''

        # print info
        log.newline()
        log.info(f'Reading data from prefix "{self.prefix}*"')

        # flush existing instance variables
        self.pc_a_df, self.pc_b_df, \
            self.hetp_df, self.miss_df, self.stat_df = \
            None, None, None, None, None

        # pc_a
        self.pc_a_df = pd.read_csv(
            f'{self.prefix}.{self.suffix_dct["pc_a"]}',
            sep='\t',
            index_col='pos',
        )

        # pc_b
        self.pc_b_df = pd.read_csv(
            f'{self.prefix}.{self.suffix_dct["pc_b"]}',
            sep='\t',
            index_col='pos',
        )

        # hetp
        self.hetp_df = pd.read_csv(
            f'{self.prefix}.{self.suffix_dct["hetp"]}',
            sep='\t',
            index_col='pos',
        )

        # # miss
        # self.miss_df = pd.read_csv(
        #     f'{self.prefix}.{self.suffix_dct["miss"]}',
        #     sep='\t',
        #     index_col='pos',
        #     dtype='Int64'
        # )

        # stat
        self.stat_df = pd.read_csv(
            f'{self.prefix}.{self.suffix_dct["stat"]}',
            sep='\t',
            index_col='pos',
        )


    def to_files(self):
        '''
        Write output (pc_a_df, pc_b_df, hetp_df, miss_df, stat_df) to files.
        '''

        # create output directory if prefix contains '/'
        if '/' in self.prefix:
            if not os.path.exists(
                f'{"/".join(self.prefix.split('/')[0:-1])}/'
            ):
                os.makedirs(f'"/".join(self.prefix.split('/')[0:-1])/')

        # pc_a
        self.pc_a_df = self.pc_a_df.round(3)
        self.pc_a_df.to_csv(
            f'{self.prefix}.{self.suffix_dct["pc_a"]}',
            sep='\t', index_label='pos', na_rep='NA')

        # pc_b
        self.pc_b_df = self.pc_b_df.round(3)
        self.pc_b_df.to_csv(
            f'{self.prefix}.{self.suffix_dct["pc_b"]}',
            sep='\t', index_label='pos', na_rep='NA')

        # hetp
        self.hetp_df = self.hetp_df.round(3)
        self.hetp_df.to_csv(
            f'{self.prefix}.{self.suffix_dct["hetp"]}',
            sep='\t', index_label='pos', na_rep='NA')

        # # miss
        # self.miss_df = self.miss_df.round(3)
        # self.miss_df.to_csv(
        #     f'{self.prefix}.{self.suffix_dct["miss"]}',
        #     sep='\t', index_label='pos', na_rep='NA')

        # stat
        self.stat_df = self.stat_df.round(3)
        self.stat_df.to_csv(
            f'{self.prefix}.{self.suffix_dct["stat"]}',
            sep='\t', index_label='pos', na_rep='NA')


    def modify_data(self, attr_name, mod_func, *args):
        '''
        Modify an attribute of Data instance with a supplied function.
        '''
        # get attribute
        attr = getattr(self, attr_name)

        # modify
        attr_mod = mod_func(attr, *args)

        # update attribute
        setattr(self, attr_name, attr_mod)
