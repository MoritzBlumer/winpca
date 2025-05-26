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
    Central container for PCA and associated data, compiles object from
    existing output files (via prefix) or directly from w_pca instance; data
    can be modified through indluded modifier methods
    '''

    def __init__(self, prefix, pcs, w_pca_obj=None):

        # input variables
        self.prefix = prefix
        self.w_pca_obj = w_pca_obj
        self.pcs = pcs

        # compile object from w_pca instance
        if w_pca_obj:
            self.from_w_pca()

        # compile object from files
        else:
            self.from_files()


    def from_w_pca(self):
        '''
        Compile from out_dct of w_pca instance
        '''

        # convert to df
        out_df = pd.DataFrame.from_dict(self.w_pca_obj.out_dct, orient='index')

        # pc_dfs
        for pc in self.pcs:
            setattr(
                self,
                f'pc_{pc}_df',
                pd.DataFrame(
                    data=out_df[f'pc_{pc}'].to_list(),
                    index=out_df['pos'],
                    columns=self.w_pca_obj.sample_lst,
                ).round(3)
            )

        # hetp
        self.hetp_df = pd.DataFrame(
            data=out_df['hetp'].to_list(),
            index=out_df['pos'],
            columns=self.w_pca_obj.sample_lst,
        ).round(3)

        # miss
        # self.miss_df = pd.DataFrame(
        #     data=out_df['n_miss'].to_list(),
        #     index=out_df['pos'],
        #     columns=self.w_pca_obj.sample_lst,
        #     dtype='Int64',
        # ).round(3)

        # stats: window index, PCs variance explained, number of variants
        self.stat_df = out_df[
            ['pos', 'w_start', 'w_stop', 'w_size', 'n_var',]
            + [f'pc_{i}_ve' for i in range(1, 11)]
        ].copy().round(3)
        self.stat_df = self.stat_df.reset_index()
        self.stat_df = self.stat_df.rename(columns={'index': 'w_idx'})
        self.stat_df = self.stat_df.set_index('pos')


    def from_files(self):
        '''
        Compile from existing files
        '''

        # print info
        log.newline()
        log.info(f'Reading data from prefix "{self.prefix}*"')

        # flush existing instance variables
        for attr in list(vars(self)):
            if attr.endswith('_df'):
                delattr(self, attr)

        # pc_dfs
        for pc in self.pcs:
            setattr(
                self,
                f'pc_{pc}_df',
                pd.read_csv(
                    f'{self.prefix}.pc_{pc}.tsv.gz',
                    sep='\t',
                    index_col='pos'
                )
            )

        # hetp
        self.hetp_df = pd.read_csv(
            f'{self.prefix}.hetp.tsv.gz',
            sep='\t',
            index_col='pos',
        )

        # # miss
        # self.miss_df = pd.read_csv(
        #     f'{self.prefix}.miss.tsv.gz}',
        #     sep='\t',
        #     index_col='pos',
        #     dtype='Int64'
        # )

        # stat
        self.stat_df = pd.read_csv(
            f'{self.prefix}.stat.tsv.gz',
            sep='\t',
            index_col='pos',
        )


    def to_files(self):
        '''
        Write output (pc_dfs, hetp_df, miss_df, stat_df) to files
        '''

        # create output directory if prefix contains '/'
        if '/' in self.prefix:
            path = f'{"/".join(self.prefix.split("/")[0:-1])}/'
            if not os.path.exists(path):
                os.makedirs(path)

        # pc_dfs
        for pc in self.pcs:
            getattr(self, f'pc_{pc}_df').to_csv(
                f'{self.prefix}.pc_{pc}.tsv.gz',
                sep='\t', index_label='pos', na_rep='NA'
            )

        # hetp
        self.hetp_df.to_csv(
            f'{self.prefix}.hetp.tsv.gz',
            sep='\t', index_label='pos', na_rep='NA')

        # # miss
        # self.miss_df.to_csv(
        #     f'{self.prefix}.miss.tsv.gz',
        #     sep='\t', index_label='pos', na_rep='NA')

        # stat
        self.stat_df.to_csv(
            f'{self.prefix}.stat.tsv.gz',
            sep='\t', index_label='pos', na_rep='NA')


    def modify_data(self, attr_name, mod_func, *args):
        '''
        Modify an attribute of Data instance with a supplied function
        '''
        # get attribute
        attr = getattr(self, attr_name)

        # modify
        attr_mod = mod_func(attr, *args)

        # update attribute
        setattr(self, attr_name, attr_mod)
