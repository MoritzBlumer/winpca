'''
Central class to store, read, write and modify PCA and associated data.
'''

# IMPORT PACKAGES
import sys
import pandas as pd


## CLASSES

class wpca_data:
    '''
    Central container for PCA and associated data. Compiles object from 
    existing output files (via prefix) or directly from w_pca instance. Data
    can be modified through indluded modifier methods.
    '''

    def __init__(self, prefix, w_pca_obj=None):

        # input variables
        self.prefix = prefix
        self.w_pca_obj = w_pca_obj

        # instance variables
        self.pc_1 = None
        self.pc_2 = None
        self.hetp = None
        self.miss = None
        self.stat = None

        # set file suffixes
        self.suffix_dct = {
            'pc_1': 'pc_1.tsv.gz',
            'pc_2': 'pc_2.tsv.gz',
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

        # pc_1
        self.pc_1 = pd.DataFrame(
            data=out_df['pc_1'].to_list(),
            index=out_df['pos'],
            columns=self.w_pca_obj.sample_lst,
        )

        # pc_2
        self.pc_2 = pd.DataFrame(
            data=out_df['pc_2'].to_list(),
            index=out_df['pos'],
            columns=self.w_pca_obj.sample_lst,
        )

        # hetp
        self.hetp = pd.DataFrame(
            data=out_df['hetp'].to_list(),
            index=out_df['pos'],
            columns=self.w_pca_obj.sample_lst,
        )

        # miss
        self.miss = pd.DataFrame(
            data=out_df['n_miss'].to_list(),
            index=out_df['pos'],
            columns=self.w_pca_obj.sample_lst,
            dtype='Int64',
        )

        # stats: window index, pc_1/pc_2 variance explained, number of variants
        self.stat = out_df[['pos', 'n_var', 'pc_1_ve', 'pc_2_ve', ]].copy()
        self.stat = self.stat.reset_index().rename(columns={'index': 'w_idx'})
        self.stat = self.stat.set_index('pos')


    def from_files(self):
        '''
        Compile from out_dct of w_pca instance.
        '''

        # print info
        print(
            f'\n[INFO] Reading data from prefix "{self.prefix}*".',
            file=sys.stderr, flush=True,
        )

        # flush existing instance variables
        self.pc_1, self.pc_2, self.hetp, self.miss, self.stat = \
            None, None, None, None, None

        # pc_1
        self.pc_1 = pd.read_csv(
            f'{self.prefix}.{self.suffix_dct["pc_2"]}',
            sep='\t', 
            index_col='pos',
        )

        # pc_2
        self.pc_2 = pd.read_csv(
            f'{self.prefix}.{self.suffix_dct["pc_1"]}',
            sep='\t', 
            index_col='pos',
        )

        # hetp
        self.hetp = pd.read_csv(
            f'{self.prefix}.{self.suffix_dct["hetp"]}',
            sep='\t', 
            index_col='pos',
        )

        # miss
        self.miss = pd.read_csv(
            f'{self.prefix}.{self.suffix_dct["miss"]}',
            sep='\t', 
            index_col='pos',
            dtype='Int64'
        )

        # stat
        self.stat = pd.read_csv(
            f'{self.prefix}.{self.suffix_dct["stat"]}',
            sep='\t', 
            index_col='pos',
        )

        # [...]
        pass


    def to_files(self):
        '''
        Write output (pc_1, pc_2, hetp, miss, stat) to files.
        '''

        # pc_1
        self.pc_1.to_csv(
            f'{self.prefix}.{self.suffix_dct["pc_1"]}', 
            sep='\t', index_label='pos', na_rep='NA')
        
        # pc_2
        self.pc_2.to_csv(
            f'{self.prefix}.{self.suffix_dct["pc_2"]}', 
            sep='\t', index_label='pos', na_rep='NA')

        # hetp
        self.hetp.to_csv(
            f'{self.prefix}.{self.suffix_dct["hetp"]}', 
            sep='\t', index_label='pos', na_rep='NA')

        # miss
        self.miss.to_csv(
            f'{self.prefix}.{self.suffix_dct["miss"]}', 
            sep='\t', index_label='pos', na_rep='NA')

        # stat
        self.stat.to_csv(
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