'''
Modifications of PC data.
'''

# IMPORT PACKAGES
import numpy as np


# CLASSES

class Polarize:
    '''
    Polarize PC data: 
    (1) Adaptive: using the n samples with highest local absolute values for 
        polarization
    (2) Guide Samples: Use one ore more user-specified guide-samples to 
        polarize along the entire sequence.
    '''

    def adaptive(self, pc_df, n_prev_windows):
        '''
        Polarize a dataframe using the signs of the most frequent polarizer 
        (individual with largest absolute value) across n_prev_windows to inform
        whether to flip the current window.
        '''

        def is_positive(x):
            return 1 if x > 0 else 0

        def invert(polarizer_arr):
            '''
            Evaluate an array with an arbitrary number of input windows (rows), 
            columns being the PC values of each sample. The first window is the
            current window, all others are previous windows that are already
            polarized. First, determine the polarizer (largest absolute PC 
            value) of the current window, then the most frequent polarizer
            across the previous windows. Then record the polarization decision
            of the current window against the included previous windows and 
            output the most frequent decision.
            '''

            # get polarizer individual of the current window
            polarizer = np.argmax(np.abs(polarizer_arr[0]))

            # get polarizer individual of the previous window
            prev_polarizers = [
                np.argmax(np.abs(polarizer_arr[x])) \
                    for x in range(1, len(polarizer_arr))
            ]
            prev_polarizer = max(
                set(prev_polarizers), key=prev_polarizers.count
            )

            # initiate polarization list
            pol_lst = []

            # case: polarizer remains the same; record for each included
            # previous window, if it has the same sign as the current
            if polarizer == prev_polarizer:

                # get sign of current window
                win_pos = is_positive(polarizer_arr[0][polarizer])

                # get signs of previous windows and append 0 to pol_lst if a 
                # window has the same sign as the current window, else append 1
                for i in range(1, len(polarizer_arr)):
                    prev_win_pos = is_positive(polarizer_arr[i][polarizer])
                    if prev_win_pos == win_pos:
                        pol_lst.append(0)  
                    else:
                        pol_lst.append(1)

            # case: polarizer changes; record for each included
            # previous window, if it has the same sign as the current
            if polarizer != prev_polarizer:
                win_pos = is_positive(polarizer_arr[0][prev_polarizer])

                # get signs of previous windows and append 0 to pol_lst if a 
                # window has the same sign as the current window, else append 1
                for i in range(1, len(polarizer_arr)):
                    prev_win_pos = is_positive(polarizer_arr[i][prev_polarizer])
                    if prev_win_pos == win_pos:
                        pol_lst.append(0)                        
                    else:
                        pol_lst.append(1)
            
            # evaluate: return True if more previous windows had the opposite
            # sign than the current window, else return False
            inv = True if sum(pol_lst) > len(pol_lst)/2 else False

            return(inv)

        # iterate across pc_df, evaluating rolling meta-window with previous 
        # windows one by one
        for i, (idx, row) in enumerate(pc_df.iterrows()):

            # window 1: initiate polarizer_arr
            if i == 0:
                polarizer_arr = np.array(
                    [
                    np.array(pc_df.head(1))[0],
                    ]
                )

                continue

            # windows 2 - $n_prev_windows: append rows and polarize with 
            # available rows
            if i >= 1 and i <= n_prev_windows:
                polarizer_arr = np.vstack([np.array(row), polarizer_arr])
                if invert(polarizer_arr):
                    polarizer_arr[0] = polarizer_arr[0] * -1
                    pc_df.loc[idx] = row * -1
            
            # all remaining windows: append and drop a row with each iteration
            else:
                polarizer_arr = np.vstack([np.array(row), polarizer_arr])[:-1]
                if invert(polarizer_arr):
                    polarizer_arr[0] = polarizer_arr[0] * -1
                    pc_df.loc[idx] = row * -1
        
        return pc_df


    def guide_samples(self, pc_df, guide_sample_lst):
        '''
        Use fixed guide samples to polarize along a chromosome.
        '''

        # list that will hold per-window flip instructions for each guide sample
        gs_flip_lst = []

        # for each guide_sample, get a list of pc values per window
        for gs in guide_sample_lst:
            gs_window_lst =  list(pc_df[[gs]][gs])

            # first window has no reference/prev_window --> If second window is ### REVIEW THIS LOGIC
            # None set prev_window window to 0 for numerical comparison
            flip = [0]
            prev_window = \
                gs_window_lst[0] if not gs_window_lst[1] == None  else 0 
            
            # parse each window separately and compare to prev_window
            for window in gs_window_lst[1:]:

                # window has no value --> 0 (=no effect on balance)
                if window == None:
                    flip.append(0)
                    continue
                
                # window closer to flipped prev_window --> 1 (=flip)
                elif abs(window - prev_window) > abs(window - (prev_window*-1)):
                    flip.append(1)
                    prev_window = (window*-1) # flip prev_window
                
                # window closer to 'unflipped' prev_window --> -1 (=don't flip)
                else:
                    flip.append(-1)
                    prev_window = window # don't flip prev_window
            
            # append to gs_flip_lst
            gs_flip_lst.append(flip)

        # create array of the per sample flip instructions
        gs_flip_arr = np.array(gs_flip_lst, dtype=int).transpose()

        # get flip consensus from all guide samples (positive row sum --> flip)
        consensus_flip_lst = list(gs_flip_arr.sum(axis=1))

        # flip windows according to consensus_flip_lst (flip if positive)
        for idx, val in enumerate(consensus_flip_lst):
            if val > 0:
                pc_df.iloc[idx] *= -1
        
        return pc_df



class Flip:
    '''
    Flip all PC values of a chromosome or specific windows.
    '''

    def flip_chrom(self, pc_df):
        '''
        Flip all PC values for a chromosome.
        '''
        
        return pc_df * -1
    
    def flip_windows(self, pc_df, coord_lst):
        '''
        Flip specifiec windows.
        '''
        
        pass

    
        



