'''
Modifications of PC data.
'''


## IMPORT PACKAGES
import numpy as np

## MODULES
from .log import Log

## INSTANTIATE LOGGER
log = Log()


## CLASSES

class Polarize:
    '''
    Polarize PC data:
    (1) Adaptive: using the n samples with highest local absolute values for
        polarization
    (2) Guide Samples: use one ore more user-specified guide-samples to
        polarize along the entire sequence
    (3) Local ML model: transformer trained for optimal local polarization,
        applied after adaptive by default
    (4) Global ML model: Convolutional network trained for optimal global
        polarization, applied after local ML model
    '''

    @staticmethod
    def adaptive(pc_df, n_prev_windows):
        '''
        Polarize a dataframe using the signs of the most frequent polarizer
        (individual with largest absolute value) across n_prev_windows to
        inform whether to flip the current window
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
            output the most frequent decision
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
                    prev_win_pos = is_positive(
                        polarizer_arr[i][prev_polarizer]
                    )
                    if prev_win_pos == win_pos:
                        pol_lst.append(0)
                    else:
                        pol_lst.append(1)

            # evaluate: return True if more previous windows had the opposite
            # sign than the current window, else return False
            inv = sum(pol_lst) > len(pol_lst) / 2

            return inv

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
            if i in range(1, n_prev_windows + 1):
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


    @staticmethod
    def guide_samples(pc_df, guide_sample_lst):
        '''
        Use fixed guide samples to polarize along a chromosome
        '''

        # list that will hold per-window flip instructions per guide sample
        gs_flip_lst = []

        # for each guide_sample, get a list of pc values per window
        for gs in guide_sample_lst:
            gs_window_lst =  list(pc_df[[gs]][gs])

            # first window has no reference/prev_window --> If second window is
            # None set prev_window window to 0 for numerical comparison
            flip = [0]
            prev_window = \
                gs_window_lst[0] if not gs_window_lst[1] is None  else 0

            # parse each window separately and compare to prev_window
            for window in gs_window_lst[1:]:

                # window has no value --> 0 (=no effect on balance)
                if window is None:
                    flip.append(0)
                    continue

                # window closer to flipped prev_window --> 1 (=flip)
                if abs(window - prev_window) > abs(window - (prev_window*-1)):
                    flip.append(1)
                    prev_window = window*-1 # flip prev_window

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


    @staticmethod
    def local_ml(pc_df, n_draws, flip_threshold, weights_path, epsilon):
        '''
        Transformer that sees up to 2048 windows at a time and decides based on
        flip_threshold whether any two neighboring windows agree or disagree
        in sign polarity - if not, all windows to the right are flipped. The
        model is run n_draws times for 64 randomly drawn individuals, and
        each flip decision is based on the average of these runs. If the
        dataset has <64 individuals, samples can be included more than once. If
        a chromosome has > 2048 windows, this is done in 2048 window chunks
        with 1024 window overlap.
        '''

        def normalize(w_arr, scale_arr, shift_arr, epsilon):
            '''
            Normalize each row to mean 0 and variance 1, then multiply element-
            wise with provided scales/weights and add shifts/biases
            (reimplementing torch.nn.LayerNorm; epsilon prevents division by
            zero)
            '''
            mean_arr = w_arr.mean(-1, keepdims=True)
            var_arr = w_arr.var(-1, keepdims=True)
            norm_arr = (w_arr - mean_arr) / np.sqrt(var_arr + epsilon) \
                * scale_arr + shift_arr

            return norm_arr

        def edge_logits(draw_arr, epsilon, n_heads=8):
            '''
            One pass of the transformer over draw_arr (windows x 64 normalized
            PC values of one draw): outputs one logit per pair of neighbouring
            windows (positive = sign break more likely than not)
            '''

            # compute cosine similarity of neighbouring windows, where
            # right_arr is shifted one window compared to left_arr to
            # facilitate neighbor comparisons
            left_arr, right_arr = draw_arr[:-1], draw_arr[1:]
            cos_arr = (left_arr * right_arr).sum(1) \
                / (np.linalg.norm(left_arr, axis=1)
                   * np.linalg.norm(right_arr, axis=1) + epsilon)

            # normalize and apply learned weights and bias
            norm_arr = normalize(
                draw_arr,
                weight_dct['input_norm.weight'],
                weight_dct['input_norm.bias'],
                epsilon,
            )

            # compute 128 input features per window from normalized PC data:
            # - matrix multiplication of norm_arr (up to 2048 n_windows x 64
            #   individuals) with learned weights (128 x 64) --> n_windows x
            #   128 features
            # - adding the same 128 learned biases per window
            feature_arr = \
                norm_arr @ weight_dct['input_proj.weight'].T \
                + weight_dct['input_proj.bias']

            # add a deterministic/algorithmically defined pattern to each
            # window (weight_dct['pos_enc.pe']) that encodes/imprints its
            # position within the current chromosome chunk to embed relative
            # position into the model (scaled by a single learned factor that
            # is applied to all windows: weight_dct['pos_enc.scale'])
            feature_arr = \
                feature_arr + weight_dct['pos_enc.scale'] \
                * weight_dct['pos_enc.pe'][0, :len(feature_arr)]

            # extract number of features
            n_features = feature_arr.shape[1]

            # each head is one of n_heads=8 components into which 128 features
            # are divided during attention, which allows the model to learn
            # different signals/properties of the input data - head size is
            # thus 16 features (128/8)
            head_size = n_features // n_heads

            # 4 layers with self-attention (8 heads) followed by feed-forward:
            # (1) self-attention (individual windows become aware of their
            #     context): per layer, windows are first normalized then
            #     transformed into query/key/value arrays, split into 8 heads;
            #     per head, each window collects information from all windows;
            #     heads are unified and blended and features are updated to
            #     reflect the new context
            # (2) feed-forward (weighting new context information): the 128
            #     updated features are normalized, processed by multiplying the
            #     128 input features with 512 x 128 learned weights, 512
            #     learned biases are added and ReLU applied; the 512 features
            #     are converted back into 128 features by multiplying them with
            #     128x512 learned weights and 128 biases
            # (3) These 128 new features are added to the input 128 feature
            #     array
            for i in range(n_layers):

                # helper lookup function for layer weights
                layer_prefix = f'transformer.layers.{i}.'
                layer_weight = lambda name: weight_dct[layer_prefix + name]

                # normalize copy of current 128 features and apply the layer's
                # scales and shifts
                norm_feature_arr = normalize(
                    feature_arr,
                    layer_weight('norm1.weight'),
                    layer_weight('norm1.bias'),
                    epsilon,
                )

                # process each window and produce query, key, value, applying
                # the layer's learned weights and biases
                att_query_arr, att_key_arr, att_value_arr = np.split(
                    norm_feature_arr @ \
                    layer_weight('self_attn.in_proj_weight').T
                    + layer_weight('self_attn.in_proj_bias'),
                    3,
                    axis=1,
                )

                # reshape into 8x16 heads
                att_query_arr = att_query_arr.reshape(-1, n_heads, head_size) \
                    .transpose(1, 0, 2)
                att_key_arr = att_key_arr.reshape(-1, n_heads, head_size) \
                    .transpose(1, 0, 2)
                att_value_arr = att_value_arr.reshape(-1, n_heads, head_size) \
                    .transpose(1, 0, 2)

                # compute attention scores as dot product of every window's
                # query with every window's key (per head) and scale the result
                # down by sqrt(head_size)
                att_score_arr = \
                    att_query_arr @ att_key_arr.transpose(0, 2, 1) \
                    / np.sqrt(head_size)

                # rescale attention scores and give large values more weight:
                # subtract max score, exponentiate and divide by the sum
                # (softmaxing)
                att_score_arr = np.exp(
                    att_score_arr - att_score_arr.max(-1, keepdims=True)
                )
                att_score_arr /= att_score_arr.sum(-1, keepdims=True)

                # process attention through matrix multiplication of attention
                # scores with the corresponding attention values, then re-join
                # the 8 heads
                att_arr = \
                    (att_score_arr @ att_value_arr).transpose(1, 0, 2) \
                    .reshape(-1, n_features)

                # apply the layer's weights and biases to the attention
                # attributes (which combines the information from the re-joined
                # heads) and add them to the current features
                feature_arr = feature_arr + att_arr \
                    @ layer_weight('self_attn.out_proj.weight').T \
                    + layer_weight('self_attn.out_proj.bias')

                # initiate feed-forward: normalize copy of current 128 features
                # and apply the layer's scales and shifts
                ff_arr = normalize(
                    feature_arr,
                    layer_weight('norm2.weight'),
                    layer_weight('norm2.bias'),
                    epsilon,
                )

                # expand to 512 features by matrix multiplication with the
                # layer's weights and adding biases, then ReLU (negative values
                # are set to 0) - effectively this is another feature-rich
                # digestion of the newly gained context
                ff_arr = np.maximum(
                    ff_arr @ layer_weight('linear1.weight').T
                    + layer_weight('linear1.bias'),
                    0,
                )

                # shrink back to 128 features by applying the layer's weights
                # and biases, then add the result to the current features
                feature_arr = feature_arr \
                    + ff_arr @ layer_weight('linear2.weight').T \
                    + layer_weight('linear2.bias')

            # edge head - compile edge array containing 385 features per edge:
            # (1) feature array except the last window (n=128)
            # (2) feature array except the first window (n=128)
            # (3) edge-wise difference (1) - (2) (n=128)
            # (4) edge-wise cosine similarity array (n=1)
            edge_arr = np.concatenate(
                [feature_arr[:-1], feature_arr[1:],
                 feature_arr[:-1] - feature_arr[1:], cos_arr[:, None]], axis=1
            )

            # shrink 385 to 64 features by applying the edge head's weights
            # and biases, then ReLU (negative values are set to 0)
            edge_arr = np.maximum(
                edge_arr @ weight_dct['edge_head.0.weight'].T
                + weight_dct['edge_head.0.bias'], 0
            )

            # assessment: shrink to one value per edge by applying learned
            # weights and bias, and return single positive (flip) or negative
            # (don't flip) logit per edge
            logit_arr = (edge_arr @ weight_dct['edge_head.2.weight'].T
                    + weight_dct['edge_head.2.bias'])[:, 0]

            return logit_arr

        def edge_probs(pc_arr, seed, epsilon, n_inds=64, max_len=2048):
            '''
            Normalize input data, draw samples and run transformer, convert
            edge logits into edge probabilities
            '''

            # normalize to mean 0 and variance 1 (epsilon prevents division by
            # zero)
            norm_arr = ((pc_arr - pc_arr.mean()) / (pc_arr.std() + epsilon)) \
                .astype(np.float32)

            # draw random samples
            n_windows, n_samples = norm_arr.shape
            sample_idx_arr = np.random.default_rng(seed).choice(
                n_samples, n_inds, replace=n_samples < n_inds
            )
            draw_arr = norm_arr[:, sample_idx_arr]

            # run transformer

            # single run for up to 2,048 windows
            if n_windows <= max_len:
                logit_arr = edge_logits(draw_arr, epsilon)

            # multiple runs if >2,048 windows
            else:

                # calculate stride (step) size
                stride = max_len // 2

                # calculate start positions
                start_lst = list(
                    range(
                        0,
                        n_windows - max_len + 1,
                        stride,
                    )
                )

                # add one more chunk if last window is not covered
                if not start_lst or start_lst[-1] + max_len < n_windows:
                    start_lst.append(max(0, n_windows - max_len))

                # initiate array which will hold the sum of logits per edge
                logit_sum_arr = np.zeros(n_windows - 1)

                # initiate array that will hold the number of logits per edge
                n_logit_arr = np.zeros(n_windows - 1)

                # run transformer chunk-wise, updating logit_sum_arr and
                # n_logit_arr
                for start in start_lst:
                    logit_sum_arr[start:start + max_len - 1] += \
                        edge_logits(draw_arr[start:start + max_len], epsilon)
                    n_logit_arr[start:start + max_len - 1] += 1

                # calculate average logit per edge
                logit_arr = logit_sum_arr / np.maximum(n_logit_arr, 1)

            # convert logits into probabilities
            prob_arr = 1.0 / (1.0 + np.exp(-logit_arr))

            return prob_arr

        # load model
        with np.load(weights_path) as npz_file:
            weight_dct = {key: npz_file[key] for key in npz_file.files}

        # extract number of layers
        n_layers = len({key.split('.')[2] for key in weight_dct
                        if key.startswith('transformer.layers.')})

        # convert pc values to array
        pc_arr = pc_df.fillna(0).values.astype(np.float32)

        # run transformer n_draws times using draw index as seed
        prob_arr_lst = []
        for i in range(n_draws):
            prob_arr = edge_probs(
                pc_arr,
                i,
                epsilon,
            )
            prob_arr_lst.append(
                prob_arr,
            )

        # calculate mean across draws
        mean_prob_arr = np.mean(
            prob_arr_lst,
            axis=0,
        )

        # evaluate edges where flip threshold is exceeded
        trigger_arr = mean_prob_arr > flip_threshold

        # calculate cumulative sum of preceding sign breaks per window
        n_trigger_arr = np.concatenate(
            [
                [0],
                np.cumsum(trigger_arr),
            ],
        )

        # derive windows to flip, i.e. those with an odd number of preceding
        # sign flips
        flip_arr = n_trigger_arr % 2

        # flip windows according to flip_arr (flip if 1)
        for idx, val in enumerate(flip_arr):
            if val == 1:
                pc_df.iloc[idx] *= -1

        return pc_df


    @staticmethod
    def global_ml(pc_df, weights_path, epsilon):
        '''
        Convolutional network that sees the entire chromosome at once, computes
        20 features per non-empty window, applies 9 convolutional blocks with
        3 taps each at increasing dilation (1, 2, 4, 8, 16, 32, 64, 128, 256)
        to find the best flip option per window (considering windows and edges)
        using viterbi
        '''

        def cos_sim(left_arr, right_arr, epsilon):
            '''
            Cosine similarity of corresponding rows of left_arr and right_arr
            (one value per row: +1 = same pattern, -1 = flipped, 0 =
            unrelated)
            '''

            # calculate dot product per row
            dot_arr = (left_arr * right_arr).sum(1)

            # multiply the magnitudes of both rows (epsilon prevents division
            # by zero)
            magnitude_product_arr = np.linalg.norm(left_arr, axis=1) \
                * np.linalg.norm(right_arr, axis=1) + epsilon

            # calculate cosine similarity
            cos_sim_arr = dot_arr / magnitude_product_arr

            return cos_sim_arr

        def get_consensus(pc_arr, epsilon):
            '''
            Compute a reference consensus pattern across all individuals based
            on all windows (in a way, a scaled, single representative window)
            '''

            # calculate mean PC score per individual across all windows
            mean_arr = pc_arr.mean(axis=0)

            # calculate the magnitude of the mean across all individuals
            mean_magnitude = np.linalg.norm(mean_arr)

            # scale PC mean by mean magnitude (epsilon prevents division by
            # zero)
            cons_arr = mean_arr / (mean_magnitude + epsilon)

            # refine consensus pattern 8 times
            for _ in range(8):

                # score current consensus by matrix multiplication with raw PC
                # values
                score_arr = pc_arr @ cons_arr

                # update consensus pattern by matrix multiplication of scores
                # with transposed PC values
                cons_arr = pc_arr.T @ score_arr

                # rescale to magnitude 1 (epsilon prevents division by zero)
                cons_arr /= np.linalg.norm(cons_arr) + epsilon

            # flip sign if they disagree with raw mean (epsilon prevents
            # sign to be exactly zero)
            cons_arr = cons_arr * np.sign(
                cons_arr @ mean_arr + epsilon,
            )

            return cons_arr

        def mean_cons_similarity(cons_sim_arr, k):
            '''
            Running mean of local consensus similarity over +/- k windows
            '''

            # calculate running cumulative sum of preceding windows and add a
            # leading 0.0 as offset
            cumsum_arr = np.concatenate(
                [[0.0], np.cumsum(cons_sim_arr)],
            )

            # fetch window indices
            idx_arr = np.arange(len(cons_sim_arr))

            # fetch region start and end, capped at sequence ends
            start_arr = np.maximum(0, idx_arr - k)
            end_arr = np.minimum(len(idx_arr), idx_arr + k + 1)

            # calculate mean similarity per window
            local_cons_sim_arr = \
                (cumsum_arr[end_arr] - cumsum_arr[start_arr]) \
                / (end_arr - start_arr)

            return local_cons_sim_arr

        def local_similarity(scaled_pc_arr, k, epsilon):
            '''
            Similarity of each window to the summed pattern of +/- k
            surrounding windows
            '''

            # fetch number of windows and samples
            n_windows, n_samples = scaled_pc_arr.shape

            # calculate running cumulative sum per individual and add a
            # leading row of zeros as offset (row i = sum of all windows
            # before window i)
            zero_row_arr = np.zeros((1, n_samples))
            cumsum_arr = np.vstack(
                [zero_row_arr, np.cumsum(scaled_pc_arr, axis=0)],
            )

            # fetch window indices
            idx_arr = np.arange(n_windows)

            # fetch region start and end, capped at sequence ends
            start_arr = np.maximum(0, idx_arr - k)
            end_arr = np.minimum(n_windows, idx_arr + k + 1)

            # calculate summed pattern per window (one row of n_samples values
            # per window)
            region_arr = cumsum_arr[end_arr] - cumsum_arr[start_arr]

            # compare each window's similarity with surrounding patterns
            local_sim_arr = cos_sim(scaled_pc_arr, region_arr, epsilon)

            return local_sim_arr

        def block_similarity(scaled_pc_arr, k, epsilon):
            '''
            Similarity of the mean pattern of the k windows left of each edge
            to the mean pattern of the k windows right of it (one value per
            window; 0 for the last window, which has no edge)
            '''

            # fetch number of windows and samples
            n_windows, n_samples = scaled_pc_arr.shape

            # calculate running cumulative sum per individual and add a
            # leading row of zeros as offset (row i = sum of all windows
            # before window i)
            zero_row_arr = np.zeros((1, n_samples))
            cumsum_arr = np.vstack(
                [zero_row_arr, np.cumsum(scaled_pc_arr, axis=0)],
            )

            # fetch edge indices (edge i lies between window i and i + 1)
            edge_arr = np.arange(n_windows - 1)

            # fetch start of the left block and end of the right block, capped
            # at sequence ends (the edge itself is at edge_arr + 1)
            start_arr = np.maximum(0, edge_arr + 1 - k)
            end_arr = np.minimum(n_windows, edge_arr + 1 + k)

            # calculate mean pattern of the left block (windows start to edge)
            left_count_arr = edge_arr + 1 - start_arr
            left_mean_arr = \
                (cumsum_arr[edge_arr + 1] - cumsum_arr[start_arr]) \
                / left_count_arr[:, None]

            # calculate mean pattern of the right block (windows after the
            # edge to end)
            right_count_arr = end_arr - edge_arr - 1
            right_mean_arr = \
                (cumsum_arr[end_arr] - cumsum_arr[edge_arr + 1]) \
                / right_count_arr[:, None]

            # compare left and right block, add 0 for the last window
            block_sim_arr = cos_sim(left_mean_arr, right_mean_arr, epsilon)
            block_sim_arr = np.concatenate([block_sim_arr, [0.0]])

            return block_sim_arr

        def compute_features(pc_arr, data_idx_arr, epsilon, max_samples=512):
            '''
            Compute 20 features per window and consensus similarity per window
            '''

            # if more than max_samples, draw max_samples random individuals
            n_windows, n_samples = pc_arr.shape
            if n_samples > max_samples:
                pc_arr = pc_arr[
                    :,
                    np.random.default_rng(0).choice(
                        n_samples, max_samples, replace=False
                    ),
                ]

            # compute magnitude of each window's vector of PC values
            magnitude_arr = np.linalg.norm(pc_arr, axis=1)

            # scale PC values per window to magnitude = 1 (epsilon prevents
            # division by zero)
            scaled_pc_arr = pc_arr / (magnitude_arr[:, None] + epsilon)

            # fetch reference consensus pattern
            cons_arr = get_consensus(pc_arr, epsilon)

            # score each window's match with the consensus by matrix
            # multiplication with raw PC values (one dot product per window)
            score_arr = pc_arr @ cons_arr

            # calculate the magnitude of the consensus across all individuals
            cons_magnitude = np.linalg.norm(cons_arr)

            # multiply PC magnitudes with the consensus magnitude (epsilon
            # prevents division by zero)
            magnitude_product_arr = magnitude_arr * cons_magnitude + epsilon

            # compute cosine similarity per window: +1 = same pattern,
            # -1 = flipped, 0 = unrelated
            cons_sim_arr = score_arr / magnitude_product_arr

            # calculate median magnitude across non-empty windows and set to
            # 1.0 if all windows are empty
            nonempty_magnitude_arr = magnitude_arr[magnitude_arr > 0]
            if len(nonempty_magnitude_arr) > 0:
                median_magnitude = np.median(nonempty_magnitude_arr)
            else:
                median_magnitude = 1.0

            # below, calculate 20 features (properties of the data) and add
            # them to a list
            feature_lst = []

            # (1) absolute per-window similarity with the consensus (ignoring
            # orientation)
            abs_sim_arr = np.abs(cons_sim_arr)
            feature_lst.append(abs_sim_arr)

            # (2, 3, 4) sign consistency metric at three local scales: among
            # 16, 128 and 1024 surrounding windows
            for k in (8, 64, 512):

                # calculate mean local sign consistency over +/- k windows
                local_cons_sim_arr = mean_cons_similarity(cons_sim_arr, k)

                # multiply per window with local mean consensus similarity
                cons_sim_agree_arr = cons_sim_arr * local_cons_sim_arr

                # per window score from +1 to -1, where +1 is perfect agreement
                # with local context, -1 perfect disagreement with local
                # context
                feature_lst.append(cons_sim_agree_arr)

            # (5) log magnitude relative to the median magnitude: 0 = typical,
            # positive = stronger, negative = weaker (offset prevents log(0)
            # for empty windows)
            offset = 1e-6
            rel_magnitude_arr = (magnitude_arr + offset * median_magnitude) \
                / median_magnitude
            feature_lst.append(np.log(rel_magnitude_arr))

            # (6) flag empty windows as 1, others as 0 (this was implemented
            # before the model was restricted to non-empty windows and will
            # thus always be 0, but is kept to match the trained model)
            empty_arr = (magnitude_arr == 0).astype(float)
            feature_lst.append(empty_arr)

            # (7, 8) similarity to the summed pattern of surrounding 64 and 512
            # windows
            for k in (32, 256):
                local_sim_arr = local_similarity(scaled_pc_arr, k, epsilon)
                feature_lst.append(local_sim_arr)

            # (9 to 12) position along the chromosome and size of the dataset

            # fetch window indices
            idx_arr = np.arange(n_windows)

            # calculate relative position along the chromosome: 0 = first
            # window, 1 = last window (max prevents division by zero for a
            # single window)
            rel_pos_arr = idx_arr / max(n_windows - 1, 1)

            # (9) relative distance to the nearest chromosome end: 0 at both
            # ends, 0.5 in the middle
            rel_end_dist_arr = np.minimum(rel_pos_arr, 1 - rel_pos_arr)
            feature_lst.append(rel_end_dist_arr)

            # (10) distance to the nearest chromosome end in windows, as log
            # (log1p = log(1 + x), so 0 at the ends)
            end_dist_arr = np.minimum(idx_arr, n_windows - 1 - idx_arr)
            log_end_dist_arr = np.log1p(end_dist_arr)
            feature_lst.append(log_end_dist_arr)

            # (11) number of windows, as log (same value for every window)
            log_n_windows_arr = np.full(n_windows, np.log(n_windows))
            feature_lst.append(log_n_windows_arr)

            # (12) number of samples, as log (same value for every window)
            log_n_samples_arr = np.full(n_windows, np.log(n_samples))
            feature_lst.append(log_n_samples_arr)

            # (13) similarity of each window's pattern to the next window's
            # pattern (one value per edge, adding 0 for the last window)
            adj_sim_arr = cos_sim(
                scaled_pc_arr[:-1],
                scaled_pc_arr[1:],
                epsilon,
            )
            adj_sim_arr = np.concatenate([adj_sim_arr, [0.0]])
            feature_lst.append(adj_sim_arr)

            # (14-17) similarity of the mean patterns of the 4, 16, 64, 256
            # windows left compared to right of each edge
            for k in (4, 16, 64, 256):
                block_sim_arr = block_similarity(scaled_pc_arr, k, epsilon)
                feature_lst.append(block_sim_arr)

            # (18) product of the consensus similarities of each window and
            # the next window: positive = same orientation relative to the
            # consensus, negative = opposite (one value per edge), adding 0 for
            # the last window
            adj_cons_sim_arr = cons_sim_arr[:-1] * cons_sim_arr[1:]
            adj_cons_sim_arr = np.concatenate([adj_cons_sim_arr, [0.0]])
            feature_lst.append(adj_cons_sim_arr)

            # (19) flag windows that have a next window (edge) with 1, the
            # last window with 0
            has_edge_arr = np.ones(n_windows - 1)
            has_edge_arr = np.concatenate([has_edge_arr, [0.0]])
            feature_lst.append(has_edge_arr)

            # (20) number of empty windows between each data window and the
            # next data window, as log (0 for the last data window)
            n_gap_arr = np.diff(data_idx_arr) - 1
            gap_arr = np.concatenate([np.log1p(n_gap_arr), [0.0]])
            feature_lst.append(gap_arr)

            # stack the 20 features into one table (windows x features)
            feature_arr = np.stack(
                feature_lst,
                axis=1,
            )

            # round to float16 precision and convert back to float32 (since
            # training used float16)
            feature_arr = feature_arr.astype(np.float16).astype(np.float32)

            return cons_sim_arr, feature_arr

        def conv1d(in_arr, weight_arr, bias_arr, dilation=1):
            '''
            1D convolution of in_arr (channels x windows): per output channel,
            a weighted sum over all input channels of the windows at
            distance-dilation, 0 and +dilation (for 3 taps), plus bias
            '''

            # fetch number of taps (=number of windows "seen", 1 for input
            # layer and head, three for core blocks)
            n_taps = weight_arr.shape[2]

            # calculate number of missing windows beyond sequence end +/-
            # dilation away for padding
            pad_size = dilation * (n_taps - 1) // 2

            # add pad_size columns of zeros at both ends
            pad_arr = np.pad(
                in_arr,
                ((0, 0),
                (pad_size, pad_size)),
            )

            # initialize empty output
            n_windows = in_arr.shape[1]
            out_arr = np.zeros((weight_arr.shape[0], n_windows))

            # per tap: take the shifted copy of the input that lines up each
            # window with its neighbour (tap 0 = left, 1 = itself, 2 = right),
            # combine its channels by matrix multiplication with the tap's
            # weights and add to the output
            for tap in range(n_taps):
                start = tap * dilation
                shifted_arr = pad_arr[:, start:start + n_windows]
                out_arr += weight_arr[:, :, tap] @ shifted_arr

            # add learned biases
            out_arr = out_arr + bias_arr[:, None]

            return out_arr

        def gelu(x_arr):
            '''
            GELU activation (smooth version of ReLU): each value is multiplied
            by the probability that a standard normal value lies below it -
            effectively this gives more weight to large positive values and
            less to negative values
            '''

            # use absolute values |x| / sqrt(2)
            y_arr = np.abs(x_arr) / np.sqrt(2.0)

            # helper value of the approximation (between 0 and 1)
            t_arr = 1.0 / (1.0 + 0.3275911 * y_arr)

            # polynomial in t with published constants (see Abramowitz &
            # Stegun, 1964)
            poly_arr = 1.061405429 * t_arr - 1.453152027
            poly_arr = poly_arr * t_arr + 1.421413741
            poly_arr = poly_arr * t_arr - 0.284496736
            poly_arr = poly_arr * t_arr + 0.254829592
            poly_arr = poly_arr * t_arr

            # error function of |x| / sqrt(2)
            erf_arr = 1.0 - poly_arr * np.exp(-y_arr * y_arr)

            # probability that a standard normal value lies below x (np.sign
            # restores the sign for negative x)
            cdf_arr = 0.5 * (1.0 + np.sign(x_arr) * erf_arr)

            # GELU: x times that probability
            gelu_arr = x_arr * cdf_arr

            return gelu_arr

        def viterbi(window_evidence_arr, edge_evidence_arr):
            '''
            Find the best global path across all flipping options considering
            all windows and edges; returns flip list per window (0=keep,
            1=flip)
            '''

            # fetch number of windows
            n_windows = len(window_evidence_arr)

            # initialize array to record the state of the previous window on
            # the best path leading there for each window and each of its two
            # states (0 = keep, 1 = flip)
            prev_state_arr = np.zeros((n_windows, 2), dtype=np.int8)

            # best scores so far if the first window is kept / flipped
            score_keep = window_evidence_arr[0]
            score_flip = -window_evidence_arr[0]

            # walk along the chromosome starting from the second window
            for i in range(1, n_windows):

                # evidence of the edge between previous and current window
                edge_evidence = edge_evidence_arr[i - 1]

                # fetch the four ways to reach this window's states:
                # same state as the previous window gains the edge evidence, a
                # change loses it
                keep_after_keep = score_keep + edge_evidence
                keep_after_flip = score_flip - edge_evidence
                flip_after_keep = score_keep - edge_evidence
                flip_after_flip = score_flip + edge_evidence

                # record which previous state gave the better score for
                # keeping current window
                if keep_after_keep >= keep_after_flip:
                    prev_state_arr[i, 0] = 0
                else:
                    prev_state_arr[i, 0] = 1

                # record which previous state gave the better score for
                # flipping current window
                if flip_after_keep >= flip_after_flip:
                    prev_state_arr[i, 1] = 0
                else:
                    prev_state_arr[i, 1] = 1

                # fetch best scores up to and including current window
                score_keep = max(keep_after_keep, keep_after_flip) \
                    + window_evidence_arr[i]
                score_flip = max(flip_after_keep, flip_after_flip) \
                    - window_evidence_arr[i]

            # last window: choose the state with the better total score
            if score_keep >= score_flip:
                state = 0
            else:
                state = 1

            # reconstruct the best-scoring path from last to first window
            flip_arr = np.zeros(n_windows, dtype=int)
            for i in range(n_windows - 1, -1, -1):
                flip_arr[i] = state
                state = prev_state_arr[i, state]

            return flip_arr

        # load model
        with np.load(weights_path) as npz_file:
            weight_dct = {key: npz_file[key] for key in npz_file.files}

        # extract number of blocks
        n_blocks = len(
            {key.split('.')[1] for key in weight_dct
             if key.startswith('blocks.')}
        )

        # convert pc values to array
        pc_arr = pc_df.fillna(0).values.astype(np.float64)

        # find windows that are not gaps
        data_idx_arr = np.flatnonzero(
            np.abs(pc_arr).sum(1) > 0,
        )

        # initialize flip array which is 0 for every window (flip will be 1)
        flip_arr = np.zeros(len(pc_arr), dtype=int)

        # if at least two data windows in input, compute per window
        # - 1 consensus similarity value
        # - 20 features
        if len(data_idx_arr) >= 2:
            cons_sim_arr, feature_arr = compute_features(
                pc_arr[data_idx_arr],
                data_idx_arr,
                epsilon,
            )

            # fetch neighbour similarity (feature 13, 'adj') before
            # standardizing
            adj_idx = list(weight_dct['names']).index('adj')
            adj_arr = feature_arr[:, adj_idx]

            # standardize each feature with mean and standard deviation
            # inferred from training data
            std_feature_arr = (feature_arr - weight_dct['mu']) \
                / weight_dct['sd']

            # apply the input layer by combining the 20 standardized features
            # into 32 channels per window
            channel_arr = conv1d(
                std_feature_arr.T,
                weight_dct['inp.weight'],
                weight_dct['inp.bias'],
            )

            # apply the blocks: each block reads all 32 channels at the focal
            # window and +/- dilation windows and adds its findings to the
            # channels
            for i in range(n_blocks):

                # normalize a copy of the channels: mean 0 and variance 1
                # over all channels and windows together (epsilon prevents
                # division by zero)
                mean = channel_arr.mean()
                var = channel_arr.var()
                norm_channel_arr = \
                    (channel_arr - mean) / np.sqrt(var + epsilon)

                # apply learned scale and shift per channel
                norm_channel_arr = norm_channel_arr \
                    * weight_dct[f'norms.{i}.weight'][:, None] \
                    + weight_dct[f'norms.{i}.bias'][:, None]

                # GELU activation: keep positive values, damp negative ones
                act_channel_arr = gelu(norm_channel_arr)

                # dilated convolution: read all 32 channels at the focal window
                # and +/- dilation windows, mix them into 32 new channels
                dilation = int(weight_dct['dilations'][i])
                block_out_arr = conv1d(
                    act_channel_arr,
                    weight_dct[f'blocks.{i}.weight'],
                    weight_dct[f'blocks.{i}.bias'],
                    dilation=dilation,
                )

                # add the block's findings to the channels
                channel_arr = channel_arr + block_out_arr

            # head: GELU, then combine the 32 channels into 2 numbers per
            # window: row 0 = adjustment of the trust in the window's
            # consensus similarity, row 1 = adjustment of the tie to the
            # next window
            act_channel_arr = gelu(channel_arr)
            head_arr = conv1d(
                act_channel_arr,
                weight_dct['head.weight'],
                weight_dct['head.bias'],
            )

            # fetch trained scaling numbers: temperature (overall strength),
            # a (base trust in the consensus similarity) and b (base weight
            # of the neighbour similarity)
            temperature = np.exp(float(weight_dct['log_temp']))
            a = float(weight_dct['a'])
            b = float(weight_dct['b'])

            # evidence per window: positive = keep, negative = flip, size
            # is how confident (consensus similarity, scaled by the trust
            # from head row 0)
            window_evidence_arr = \
                temperature * cons_sim_arr * (a + head_arr[0])

            # edge evidence from each window to the next: positive = same sign,
            # negative = opposite sign (neighbour similarity, adjusted by head
            # row 1)
            edge_evidence_arr = temperature * (b * adj_arr + head_arr[1])

            # find the best flip decision per window considering all windows
            # and edges
            nonempty_flip_arr = viterbi(window_evidence_arr, edge_evidence_arr)

            # record flip decisions for non-empty windows
            flip_arr[data_idx_arr] = nonempty_flip_arr

        # flip windows marked 1
        for idx, flip in enumerate(flip_arr):
            if flip:
                pc_df.iloc[idx] *= -1

        return pc_df


class Flip:
    '''
    Flip all PC values of a chromosome or specific windows
    '''

    @staticmethod
    def flip_chrom(pc_df):
        '''
        Flip all PC values for a chromosome
        '''

        return pc_df * -1


    @staticmethod
    def flip_windows(pc_df, coord_lst):
        '''
        Flip specific windows
        '''

        # infer step size
        w_lst = list(pc_df.index)

        # iterate specified flip regions
        w_flip_lst = []
        for record in coord_lst:

            # strip accidental white spaces
            record = record.strip()

            # ignore accidental empty lines
            if record == '':
                continue

            # expand regions
            if '-' in record:
                start = int(record.split('-')[0])
                end = int(record.split('-')[1])
                w_flip_lst += [x for x in w_lst if start <= x <= end]
                #w_flip_lst += list(range(start, end, w_step)) + [end]
            else:
                w_flip_lst.append(record)

        # convert to int
        w_flip_lst = [int(x) for x in w_flip_lst]

        # check if all windows are in the index, else print error
        missing_lst = []
        for w in w_flip_lst:
            if w not in pc_df.index:
                missing_lst.append(w)
        if missing_lst:
            log.error_nl('-w/--windows: the following windows are not'
                       ' in the supplied data:'
                      f' {",".join(str(x) for x in missing_lst)}')

        # flip windows
        pc_df.loc[w_flip_lst] *= -1

        return pc_df


class Orient:
    '''
    Read PC values from a list of chromosomes and orient them consistently
    '''

    def __init__(self,
                 run_prefix,
                 run_id_lst,
                 pcs,
                 orient_pcs,
                 ):

        # input variables
        self.run_prefix = run_prefix
        self.run_id_lst = run_id_lst
        self.pcs = pcs
        self.orient_pcs = orient_pcs

        # instance variables
        self.data_lst = []

    def orient(self, epsilon):
        '''
        Read specified PC data, assess different chromosomes' relative
        orientation and orient them consistently with the most extreme samples
        pointing upwards
        '''

        # import data module
        from .data import WPCAData

        # sequentially read in data from input chromosomes
        for run_id in self.run_id_lst:
            data = WPCAData(
                f'{self.run_prefix}{run_id}',
                self.pcs,
            )
            self.data_lst.append(data)

        # apply flips per PC
        for pc in self.orient_pcs:

            # collect PC data across chromosomes
            pc_df_lst = []
            for data in self.data_lst:
                pc_df_lst.append(getattr(data, f'pc_{pc}_df'))

            # calculate per sample mean PC score per chromosome
            mean_lst = []
            for pc_df in pc_df_lst:
                mean_arr = pc_df.mean(axis=0).values.astype(float)
                mean_lst.append(np.nan_to_num(mean_arr))
            mean_arr = np.array(mean_lst)

            # scale each chromosome to magnitude 1 (epsilon prevents
            # division by zero)
            magnitude_arr = np.linalg.norm(mean_arr, axis=1)
            profile_arr = mean_arr / (magnitude_arr[:, None] + epsilon)

            # weight each sample by how extreme its PC values are across
            # chromosomes
            weight_arr = np.median(np.abs(profile_arr), axis=0) ** 2
            weighted_arr = profile_arr * weight_arr

            # rescale weighted profiles to magnitude 1 (epsilon prevents
            # division by zero)
            magnitude_arr = np.linalg.norm(weighted_arr, axis=1)
            weighted_arr = weighted_arr / (magnitude_arr[:, None] + epsilon)

            # compute PCA across chromosomes
            chrom_score_arr = np.linalg.svd(weighted_arr)[0]

            # fetch sign from chromosome PCA PC1 score
            sign_arr = np.where(chrom_score_arr[:, 0] >= 0, 1, -1)

            # flip chromosomes with sign -1 (relative orientation)
            for data, sign in zip(self.data_lst, sign_arr):
                if sign < 0:
                    data.modify_data(f'pc_{pc}_df', Flip.flip_chrom)

            # genome-wide direction: concatenate all chromosomes and
            # calculate per sample mean PC score
            genome_pc_lst = []
            for data in self.data_lst:
                pc_df = getattr(data, f'pc_{pc}_df')
                genome_pc_lst.append(pc_df.values.astype(float))
            genome_arr = np.nanmean(np.concatenate(genome_pc_lst), axis=0)
            genome_arr = np.nan_to_num(genome_arr)

            # fetch the 5 % most extreme samples genome-wide
            n_top = max(1, int(round(0.05 * len(genome_arr))))
            top_idx_arr = np.argsort(-np.abs(genome_arr))[:n_top]

            # if most extreme samples point downwards, flip all chromosomes
            if genome_arr[top_idx_arr].sum() < 0:
                for data in self.data_lst:
                    data.modify_data(f'pc_{pc}_df', Flip.flip_chrom)

    def to_files(self):
        '''
        Write the oriented data of all chromosomes to files
        '''

        for data in self.data_lst:
            data.to_files()
