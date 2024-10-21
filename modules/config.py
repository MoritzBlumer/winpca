'''
Configuration.
'''


## DEFAULTS
#  (overridden by CLI arguments, see documentation)

# winpca pca
VAR_FMT = 'GT'             # variant type (one of 'GT', 'GL', 'PL')
MIN_MAF = 0.01             # minor allele frequency threshold per window
W_SIZE = 1000000           # window size in bp
W_STEP = 10000             # step size in bp
N_THREADS = 2              # # of threads – only affects PCAngsd
POL_PC = 'both'            # which PC to polarize (one of '1', '2', 'both')

# winpca polarize
POL_MODE = 'auto'          # polarization mode ('auto' or 'guide_samples')

# winpca flip
FLIP_PC = '1'              # which PC to modify (one of '1', '2', 'both')

# chromplot + genomeplot
PLOT_VAR = 'PC1'           # values to plot (one of 'PC1', 'PC2', 'het')
PLOT_FMT = 'html'          # output  format (one of 'HTML', 'PDF', 'SVG', 'PNG)
PLOT_INTERVAL = None       # plot only every nth value (10th if specifying 10)



## SETTINGS
#  (these can only be changed here, i.e. no CLI arguments)

# pca
MIN_VAR_PER_W = 25        # min # of variants per window, otherwise empty output
GT_MIN_VAR_PER_W = 25     # min # of variants per window, otherwise empty output
GL_PL_MIN_VAR_PER_W = 100 # min # of variants per window, otherwise empty output
SKIP_MONOMORPHIC = True   # skip invariant sites
VCF_PASS_FILTER = True    # include only sites where FILTER is set to PASS
GT_MEAN_IMPUTE = True     # mean-impute missing genotypes (check documentation)

# polarize
N_PREV_WINDOWS = 5        # # of previous windows to consider for polarization

# chromplot
CHROMPLOT_W = 1200        # plot width in pixels
CHROMPLOT_H = 400         # plot height in pixels

# genomeplot
GENOMEPLOT_W = 1200       # plot width in pixels
GENOMEPLOT_H = 500        # plot height in pixels
