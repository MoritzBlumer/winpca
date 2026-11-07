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
N_THREADS = 4              # # of threads

# winpca polarize
POL_MODE = 'global'        # polarization mode ('adaptive', 'guide_samples',
                           # 'local' or 'global')

# chromplot + genomeplot
PLOT_FMT = 'HTML'          # output  format (one of 'HTML', 'PDF', 'SVG', 'PNG)
PLOT_INTERVAL = 5          # plot only every nth value (5th if specifying 5)



## SETTINGS
#  (these can only be changed here, i.e. no CLI arguments)

# pca
PCS = [1, 2]               # PCs to output
PCANGSD_EM_EIG = 2         # sets PCAngsd '-eig' parameter (should be >= PCS)
PCANGSD_EM_ITER = 100      # max EM iterations to perform (0 --> ngsTools like)
GT_MIN_VAR_PER_W = 20      # min # of variants per window
GL_PL_MIN_VAR_PER_W = 25   # min # of variants per window
VCF_PASS_FILTER = True     # include only PASS sites (disable with --np)
SKIP_MONOMORPHIC = True    # skip invariant sites (uninformative for PCA)
GT_MEAN_IMPUTE = True      # mean-impute missing genotypes
PROC_TRAIL_TRUNC_W = True  # process truncated trailing window if present

# polarize
N_PREV_WINDOWS = 5         # # of previous windows used in adaptive polarization
N_DRAWS_LOCAL = 10         # # of draws of 64 individuals in local ML model
FLIP_THRESHOLD_LOCAL = 0.5 # threshold to trigger sign flip in local ML model
WEIGHTS_LOCAL_PATH = \
  'ml_weights/local.npz'   # local ML model weights
WEIGHTS_GLOBAL_PATH = \
  'ml_weights/global.npz'  # global ML model weights
EPSILON = 1e-12            # small value preventing division-by-zero errors

# chromplot
CHROMPLOT_W = 1200         # plot width in pixels
CHROMPLOT_H = 500          # plot height in pixels

# genomeplot
GENOMEPLOT_W = 1200       # plot width in pixels
GENOMEPLOT_H = 500        # plot height in pixels
