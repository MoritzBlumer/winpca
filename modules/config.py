'''
Default settings.
'''

# pca + pcangsd
SKIP_MONOMORPHIC=False
MIN_VAR_PER_W = 25
MIN_MAF = 0.01
W_SIZE = 1000000
W_STEP = 10000

# pcangsd
N_THREADS = 2

# polarize
N_PREV_WINDOWS = 5
POL_MODE = 'auto'
POL_PC = 'both'

# flip
FLIP_PC = '1'

# chromplot + genomeplot
PLOT_FMT = 'html'
PLOT_INTERVAL = None
CHROMPLOT_W = 1200
CHROMPLOT_H = 400
GENOMEPLOT_W = 1200
GENOMEPLOT_H = 300
