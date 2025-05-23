'''
Command line interface.
'''

## IMPORT CONFIG AND VERSION
from . import config
from . import __version__

## IMPORT PACKAGES
import os
import argparse

# MODULES
from modules.log import Log

## INSTANTIATE LOGGER
log = Log()



## CLASSES

class CLI:
    '''
    Command line interface and argument parser
    '''
    
    def __init__(self):

        # initiate argument parser
        self.parser = argparse.ArgumentParser(
            description='WinPCA v1.0',
            epilog='contact: lmb215@cam.ac.uk',
            formatter_class=argparse.RawDescriptionHelpFormatter
        )

        # add -v/--version
        self.parser.add_argument(
            '-v', '--version',
            action='version',
            version=f'WinPCA {__version__}'
        )

        # initiate subparsers
        self.subparsers = self.parser.add_subparsers(
            dest='winpca', required=True
        )

        # variables
        self.args = None

        # allowed variant file suffixes
        self.variant_file_suffixes = [
            '.vcf', '.vcf.gz', '.tsv', '.tsv.gz', '.beagle', '.beagle.gz',
        ]

        self.plot_file_suffixes = [
            '.html', '.pdf', '.svg', '.png',
        ]

        # define magick tab for argparse formatting
        self.tab = '\u200B \u200B \u200B \u200B '

    @staticmethod
    def shared_arguments(subparser):
        '''
        Arguments shared between all sub-commands
        '''
        # define arguments that are shared across all subparsers
        subparser.add_argument(
            dest='prefix', metavar='<PREFIX>',
            type=str, default=None,
            help='prefix for this WinPCA analysis')


    def pca(self):
        '''
        Windowed PCA using scikit-allel when working with called genotypes(GT)
        and PCAngsd when working with genotype likelihoods (GL/PL)
        '''

        # add subparser
        pca_parser = self.subparsers.add_parser(
            'pca',
            help='Perform windowed PCA on called genotypes (GT) or on'
              ' genotype likelihoods (GL/PL)'
        )

        # positional arguments
        self.shared_arguments(pca_parser)
        pca_parser.add_argument(
            dest='variant_file_path', metavar='<VARIANT_FILE>',
            type=str, default=None,
            help='path to variant file: VCF(.gz), TSV(.gz) or BEAGLE(.gz)')
        pca_parser.add_argument(
            dest='region', metavar='<REGION>',
            type=str, default=None,
            help='genomic region in format "chrom:start-end"')

        # optional arguments
        pca_parser.add_argument(
            '-s', '--samples', dest='samples', metavar='\b',
            required=False, type=str, default=None,
            help=f'{self.tab}comma-separated list of samples OR file with one'
              ' sample ID per line')
        pca_parser.add_argument(
            '-w', '--window_size', dest='w_size', metavar='\b',
            required=False, type=int, default=int(config.W_SIZE),
            help=f'window size in base pairs [{config.W_SIZE}]')
        pca_parser.add_argument(
            '-i', '--increment', dest='w_step', metavar='\b',
            required=False, type=int, default=int(config.W_STEP),
            help=f'{self.tab}step size in base pairs [{config.W_STEP}]')
        pca_parser.add_argument(
            '-m', '--min_maf', dest='min_maf', metavar='\b',
            required=False, type=float, default=float(config.MIN_MAF),
            help=f'{self.tab}minor allele frequency threshold'
             f' [{config.MIN_MAF}]')
        pca_parser.add_argument(
            '-p', '--polarize', dest='polarize', metavar='\b',
            required=False, type=str, default=str(config.POL_MODE),
            choices=['auto', 'guide_samples', 'skip'],
            help=f'{self.tab}sign polarization strategy: auto,guide_samples or'
             f' skip [{config.POL_MODE}]')
        pca_parser.add_argument(
            '-g', '--guide_samples', dest='guide_samples', metavar='\b',
            required=False, type=str, default=None,
            help='one or more (comma-separated) samples to guide PC'
              ' polarization (used only if "guide_samples" is selected for'
              ' -p/--polarize)')
        pca_parser.add_argument(
            '-v', '--var_format', dest='var_fmt', metavar='\b',
            required=False, default=str(config.VAR_FMT),
            choices=['GT', 'GL', 'PL'],
            help=f'{self.tab}variant format: GT, GL or PL (GL/PL invoke'
             f' PCAngsd) [{config.VAR_FMT}]')
        pca_parser.add_argument(
            '-t', '--threads', dest='n_threads', metavar='\b',
            required=False, type=int, default=int(config.N_THREADS),
            help=f'{self.tab}number of threads (multi-threading is only'
             f' available with PCAngsd, i.e. GL/PL) [{config.N_THREADS}]')
        pca_parser.add_argument(
            '--np', dest='no_pass_filter',
            required=False, action='store_true', default=False,
            help='disable VCF PASS filter (overrides modules/config.py)')


    def polarize(self):
        '''
        (Re)-polarize windowed PC data from a previous run.
        '''

        # add subparser
        polarize_parser = self.subparsers.add_parser(
            'polarize',
            help='(Re)-polarize windowed PC data from a previous run'
              ' (overwrites input data)'
        )

        # positional arguments
        self.shared_arguments(polarize_parser)

        # optional arguments
        polarize_parser.add_argument(
            '-c', '--principal_component', dest='pol_pc', metavar='\b',
            required=False, type=str, default=str(config.POL_PC),
            choices=[str(config.PC_A), str(config.PC_B), 'both'], 
            help='specify which PC to re-polarize (e.g. "1", "2" or "both")'
             f' [{config.POL_PC}]')
        polarize_parser.add_argument(
            '-p', '--polarize', dest='polarize', metavar='\b',
            required=False, type=str, default=str(config.POL_MODE),
            choices=['auto', 'guide_samples'],
            help=f'{self.tab}sign polarization strategy ("auto" or'
             ' "guide_samples") [{config.POL_MODE}]')
        polarize_parser.add_argument(
            '-g', '--guide_samples', dest='guide_samples', metavar='\b',
            required=False, type=str, default=None,
            help='one or more (comma-separated) samples to guide PC'
              ' polarization (used only if "guide_samples" is selected for'
              ' -p/--polarize)')

    def flip(self):
        '''
        Flip/reflect windowed PC data from a previous run (multiply values by
        -1)
        '''

        # add subparser
        flip_parser = self.subparsers.add_parser(
            'flip',
            help='Flip/reflect windowed PC data from a previous run, i.e.'
              ' multiply some or all PC scores values by -1 (overwrites'
              ' input data)'
        )

        # positional arguments
        self.shared_arguments(flip_parser)

        # optional arguments
        flip_parser.add_argument(
            '-w', '--windows', dest='flip_windows', metavar='\b',
            required=False, type=str, default=None,
            help=f'{self.tab}comma-separated list of positions (e.g. 100000)'
              ' and/or regions (e.g. 100000-250000) to flip OR file with one'
              ' position/region per line')
        flip_parser.add_argument(
            '--r', '--reflect', dest='reflect',
            required=False, action='store_true', default=None,
            help='set flag to reflect the entire chromosome, i.e. to flip all'
              ' windows (--r/--reflect is applied independently from'
              ' -w/--windows and they can be used together')
        flip_parser.add_argument(
            '-c', '--principal_component', dest='flip_pc', metavar='\b',
            required=False, type=str, default=str(config.FLIP_PC),
            choices=[str(config.PC_A), str(config.PC_B), 'both'],
            help='specify which PC to flip (e.g. "1", "2" or "both")'
             f' [{config.FLIP_PC}]')


    def chromplot(self):
        '''
        Plot PC or heterozygosity and per window stats for a specified input
        chromosome
        '''

        # add subparser
        chromplot_parser = self.subparsers.add_parser(
            'chromplot',
            help='Plot principal component or heterozygosity and per window'
              ' stats for a specified input chromosome'
        )

        # positional arguments
        self.shared_arguments(chromplot_parser)
        chromplot_parser.add_argument(
            dest='region', metavar='<REGION>',
            help='genomic region in format "chrom:start-end"')

        # optional arguments
        chromplot_parser.add_argument(
            '-p', '--plot_variable', dest='plot_var', metavar='\b',
            required=False, type=str, default=str(config.PC_A),
            choices=[str(config.PC_A), str(config.PC_B), 'het'], 
            help='specify what to plot, e.g. "1" for PC 1 or "het" for SNP'
              ' heterozygosity (default: first specified PC)'
             f' [{config.PC_A}]')
        chromplot_parser.add_argument(
            '-m', '--metadata', dest='metadata_path', metavar='\b',
            required=False, type=str, default=None,
            help=f'{self.tab}path to metadata TSV; first column must contain'
              ' sample IDs, all other columns are used to annotate data in HTML'
              ' plot')
        chromplot_parser.add_argument(
            '-g', '--groups', dest='color_by', metavar='\b',
            required=False, type=str, default=None,
            help=f'{self.tab}metadata column for color-grouping (requires'
              ' -m/--metadata)')
        chromplot_parser.add_argument(
            '-c', '--colors', dest='hex_codes', metavar='\b',
            required=False, type=str, default=None,
            help=f'{self.tab}HEX codes (drop "#") for color groups; check'
              ' documentation for formatting instructions (requires'
              ' -g/--groups)')
        chromplot_parser.add_argument(
            '-i', '--interval', dest='interval', metavar='\b',
            required=False, type=int, default=int(config.PLOT_INTERVAL),
            help=f'{self.tab}plot every nth window (e.g. "5" --> every 5th'
             f' window) [{config.PLOT_INTERVAL}]')
        chromplot_parser.add_argument(
            '-f', '--format', dest='plot_fmt', metavar='\b',
            required=False, type=str, default=str(config.PLOT_FMT),
            help=f'{self.tab}output plot file format ("HTML", "PDF", "SVG" or'
             f' "PNG"); specify a single one or a comma-separated list'
              ' [{config.PLOT_FMT}]')
        chromplot_parser.add_argument(
            '--n', '--numeric', dest='numeric',
            required=False, action='store_true', default=False,
            help='set flag to apply a continuous color scale (requires'
              ' numerical data for -g/--groups)')
        chromplot_parser.add_argument(
            '--r', '--reverse', dest='reverse',
            required=False, action='store_true', default=False,
            help='set flag to reverse the plotting order (requires'
              ' -g/--groups)')


    def genomeplot(self):
        '''
        Plot PC or heterozygosity for the specified input chromosomes
        '''

        # add subparser
        genomeplot_parser = self.subparsers.add_parser(
            'genomeplot',
            help='Plot PC or heterozygosity for the specified input'
              ' chromosomes'
        )

        # positional arguments
        genomeplot_parser.add_argument(
            dest='run_prefix', metavar='<RUN_PREFIX>',
            type=str, default=None,
            help='prefix shared by all chromosomes runs to include in'
              ' genome-wide plot')
        genomeplot_parser.add_argument(
            dest='run_ids', metavar='<RUN_IDS>',
            type=str, default=None,
            help=f'comma-separated list of run IDs to include, format: e.g.'
              ' {prefix}.{run_id}.pc_1.tsv.gz; also used to determine plotting;'
              ' order')

        # positional arguments
        genomeplot_parser.add_argument(
            '-p', '--plot_variable', dest='plot_var', metavar='\b',
            required=False, type=str, default=str(config.PC_A),
            choices=[str(config.PC_A), str(config.PC_B), 'het'],
            help='specify what to plot, e.g. "1" for PC 1 or "het" for SNP'
              ' heterozygosity (default: first specified PC)'
             f' [{config.PC_A}]')
        genomeplot_parser.add_argument(
            '-m', '--metadata', dest='metadata_path', metavar='\b',
            required=False, type=str, default=None,
            help=f'{self.tab}path to metadata TSV; first column must contain'
              ' sample IDs, all other columns are used to annotate data in HTML'
              ' plot')
        genomeplot_parser.add_argument(
            '-g', '--groups', dest='color_by', metavar='\b',
            required=False, type=str, default=None,
            help=f'{self.tab}metadata column for color-grouping (requires'
              ' -m/--metadata)')
        genomeplot_parser.add_argument(
            '-c', '--colors', dest='hex_codes', metavar='\b',
            required=False, type=str, default=None,
            help=f'{self.tab}HEX codes (drop "#") for color groups; check'
              ' documentation for formatting instructions (requires'
              ' -g/--groups)')
        genomeplot_parser.add_argument(
            '-i', '--interval', dest='interval', metavar='\b',
            required=False, type=int, default=int(config.PLOT_INTERVAL),
            help=f'{self.tab}plot every nth window (e.g. "5" --> every 5th'
             f' window) [{config.PLOT_INTERVAL}]')
        genomeplot_parser.add_argument(
            '-f', '--format', dest='plot_fmt', metavar='\b',
            required=False, type=str, default=str(config.PLOT_FMT),
            help=f'{self.tab}output plot file format ("HTML", "PDF", "SVG" or'
             f' "PNG"); specify a single one or a comma-separated list'
              ' [{config.PLOT_FMT}]')
        genomeplot_parser.add_argument(
            '--n', '--numeric', dest='numeric',
            required=False, action='store_true', default=False,
            help='set flag to apply a continuous color scale (requires'
              ' numerical data for -g/--groups)')
        genomeplot_parser.add_argument(
            '--r', '--reverse', dest='reverse',
            required=False, action='store_true', default=False,
            help='set flag to reverse the plotting order (requires'
              ' -g/--groups)')


    def parse_args(self):
        '''
        Parse and sanity check command-line arguments
        '''

        # parse arguments
        self.args = vars(
            self.parser.parse_args()
        )

        # handle interdependent options
        if self.args.get('polarize') == 'guide_samples' \
            and not self.args.get('guide_samples'):
            log.error_nl(
                 '-p/--polarize: please provide at least one guide sample via'
                 ' -g/--guide_samples'
            )
        if self.args.get('winpca') == 'flip' \
            and not self.args.get('reflect') \
            and not self.args.get('flip_windows'):
            log.error_nl(
                 'One of --r/--reflect and -w/--windows must be set'
            )
        if self.args.get('color_by') \
            and not self.args.get('metadata_path'):
            log.error_nl(
                 '-g/--group: -m/--metadata is required to infer groups'
            )
        if self.args.get('hex_codes') \
            and not self.args.get('color_by'):
            log.error_nl(
                 '-c/--colors: -g/--groups is required to set specify'
                 ' group-specificcolors'
            )


        # variant_file_path
        if self.args.get('variant_file_path'):
            path = self.args['variant_file_path']
            if path.lower().endswith('.gz'):
                suffix = '.' + '.'.join(path.split('.')[-2:])
            else:
                suffix = '.' + '.'.join(path.split('.')[-1:])
            if suffix.lower() not in self.variant_file_suffixes:
                log.error_nl(
                    f'VARIANT_FILE: {self.args["variant_file_path"]} has an'
                    f' unexpected suffix ({suffix} is not a supported file'
                    f' format or compression)'
                )
            else:
                if not os.path.exists(path):
                    log.error_nl(
                        f'VARIANT_FILE: {self.args["variant_file_path"]}:'
                         ' file does not exist'
                    )

        # region --> chrom, start, end
        if self.args.get('region'):
            if ':' not in self.args['region'] or '-' not in self.args['region']:
                log.error_nl(
                    f'REGION: {self.args["region"]} is formatted incorrectly,'
                     ' please make sure to specify genomic coordinates:'
                     ' correctly chrom:start-end (start and end must always be'
                     ' specified)'
                )
            else:
                self.args['chrom'] = self.args['region'].split(':')[0]
                self.args['start'] = \
                    int(self.args['region'].split(':')[1].split('-')[0])
                self.args['end'] = \
                    int(self.args['region'].split(':')[1].split('-')[1])

        # samples --> sample_lst
        if self.args.get('samples'):
            if ',' in self.args['samples']:
                self.args['sample_lst'] = self.args['samples'].split(',')
            else:
                if not os.path.exists(self.args['samples']):
                    log.error_nl(
                        f'-s/--samples: {self.args["samples"]}: file does not'
                         ' exist'
                    )
                else:
                    self.args['sample_lst'] = []
                    with open(self.args['samples'], 'r') as sample_file:
                        for line in sample_file:
                            self.args['sample_lst'].append(
                                line.strip().split('\t')[0]
                            )
            if len(set(self.args['sample_lst'])) < len(self.args['sample_lst']):
                seen, dups = set(), set()
                for id in self.args['sample_lst']:
                    dups.add(id) if id in seen else seen.add(id)
                log.error_nl(
                    f'-s/--samples: duplicate IDs found: {", ".join(dups)}'
                )
        else: 
            self.args['sample_lst'] = None

        # guide_samples --> guide_sample_lst
        if self.args.get('guide_samples'):
            self.args['guide_sample_lst'] = \
                self.args['guide_samples'].split(',')
        else:
            self.args['guide_sample_lst'] = None

        # flip_windows --> flip_window_lst
        if self.args.get('flip_windows'):
            if ',' in self.args['flip_windows']:
                self.args['flip_window_lst'] = \
                    self.args['flip_windows'].split(',')
            else:
                if os.path.exists(self.args['flip_windows']):
                    self.args['flip_window_lst'] = []
                    with open(self.args['flip_windows'], 'r') as flip_win_file:
                        for line in flip_win_file:
                            self.args['flip_window_lst'].append(
                                line.strip().split('\t')[0]
                            )
                else:
                    self.args['flip_window_lst'] = [self.args['flip_windows']]
        else:
            self.args['flip_window_lst'] = None

        # plot_var
        if self.args.get('plot_var'):
            if self.args['plot_var'] == 'het':
                self.args['plot_var'] = 'hetp'

        if self.args.get('metadata_path'):
            if not os.path.exists(self.args['metadata_path']):
                log.error_nl(
                    f'-m/--metadata: {self.args["metadata_path"]}: file does'
                     ' not exist'
                )

        # hex codes --> hex_code_dct
        if self.args.get('hex_codes'):       
            self.args['hex_code_dct'] = {}
            for i in self.args['hex_codes'].split(','):
                group, hex_code = i.split(':')[0], i.split(':')[1]
                self.args['hex_code_dct'][group] = f'#{hex_code}'
        else:
            self.args['hex_code_dct']= None

        # plot_fmt
        if self.args.get('plot_fmt'):
            if ',' in self.args['plot_fmt']:
                self.args['plot_fmt_lst'] = \
                    [x.lower() for x in self.args['plot_fmt'].split(',')]
            else:
                self.args['plot_fmt_lst'] = [self.args['plot_fmt'].lower()]
            for fmt in self.args['plot_fmt_lst']:
                if f'.{fmt}' not in self.plot_file_suffixes:
                    log.error_nl(
                        f'-f/--format: {fmt} is not supported as output format'
                    )

        if self.args.get('run_ids'):
            if ',' not in self.args['run_ids']:
                log.error_nl(
                    'RUN_IDS: provide at least two sequences as a'
                    ' comma-separated list'
                )
            else:
                self.args['run_id_lst'] = self.args['run_ids'].split(',')
            for r_id in self.args['run_id_lst']:
                path = \
                    f'{self.args["run_prefix"]}{r_id}.pc_{config.PC_A}.tsv.gz'
                if not os.path.exists(path):
                    log.error_nl(
                        f'RUN_PREFIX: {self.args["run_prefix"]}{r_id} does'
                         ' not exist or is incomplete'
                    )

        # add config settings
        self.args['skip_monomorphic'] = bool(config.SKIP_MONOMORPHIC)
        self.args['gt_min_var_per_w'] = int(config.GT_MIN_VAR_PER_W)
        self.args['gl_pl_min_var_per_w'] = int(config.GL_PL_MIN_VAR_PER_W)
        self.args['n_prev_windows'] = int(config.N_PREV_WINDOWS)
        self.args['gt_mean_impute'] = bool(config.GT_MEAN_IMPUTE)
        self.args['n_pcs'] = str(config.N_PCS)
        self.args['pc_a'] = str(config.PC_A)
        self.args['pc_b'] = str(config.PC_B)
        self.args['chromplot_w'] = int(config.CHROMPLOT_W)
        self.args['chromplot_h'] = int(config.CHROMPLOT_H)
        self.args['genomeplot_w'] = int(config.GENOMEPLOT_W)
        self.args['genomeplot_h'] = int(config.GENOMEPLOT_H)

        # add default values from config if unset
        if not self.args.get('pol_pc'):
            self.args['pol_pc'] = str(config.POL_PC)
        if self.args.get('no_pass_filter') == True:
            self.args['vcf_pass_filter'] = False
        else:
            self.args['vcf_pass_filter'] = bool(config.VCF_PASS_FILTER)

        # check if N_PCS is at larger than PC_1 and PC_2 setting or 0 if GL/PL
        if self.args.get('var_fmt'):
            if not ((int(config.PC_A) < int(config.PC_B) <= int(config.N_PCS)) \
                or int(config.N_PCS) == 0 and (self.args['var_fmt'] in ['GL', 'PL'])):
                            log.error_nl(
                                'modules/config.py: N_PCS must be >= PC_A/PC_B ' \
                                ' setting'
                            )

        # check if specified PC exists in data if requested
        if self.args.get('pol_pc'):
            if self.args['pol_pc'] \
                not in [self.args['pc_a'], self.args['pc_b'], 'both']:
                log.error_nl(
                    f'-c/--principal_component: no data found for PC'
                    f' {self.args["pol_pc"]}'
                )
