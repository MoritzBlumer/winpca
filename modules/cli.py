'''
Command line interface.
'''

## IMPORT CONFIG AND VERSION
from . import config
from . import __version__

## IMPORT PACKAGES
import os
import argparse


## CLASSES

class CLI:
    '''
    Command line interface and argument parser.
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
        self.args_dct = None

        # allowed variant file suffixes
        self.variant_file_suffixes = [
            '.vcf', '.vcf.gz', '.tsv', '.tsv.gz', '.beagle', '.beagle.gz',
        ]

        self.plot_file_suffixes = [
            '.html', '.pdf', '.svg', '.png',
        ]

    @staticmethod
    def shared_arguments(subparser):
        '''
        Arguments shared between all sub-commands.
        '''
        # define arguments that are shared across all subparsers
        subparser.add_argument(
            dest='prefix', metavar='<PREFIX>',
            help='Prefix for this run, i.e. prefix for all results generated in'
            ' this WinPCA analysis.')


    def pca(self):
        '''
        Windowed PCA using scikit-allel when working with callled genotypes(GT)'
        ' and PCAngsd when working with genotype likelihoods (GL/PL).
        '''

        # add subparser
        pca_parser = self.subparsers.add_parser(
            'pca', help='Perform windowed PCA on called genotypes (GT) or on'
            ' genotype likelihoods (GL/PL).'
        )

        # positional arguments
        pca_parser.add_argument(
            dest='variant_file_path', metavar='<VARIANT_FILE>', help='Path to'
            ' variant file (optionally gzipped VCF, TSV or BEAGLE; see'
            ' documentation for input file specifications).')
        pca_parser.add_argument(
            dest='region', metavar='<REGION>', help='Genomic region in format'
            ' "chrom:start-end".')
        self.shared_arguments(pca_parser)
        
        # optional arguments
        pca_parser.add_argument(
            '-s', '--samples', dest='samples', required=False, metavar='\b',
            help='''Comma-separated list of samples to include or file with one'
            ' sample per line.''')
        pca_parser.add_argument(
            '-w', '--window_size', dest='w_size', required=False, type=int,
            default=config.W_SIZE, metavar='\b', help='Window size in base'
            f' pairs (bp) [default: {config.W_SIZE}].')
        pca_parser.add_argument(
            '-i', '--increment', dest='w_step', required=False, type=int,
            default=config.W_STEP, metavar='\b', help='Step size in base pairs'
            f' (bp). [default: {config.W_STEP}].')
        pca_parser.add_argument(
            '-m', '--min_maf', dest='min_maf', required=False, type=float,
            default=config.MIN_MAF, metavar='\b', help='Minor allele frequency'
            f' threshold [default: {config.MIN_MAF}].')
        pca_parser.add_argument(
            '-p', '--polarize', dest='polarize', required=False,
            default=config.POL_MODE, choices=['auto', 'guide_samples', 'skip'],
            metavar='\b', help='Sign polarization strategy'
            ' ("auto", "guide_samples" or "skip") [default:'
            f' "{config.POL_MODE}"].')
        pca_parser.add_argument(
            '-g', '--guide_samples', dest='guide_samples', required=False,
            metavar='\b', help='Applies only if "guide_samples" is selected for'
            ' -p/--polarize: One or more (-> comma-separated list) samples to'
            ' guide PC sign polarization.')
        pca_parser.add_argument(
            '-v', '--var_format', dest='var_fmt', required=False, 
            default=config.VAR_FMT, choices=['GT', 'GL', 'PL'], metavar='\b', 
            help='Variant format ("GT", "GL", "PL") [default:' 
            f' {config.VAR_FMT}]. GL/PL invoke PCAngsd usage.')
        pca_parser.add_argument(
            '-t', '--threads', dest='threads', required=False, type=int,
            default=config.N_THREADS, metavar='\b', help='Number of threads'
            f' [default: {config.N_THREADS}]. Multithreading is only used'
            ' when invoking PCAngsd, i.e. when inputting GL/PL variants.')


    def polarize(self):
        '''
        (Re)-polarize windowed PC data from a previous run.
        '''

        # add subparser
        polarize_parser = self.subparsers.add_parser(
            'polarize', help='(Re)-polarize windowed PC data from a previous'
            ' run. Overwrites input data.'
        )

        # positional arguments
        self.shared_arguments(polarize_parser)

        # optional arguments
        polarize_parser.add_argument(
            '-c', '--principal_component', dest='pol_pc', required=False,
            choices=['1', '2', 'both'], default=1, metavar='\b', help='Specify'
            ' which PC to re-polarize ("1", "2" or "both") [default: 1].')
        polarize_parser.add_argument(
            '-p', '--polarize', dest='polarize', required=False,
            default=config.POL_MODE, choices=['auto', 'guide_samples'],
            metavar='\b', help='Sign polarization strategy ("auto" or'
            f' "guide_samples") [default: {config.POL_MODE}].')
        polarize_parser.add_argument(
            '-g', '--guide_samples', dest='guide_samples', required=False,
            metavar='\b', help='Applies only if "guide_samples" is selected for'
            ' -p/--polarize: One or more (-> comma-separated list) samples to'
            ' guide PC sign polarization.')


    def flip(self):
        '''
        Flip/reflect windowed PC data from a previous run (multiply values by
        -1).
        '''

        # add subparser
        flip_parser = self.subparsers.add_parser(
            'flip', help='Flip/reflect windowed PC data from a previous run'
            ' (multiply values by -1). Overwrites input data.'
        )

        # positional arguments
        self.shared_arguments(flip_parser)

        # optional arguments
        flip_parser.add_argument(
            '-w', '--windows', dest='flip_windows', required=False,
            metavar='\b', help='Comma-separated list of positions'
            ' (e.g. 100000) or regions (100000-250000) to flip or file with'
            ' one position/region per line.')
        flip_parser.add_argument(
            '--r', '--reflect', dest='reflect', required=False,
            action='store_true', help='Set flag to reflect the entire'
            ' chromosome, i.e. to flip all windows. --r/--reflect is'
            ' applied independently from -w/--windows and both can be'
            ' combined.')
        flip_parser.add_argument(
            '-c', '--principal_component', dest='flip_pc', required=False,
            choices=['1', '2', 'both'], default=config.FLIP_PC, metavar='\b',
            help='Specify which PC to flip ("1", "2" or "both")'
            f' [default: {config.FLIP_PC}].')


    def chromplot(self):
        '''
        Plot PC1, PC2, heterozygosity and per window stats for a specified input
        chromosome.
        '''

        # add subparser
        chromplot_parser = self.subparsers.add_parser(
            'chromplot', help='Plot PC1, PC2, heterozygosity and per window'
            ' stats for a specified input chromosome.'
        )

        # positional arguments
        self.shared_arguments(chromplot_parser)
        chromplot_parser.add_argument(
            dest='region', metavar='<REGION>', help='Genomic region in format'
            ' "chrom:start-end".')
        
        # optional arguments
        chromplot_parser.add_argument(
            '-p', '--plot_variable', dest='plot_var', required=False,
            choices=['PC1', 'PC2', 'het'], default=config.PLOT_VAR,
            metavar='\b', help='Specify which values to plot ("PC1", "PC2" or'
            f' "het") [default: {config.PLOT_VAR}].')
        chromplot_parser.add_argument(
            '-m', '--metadata', dest='metadata_path', required=False,
            metavar='\b', help='Path to metadata TSV where first column are'
            ' sample IDs. Additional columns will be used to annotate data in'
            ' HTML plot.')
        chromplot_parser.add_argument(
            '-g', '--groups', dest='color_by', required=False, metavar='\b',
            help='Metadata column for color-grouping. Requires -m/--metadata.')
        chromplot_parser.add_argument(
            '-c', '--colors', dest='hex_codes', required=False, metavar='\b',
            help='HEX codes (drop "#") for color groups. Check documentation'
            ' for formatting instructions. Requires -g/--groups.')
        chromplot_parser.add_argument(
            '-i', '--interval', dest='interval', required=False, type=int,
            default=None, metavar='\b', help='If set, only plot values for every'
            ' nth window (10 --> 10th) [default: no interval].')
        chromplot_parser.add_argument(
            '-f', '--format', dest='plot_fmt', required=False, metavar='\b',
            default=config.PLOT_FMT, help='Output plot file format ("HTML",'
            ' "PDF", "SVG" or "PNG") [default: {config.PLOT_FMT}].')


    def genomeplot(self):
        '''
        Plot PC1, PC2 or heterozygosity for the specified input chromosomes.
        '''

        # add subparser
        genomeplot_parser = self.subparsers.add_parser(
            'genomeplot', help='Plot PC1, PC2 or heterozygosity for the'
            ' specified input chromosomes.'
        )

        # positional arguments
        genomeplot_parser.add_argument(
            dest='run_prefix', metavar='<RUN_PREFIX>', help='Prefix shared by'
            ' all chromosomes runs to include in genome-wide plot.')
        genomeplot_parser.add_argument(
            dest='run_ids', metavar='<RUN_IDS>', help='Comma-separated list of'
            ' run IDs to include, format: e.g. {prefix}.{run_id}.pc_1.tsv.gz.'
            ' Also used to determine plotting order.')
        
        # positional arguments
        genomeplot_parser.add_argument(
            '-p', '--plot_variable', dest='plot_var', required=False,
            choices=['PC1', 'PC2', 'het'], default=config.PLOT_VAR,
            metavar='\b', help='Specify which values to plot ("PC1", "PC2" or'
            f' "het") [default: {config.PLOT_VAR}].')
        genomeplot_parser.add_argument(
            '-m', '--metadata', dest='metadata_path', required=False,
            metavar='\b', help='Path to metadata TSV where first column are'
            ' sample IDs. Additional columns will be used to annotate data in'
            ' HTML plot.')
        genomeplot_parser.add_argument(
            '-g', '--groups', dest='color_by', required=False, metavar='\b',
            help='Metadata column for color-grouping. Requires -m/--metadata.')
        genomeplot_parser.add_argument(
            '-c', '--colors', dest='hex_codes', required=False, metavar='\b',
            help='HEX codes (drop "#") for color groups. Check documentation'
            ' for formatting instructions. Requires -g/--groups.')
        genomeplot_parser.add_argument(
            '-i', '--interval', dest='interval', required=False, type=int,
            default=config.PLOT_INTERVAL, metavar='\b', help='If set, only plot'
            ' values for every nth window (10 --> 10th)'
            f' [default: {config.PLOT_INTERVAL}.')
        genomeplot_parser.add_argument(
            '-f', '--format', dest='plot_fmt', required=False, metavar='\b',
            default=config.PLOT_FMT, help='Output plot file format ("HTML",'
            ' "PDF", "SVG" or "PNG") [default: {config.PLOT_FMT}].')


    def parse_args(self):
        '''
        Parse and sanity check command-line arguments.
        '''

        # parse arguments
        args = self.parser.parse_args()

        # handle interdependent options
        if hasattr(args, 'polarize') and args.polarize == 'guide_samples':
            if not args.guide_samples:
                self.parser.error(
                    '-g/--guide_samples is required when chosing'
                    ' "guide_samples" as polarization strategy (-p)')
        if hasattr(args, 'color_by') and not hasattr(args, 'metadata_path'):
            self.parser.error(
                '-m/--metadata is required to infer -g/--groups.')
        if hasattr(args, 'hex_codes') and not hasattr(args, 'color_by'):
            self.parser.error(
                '-g/--groups is required if -c/--colors is set.')
        if args.winpca == 'flip' \
            and not hasattr(args, 'reflect') \
            and not hasattr(args, 'flip_windows'):
                self.parser.error(
                    'One of --r/--reflect and -w,--windows must be set.')

        # check formatting
        if hasattr(args, 'variant_file_path'):
            if args.variant_file_path.lower().endswith('.gz'):
                suffix = '.' + '.'.join(args.variant_file_path.split('.')[-2:])
            else:
                suffix = '.' + '.'.join(args.variant_file_path.split('.')[-1:])
            if suffix.lower() not in self.variant_file_suffixes:
                self.parser.error(
                    f'{args.variant_file_path} has an unexpected file format'
                    + f' based on its suffix "{suffix}")')
        if hasattr(args, 'region'):
            if ':' not in args.region or '-' not in args.region:
                self.parser.error(
                    f'{args.region} is formatted incorrectly. Please make sure'
                    + ' to specify genomic coordinates correctly:'
                    + ' chrom:start-end (start and end must always be'
                    + ' specified).' )
        if hasattr(args, 'run_ids'):
            if not ',' in args.run_ids:
                self.parser.error(
                    'Please provide at least two sequences as a comma-separated'
                    ' list.' )

        # decompose/preprocess arguments
        if hasattr(args, 'region'):
            chrom = args.region.split(':')[0]
            start = int(args.region.split(':')[1].split('-')[0])
            end = int(args.region.split(':')[1].split('-')[1])
        if hasattr(args, 'samples') and not args.samples is None:
            if ',' in args.samples:
                sample_lst = args.samples.split(',')
            else:
                if not os.path.exists(args.samples):
                    self.parser.error(
                        args.samples + ': file does not exist.')
                else:
                    sample_lst = []
                    with open(args.samples, 'r') as sample_file:
                        for line in sample_file:
                            sample_lst.append(line.strip().split('\t')[0])
        if hasattr(args, 'flip_windows') and not args.flip_windows is None:
            if ',' in args.flip_windows:
                flip_window_lst = args.flip_windows.split(',')
            else:
                if os.path.exists(args.flip_windows):
                    flip_window_lst = []
                    with open(args.flip_windows, 'r') as flip_windows_file:
                        for line in flip_windows_file:
                            flip_window_lst.append(line.strip().split('\t')[0])
                else:
                    flip_window_lst = [args.flip_windows]
        if hasattr(args, 'guide_samples'):
            if args.guide_samples:
                guide_sample_lst = args.guide_samples.split(',')
            else:
                guide_sample_lst = None
        if hasattr(args, 'plot_var'):
            if args.plot_var == 'PC1': args.plot_var = 'pc_1'
            if args.plot_var == 'PC2': args.plot_var = 'pc_2'
            if args.plot_var == 'het': args.plot_var = 'hetp'
        if hasattr(args, 'hex_codes'):
            if args.hex_codes:
                hex_code_dct = {}
                for i in args.hex_codes.split(','):
                    group, hex_code = i.split(':')[0], i.split(':')[1]
                    hex_code_dct[group] = '#' + hex_code
            else:
                hex_code_dct = None
        if hasattr(args, 'plot_fmt'):
            if ',' in args.plot_fmt:
                plot_fmt_lst = args.plot_fmt.split(',')
            else:
                plot_fmt_lst = [args.plot_fmt]
            plot_fmt_lst = [x.lower() for x in plot_fmt_lst]
        if hasattr(args, 'run_ids'):
            run_id_lst = args.run_ids.split(',')

        # further checks
        if hasattr(args, 'plot_fmt'):
            for fmt in plot_fmt_lst:
                if '.' + fmt not in self.plot_file_suffixes:
                    self.parser.error(
                        f'"{fmt}" is not supported as output format.')

        # check if files exist
        if hasattr(args, 'variant_file_path'):
            if not os.path.exists(args.variant_file_path):
                self.parser.error(
                    args.variant_file_path + ': file does not exist.')
        if hasattr(args, 'metadata_path') and args.metadata_path:
            if not os.path.exists(args.metadata_path):
                self.parser.error(
                    args.metadata_path + ': file does not exist.')
        if hasattr(args, 'run_ids'):
            for r_id in run_id_lst:
                if not os.path.exists(f'{args.run_prefix}{r_id}.pc_1.tsv.gz'):
                    self.parser.error(
                        f'{args.run_prefix}{r_id}: run does not exist or is'
                        ' incomplete.')

        # convert args to dict and add to class as instance variable
        self.args_dct = vars(args)

        # add derived/processed arguments
        if hasattr(args, 'region'):
            self.args_dct['chrom'] = chrom
            self.args_dct['start'] = start
            self.args_dct['end'] = end
        if hasattr(args, 'samples') and not args.samples is None:
            self.args_dct['sample_lst'] = sample_lst
        if hasattr(args, 'guide_samples'):
            self.args_dct['guide_sample_lst'] = guide_sample_lst
        if hasattr(args, 'hex_codes'):
            self.args_dct['hex_code_dct'] = hex_code_dct
        if hasattr(args, 'plot_fmt'):
            self.args_dct['plot_fmt_lst'] = plot_fmt_lst
        if hasattr(args, 'run_ids'):
            self.args_dct['run_id_lst'] = run_id_lst
        if hasattr(args, 'flip_windows') and not args.flip_windows is None:
            self.args_dct['flip_window_lst'] = flip_window_lst

        # add in settings from config
        self.args_dct['skip_monomorphic'] = config.SKIP_MONOMORPHIC
        self.args_dct['min_var_per_w'] = config.MIN_VAR_PER_W
        self.args_dct['vcf_pass_filter'] = config.VCF_PASS_FILTER      
        self.args_dct['n_prev_windows'] = config.N_PREV_WINDOWS

        # add default values from config if unset
        if not 'pol_pc' in self.args_dct:
            self.args_dct['pol_pc'] = config.POL_PC
        if not 'flip_pc' in self.args_dct:
            self.args_dct['flip_pc'] = config.FLIP_PC
        if not 'chromplot_w' in self.args_dct:
            self.args_dct['chromplot_w'] = config.CHROMPLOT_W
        if not 'chromplot_h' in self.args_dct:
            self.args_dct['chromplot_h'] = config.CHROMPLOT_H
        if not 'genomeplot_w' in self.args_dct:
            self.args_dct['genomeplot_w'] = config.GENOMEPLOT_W
        if not 'genomeplot_h' in self.args_dct:
            self.args_dct['genomeplot_h'] = config.GENOMEPLOT_H
        if not 'n_threads' in self.args_dct:                                    # WHY IS THIS NOT AUTOMATICALLY IN THE args_dct ?
            self.args_dct['n_threads'] = config.N_THREADS
        
        # add sample_lst if not unset
        if not 'sample_lst' in self.args_dct:                                    # WHY IS THIS NOT AUTOMATICALLY IN THE args_dct ?
            self.args_dct['sample_lst'] = None


