'''
Command line interface.
'''

# IMPORT PACKAGES
import argparse
import os

# IMPORT MODULES
from modules import config

print(config.min_maf)
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
        )

        # initiate subparsers
        self.subparsers = self.parser.add_subparsers(
            dest='winpca', required=True
        )

        # variables
        self.args_dct = None

        # allowed variant file suffixes
        self.variant_file_suffixes = [
            '.vcf', '.vcf.gz', '.tsv', '.tsv.gz', 'beagle', 'beagle.gz',
        ]

        self.plot_file_suffixes = [
            '.html', '.pdf', '.svg', '.png',
        ]

    def shared_arguments(self, subparser):
        '''
        Arguments shared between all sub-commands.
        '''
        # define arguments that are shared across all subparsers
        subparser.add_argument(
            '-n', '--name', dest='prefix', required=True, metavar='\b', 
            help='Name of the run, i.e. prefix for all results generated in'
            ' this WinPCA analysis.')  


    def pca(self):
        '''
        Windowed PCA on called genotypes (GT) with scikit-allel.
        '''

        # add subparser
        pca_parser = self.subparsers.add_parser(
            'pca', help='Perform windowed PCA on called genotypes (GT),'
            ' accepts VCF or TSV as input.'
        )

        # add shared arguments
        self.shared_arguments(pca_parser)

        # define subparser-specific arguments
        variant_file = pca_parser.add_mutually_exclusive_group(required=True)
        variant_file.add_argument(
            '-v', '--vcf', dest='variant_file_path', metavar='\b', help='Path to'
            ' (optionally gzipped) VCF file with GT field.')
        variant_file.add_argument(
            '-t', '--tsv', dest='variant_file_path', metavar='\b', help='Path to'
            ' (optionally gzipped) TSV file with genotypes (see README for'
            ' specifications).')
        pca_parser.add_argument(
            '-r', '--region', dest='region', required=True, metavar='\b', 
            help='Genomic region in format "chrom:start-end".')
        pca_parser.add_argument(
            '-s', '--samples', dest='samples', required=False, metavar='\b',
            help='Comma-separated list of samples to include or file with one'
            ' sample per line.')
        pca_parser.add_argument(
            '-w', '--window_size', dest='w_size', required=False, type=int, 
            default=config.w_size, metavar='\b', help='Window size in base pairs'
            f' (bp) [default: {config.w_size}].')
        pca_parser.add_argument(
            '-i', '--increment', dest='w_step', required=False, type=int,
            default=config.w_step, metavar='\b', help='Step size in base pairs'
            f' (bp). [default: {config.w_step}].')
        pca_parser.add_argument(
            '-m', '--min_maf', dest='min_maf', required=False, type=float, 
            default=config.min_maf, metavar='\b', help='Minor allele frequency'
            f' threshold [default: {config.min_maf}].')
        pca_parser.add_argument(
            '-p', '--polarize', dest='polarize', required=False, 
            default=config.pol_mode, choices=['auto', 'guide_samples', 'skip'], 
            metavar='\b', help='Sign polarization strategy'
            ' ("auto"/"guide_samples" or "skip") [default:'
            f' "{config.pol_mode}"].')
        pca_parser.add_argument(
            '-g', '--guide_samples', dest='guide_samples', required=False, 
            metavar='\b', help='Applies only if "guide_samples" is selected for'
            ' -p/--polarize: One or more (-> comma-separated list) samples to'
            ' guide PC sign polarization.')

    def pcangsd(self):
        '''
        Windowed PCA on genotype likelihoods (GL or PL) with scikit-allel.
        '''

        # add subparser
        pcangsd_parser = self.subparsers.add_parser(
            'pcangsd', help='Perform windowed PCA on genotype likelihoods (GL'
            ' or PL), accepts VCF, BEAGLE or TSV as input.'
        )

        # add shared arguments
        self.shared_arguments(pcangsd_parser)

        # define subparser-specific arguments
        variant_file = pcangsd_parser.add_mutually_exclusive_group(required=True)
        variant_file.add_argument(
            '-v', '--vcf', dest='variant_file_path', metavar='\b', help='Path to'
            ' (optionally gzipped) VCF file with GL or PL field.')
        variant_file.add_argument(
            '-b', '--beagle', dest='variant_file_path', metavar='\b', help='Path'
            ' to (optionally gzipped) BEAGLE file with GL or PL values.')
        variant_file.add_argument(
            '-t', '--tsv', dest='variant_file_path', metavar='\b', help='Path to'
            ' (optionally gzipped) TSV file with GL or PL values (see README'
            ' for specifications).')
        pcangsd_parser.add_argument(
            '-f', '--format', dest='gl_format', required=True, 
            choices=['GL', 'PL'], metavar='\b', help='Genotype likelihood'
            ' format ("GL" or "PL").')
        pcangsd_parser.add_argument(
            '-r', '--region', dest='region', required=True, metavar='\b',
            help='Genomic region in format "chrom:start-end".')
        pcangsd_parser.add_argument(
            '-s', '--samples', dest='samples', required=False, metavar='\b',
            help='Comma-separated list of samples to include or file with one'
            ' sample per line.')
        pcangsd_parser.add_argument(
            '-w', '--window_size', dest='w_size', required=False, type=int, 
            default=config.w_step, metavar='\b', help='Window size in base pairs'
            f' (bp) [default: {config.w_size}].')
        pcangsd_parser.add_argument(
            '-i', '--increment', dest='w_step', required=False, type=int,
            default=config.w_step, metavar='\b', help='Step size in base pairs'
            f' (bp). [default: {config.w_step}].')
        pcangsd_parser.add_argument(
            '-m', '--min_maf', dest='min_maf', required=False, type=float, 
            default=config.min_maf, metavar='\b', help='Minor allele frequency'
            f' threshold [default: {config.min_maf}].')
        pcangsd_parser.add_argument(
            '-p', '--polarize', dest='polarize', required=False, 
            default=config.pol_mode, choices=['auto', 'guide_samples', 'skip'], 
            metavar='\b', help='Sign polarization strategy'
            ' ("auto"/"guide_samples" or "skip") [default:'
            f' "{config.pol_mode}"].')
        pcangsd_parser.add_argument(
            '-g', '--guide_samples', dest='guide_samples', required=False, 
            metavar='\b', help='Applies only if "guide_samples" is selected for'
            ' -p/--polarize: One or more (-> comma-separated list) samples to'
            ' guide PC sign polarization.')
        pcangsd_parser.add_argument(
            '-@', '--threads', dest='threads', required=False, type=int, 
            default=config.n_threads, metavar='\b', help='Number of threads'
            f' [default: {config.n_threads}].')


    def polarize(self):
        '''
        (Re)-polarize windowed PC data from a previous run.
        '''

        # add subparser
        polarize_parser = self.subparsers.add_parser(
            'polarize', help='(Re)-polarize windowed PC data from a previous'
            ' run. Overwrites input data.'
        )

        # add shared arguments
        self.shared_arguments(polarize_parser)

        # define subparser-specific arguments
        polarize_parser.add_argument(
            '-c', '--principal_component', dest='pol_pc', required=False, 
            choices=['1', '2', 'both'], default=1, metavar='\b', help='Specify'
            f' which PC to re-polarize ("1", "2" or "both") [default: 1].')
        polarize_parser.add_argument(
            '-p', '--polarize', dest='polarize', required=False,
            default=config.pol_mode, choices=['auto', 'guide_samples'], 
            metavar='\b', help='Sign polarization strategy ("auto" or'
            f' "guide_samples") [default: {config.pol_mode}].')
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

        # add shared arguments
        self.shared_arguments(flip_parser)

        # define subparser-specific arguments
        flip_parser.add_argument(
            '-c', '--principal_component', dest='flip_pc', required=True, 
            choices=['1', '2', 'both'], default=config.flip_pc, metavar='\b', 
            help='Specify which PC to flip ("1", "2" or "both")'
            f' [default: {config.flip_pc}].')


    def chromplot(self):
        '''
        Plot PC1, PC2, heterozygosity and per window stats for all input 
        chromosomes.
        '''

        # add subparser
        chromplot_parser = self.subparsers.add_parser(
            'chromplot', help='Plot PC1, PC2, heterozygosity and per window'
            ' stats for all input chromosomes.'
        )

        # add shared arguments
        self.shared_arguments(chromplot_parser)

        # define subparser-specific arguments
        chromplot_parser.add_argument(
            '-r', '--region', dest='region', required=True, metavar='\b',
            help='Genomic region in format "chrom:start-end".')
        chromplot_parser.add_argument(
            '-m', '--metadata', dest='metadata_path', required=False, 
            metavar='\b', help='Path to metadata TSV where first column are '
            ' sample names. Other columns will be used to annotate data in HTML'
            ' plot.')
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
            default=config.plot_fmt, help='Output plot file format ("HTML",'
            ' "PDF", "SVG" or "PNG") [default: {config.plot_fmt}].')


    def genomeplot(self):
        '''
        Plot PC1, PC2, heterozygosity and per window stats for all input 
        chromosomes.
        '''

        # add subparser
        genomeplot_parser = self.subparsers.add_parser(
            'genomeplot', help='Plot PC1, PC2, heterozygosity and per window'
            ' stats for all input chromosomes.'
        )

        # define subparser-specific arguments
        genomeplot_parser.add_argument(
            '-p', '--prefix', dest='run_prefix', required=True, 
            metavar='\b', help='Prefix shared by all chromosomes runs to'
            ' include in genome-wide plot.')
        genomeplot_parser.add_argument(
            '-r', '--run_ids', dest='run_ids', required=True, 
            metavar='\b', help='Comma-separated list of run IDs to include,'
            ' format: e.g. {prefix}.{run_id}.pc_1.tsv.gz. Also used to'
            ' determine plotting order.')
        genomeplot_parser.add_argument(
            '-m', '--metadata', dest='metadata_path', required=False, 
            metavar='\b', help='Path to metadata TSV where first column are'
            ' sample names. Other columns will be used to annotate data in HTML'
            ' plot.')
        genomeplot_parser.add_argument(
            '-g', '--groups', dest='color_by', required=False, metavar='\b',
            help='Metadata column for color-grouping. Requires -m/--metadata.')
        genomeplot_parser.add_argument(
            '-c', '--colors', dest='hex_codes', required=False, metavar='\b',
            help='HEX codes (drop "#") for color groups. Check documentation'
            ' for formatting instructions. Requires -g/--groups.')
        genomeplot_parser.add_argument(
            '-i', '--interval', dest='interval', required=False, type=int, 
            default=config.plot_interval, metavar='\b', help='If set, only plot'
            ' values for every nth window (10 --> 10th)'
            f' [default: {config.plot_interval}.')
        genomeplot_parser.add_argument(
            '-f', '--format', dest='plot_fmt', required=False, metavar='\b',
            default=config.plot_fmt, help='Output plot file format ("HTML",'
            ' "PDF", "SVG" or "PNG") [default: {config.plot_fmt}].')


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
        if hasattr(args, 'samples'):
            if ',' in args.samples:
                sample_lst = args.samples.split(',')
            else:
                sample_lst = []
                with open(args.samples, 'r') as sample_file:
                    for line in sample_file:
                       sample_lst.append(line.strip().split('\t')[0])
        if hasattr(args, 'guide_samples'):
            if args.guide_samples:
                guide_sample_lst = args.guide_samples.split(',')
            else:
                guide_sample_lst = None
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
        if hasattr(args, 'samples'):
            if not ',' in args.samples:
                if not os.path.exists(args.samples):
                    self.parser.error(
                        args.samples + ': file does not exist.')
        if hasattr(args, 'run_ids'):
            for r_id in run_id_lst:
                 if not os.path.exists(f'{args.run_prefix}{r_id}.pc_1.tsv.gz'):
                    self.parser.error(
                        f'{args.run_prefix}{r_id}: run does not exist or is'
                        ' incomplete.')

        print(args)                                                             ### DELETE

        # convert args to dict and add to class as instance variable
        self.args_dct = vars(args)

        # add derived/processed arguments
        if hasattr(args, 'region'):
            self.args_dct['chrom'] = chrom
            self.args_dct['start'] = start
            self.args_dct['end'] = end
        if hasattr(args, 'samples'):
            self.args_dct['sample_lst'] = sample_lst
        if hasattr(args, 'guide_samples'):
            self.args_dct['guide_sample_lst'] = guide_sample_lst
        if hasattr(args, 'hex_codes'):
            self.args_dct['hex_code_dct'] = hex_code_dct
        if hasattr(args, 'plot_fmt'):
            self.args_dct['plot_fmt_lst'] = plot_fmt_lst
        if hasattr(args, 'run_ids'):
            self.args_dct['run_id_lst'] = run_id_lst
        
        # fill in default values from config if unset
        if not 'skip_monomorphic' in self.args_dct:
            self.args_dct['skip_monomorphic'] = config.skip_monomorphic
        if not 'min_var_per_w' in self.args_dct:  
            self.args_dct['min_var_per_w'] = config.min_var_per_w
        if not 'n_prev_windows' in self.args_dct:  
            self.args_dct['n_prev_windows'] = config.n_prev_windows
        if not 'pol_pc' in self.args_dct:  
            self.args_dct['pol_pc'] = config.pol_pc
        if not 'flip_pc' in self.args_dct:  
            self.args_dct['flip_pc'] = config.flip_pc
        if not 'chrom_plot_w' in self.args_dct:  
            self.args_dct['chrom_plot_w'] = config.chrom_plot_w
        if not 'chrom_plot_h' in self.args_dct:  
            self.args_dct['chrom_plot_h'] = config.chrom_plot_h
        
        #self.args_dct['min_maf'] = config.min_maf
        #self.args_dct['w_size'] = config.w_size
        #self.args_dct['w_step'] = config.w_step
        #self.args_dct['n_threads'] = config.n_threads
        #self.args_dct['pol_mode'] = config.pol_mode
        #self.args_dct['plot_fmt'] = config.plot_fmt
        #self.args_dct['plot_chroms'] = config.plot_chroms
        #self.args_dct['plot_interval'] = config.plot_interval

    
# region
# samples
# w_size
# w_step
# min_maf
# polarize
# guide_samples
# gl_format
# threads
# pc
# polarize
# sequence
# sequences
# metadata_path
# color_by
# hex_codes
# interval
# plot_fmt



