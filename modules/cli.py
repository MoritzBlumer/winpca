import argparse

class CLI:

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


    def shared_arguments(self, subparser):
        '''
        Arguments shared between all sub-commands.
        '''
        # define arguments that are shared across all subparsers
        subparser.add_argument(
            '-n', '--name', dest='prefix', required=True, help='Name of the'
            ' run, i.e. prefix for all results generated in this WinPCA'
            ' analysis.')  


    def pca(self):
        '''
        Windowed PCA on called genotypes (GT) with scikit-allel.
        '''

        # add subparser
        pca_parser = self.subparsers.add_parser(
            'pca', help='Perform windowed PCA on called genotypes (GT),'
            'accepts VCF or TSV as input.'
        )

        # add shared arguments
        self.shared_arguments(pca_parser)

        # define subparser-specific arguments
        variant_file = pca_parser.add_mutually_exclusive_group(required=True)
        variant_file.add_argument(
            '-v', '--vcf', dest='variant_file_path', help='Path to (optionally'
            ' gzipped) VCF file with GT field.')
        variant_file.add_argument(
            '-t', '--tsv', dest='variant_file_path', help='Path to (optionally'
            ' gzipped) TSV file with genotypes (see README for '
            ' specifications).')
        pca_parser.add_argument(
            '-r', '--region', dest='region', required=True, help='Genomic'
            ' region in format chrom:start-end.')
        pca_parser.add_argument(
            '-s', '--samples', dest='samples', required=False, help='Comma-'
            ' separated list of samples to include or file with one sample per'
            ' line.')
        pca_parser.add_argument(
            '-w', '--window_size', dest='w_size', required=False, type=int, 
            default=1000000, help='Window size in base pairs (bp).')
        pca_parser.add_argument(
            '-i', '--increment', dest='w_step', required=False, type=int,
            default=10000, help='Step size in base pairs (bp).')
        pca_parser.add_argument(
            '-m', '--min_maf', dest='min_maf', required=False, type=float, 
            default=0.01, help='Minor allele frequency threshold.')
        pca_parser.add_argument(
            '-p', '--polarize', dest='polarize', required=False, default='auto', 
            choices=['auto', 'guide_samples', 'skip'], help='Sign polarization'
            ' strategy (or skip).')
        pca_parser.add_argument(
            '-g', '--guide_samples', dest='guide_samples', required=False, \
            help='Applies only if "guide_samples" is selected for'
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
            '-v', '--vcf', dest='variant_file_path', help='Path to (optionally'
            ' gzipped) VCF file with GL or PL field.')
        variant_file.add_argument(
            '-b', '--beagle', dest='variant_file_path', help='Path to'
            ' (optionally gzipped) BEAGLE file with GL or PL values.')
        variant_file.add_argument(
            '-t', '--tsv', dest='variant_file_path', help='Path to (optionally'
            ' gzipped) TSV file with GL or PL values (see README for'
            ' specifications).')
        pcangsd_parser.add_argument(
            '-f', '--format', dest='gl_format', required=True, 
            choices=['GL', 'PL'], help='Genotype likelihood format.')
        pcangsd_parser.add_argument(
            '-r', '--region', dest='region', required=True, help='Genomic'
            ' region in format chrom:start-end.')
        pcangsd_parser.add_argument(
            '-s', '--samples', dest='samples', required=False, help='Comma-'
            ' separated list of samples to include or file with one sample per'
            ' line.')
        pcangsd_parser.add_argument(
            '-w', '--window_size', dest='w_size', required=False, type=int, 
            default=1000000, help='Window size in base pairs (bp).')
        pcangsd_parser.add_argument(
            '-i', '--increment', dest='w_step', required=False, type=int, 
            default=10000,  help='Step size in base pairs (bp).')
        pcangsd_parser.add_argument(
            '-m', '--min_maf', dest='min_maf', required=False, type=float, 
            default=0.01, help='Minor allele frequency threshold.')
        pcangsd_parser.add_argument(
            '-p', '--polarize', dest='polarize', required=False, default='auto', 
            choices=['auto', 'guide_samples', 'skip'], help='Sign polarization'
            ' strategy (or skip).')
        pcangsd_parser.add_argument(
            '-g', '--guide_samples', dest='guide_samples', required=False, \
            help='Applies only if "guide_samples" is selected for'
            ' -p/--polarize: One or more (-> comma-separated list) samples to'
            ' guide PC sign polarization.')
        pcangsd_parser.add_argument(
            '-@', '--threads', dest='threads', required=False, type=int, 
            default=2, help='Number of threads.')


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
            '-c', '--principal_component', dest='pc', required=False, \
            choices=['1', '2', 'both'], default='both', help='Specify which PC'
            ' to re-polarize or select both.')
        polarize_parser.add_argument(
            '-p', '--polarize', dest='polarize', required=False, \
            default='auto', choices=['auto', 'guide_samples'], help='Sign'
            ' polarization strategy.')
        polarize_parser.add_argument(
            '-g', '--guide_samples', dest='guide_samples', required=False, \
            help='Applies only if "guide_samples" is selected for'
            ' -p/--polarize: One or more (-> comma-separated list) samples to'
            ' guide PC sign polarization.')


    def flip(self):
        '''
        Flip/reflect windowed PC data from a previous run (multiply values by
        -1).
        '''

        # add subparser
        flip_parser = self.subparsers.add_parser(
            'flip', help='Flip/reflect windowed PC data from a previous'
            ' run (multiply values by -1). Overwrites input data.'
        )

        # add shared arguments
        self.shared_arguments(flip_parser)

        # define subparser-specific arguments
        flip_parser.add_argument(
            '-c', '--principal_component', dest='pc', required=True, \
            choices=['1', '2', 'both'], default='both', help='Specify which PC'
            ' to flip select both.')


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
            '-m', '--metadata', dest='metadata_path', required=False, 
            help='Path to metadata TSV where first column are sample names and'
            ' other columns will be used to annotate data in HTML plot.')
        chromplot_parser.add_argument(
            '-s', '--sequences', dest='sequences', required=False, 
            default='all', help='Comma-separated list of reference sequences '
            ' IDs, e.g. chromosomes to include (or to define plotting order).')
        chromplot_parser.add_argument(
            '-g', '--groups', dest='color_by', required=False,
            help='Metadata column for color-grouping. Requires -m/--metadata.')
        chromplot_parser.add_argument(
            '-c', '--colors', dest='hex_codes', required=False,
            help='HEX codes (drop "#") for color groups. Check documentation'
            ' for formatting instructions. Requires -g/--groups.')
        chromplot_parser.add_argument(
            '-i', '--interval', dest='interval', required=False, type=int,
            help='If set, only plot values for every nth window (10 --> 10th).')
        chromplot_parser.add_argument(
            '-f', '--format', dest='file_format', required=False, 
            default='both', choices=['HTML', 'PDF', 'both'], help='Output'
            ' plot file format.')


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

        # add shared arguments
        self.shared_arguments(genomeplot_parser)

        # define subparser-specific arguments
        genomeplot_parser.add_argument(
            '-m', '--metadata', dest='metadata_path', required=False, 
            help='Path to metadata TSV where first column are sample names and'
            ' other columns will be used to annotate data in HTML plot.')
        genomeplot_parser.add_argument(
            '-s', '--sequences', dest='sequences', required=False, 
            default='all', help='Comma-separated list of reference sequences '
            ' IDs, e.g. chromosomes to include (or to define plotting order).')
        genomeplot_parser.add_argument(
            '-g', '--groups', dest='color_by', required=False,
            help='Metadata column for color-grouping. Requires -m/--metadata.')
        genomeplot_parser.add_argument(
            '-c', '--colors', dest='hex_codes', required=False,
            help='HEX codes (drop "#") for color groups. Check documentation'
            ' for formatting instructions. Requires -g/--groups.')
        genomeplot_parser.add_argument(
            '-i', '--interval', dest='interval', required=False, type=int,
            help='If set, only plot values for every nth window (10 --> 10th).')
        genomeplot_parser.add_argument(
            '-f', '--format', dest='file_format', required=False, 
            default='both', choices=['HTML', 'PDF', 'both'], help='Output'
            ' plot file format.')     

        

        

    def parse_args(self):
        '''
        Parse command-line arguments.
        '''

        # parse arguments
        args = self.parser.parse_args()  
        
        # handle interdependent options
        print(args)
        if args.winpca == 'pca':
            if args.polarize == 'guide_samples' and not args.guide_samples:
                self.parser.error(
                    '-g/--guide_samples is required when chosing'
                    ' "guide_samples" set as polarization strategy (-p)')
        if args.winpca == 'pcangsd':
            if args.polarize == 'guide_samples' and not args.guide_samples:
                self.parser.error(
                    '-g/--guide_samples is required when chosing'
                    ' "guide_samples" set as polarization strategy (-p)')
        if args.winpca == 'polarize':
            if args.polarize == 'guide_samples' and not args.guide_samples:
                self.parser.error(
                    '-g/--guide_samples is required when chosing'
                    ' "guide_samples" set as polarization strategy (-p)')
        if args.winpca == 'chromplot':
            if args.color_by and not args.metadata_path:
                self.parser.error(
                    '-m/--metadata is required to infer -g/--groups.')        
            if args.hex_codes_str and not args.color_by:
                self.parser.error(
                    '-g/--groups is required if -c/--colors is set.')
        if args.winpca == 'genomeplot':
            if args.color_by and not args.metadata_path:
                self.parser.error(
                    '-m/--metadata is required to infer -g/--groups.')        
            if args.hex_codes_str and not args.color_by:
                self.parser.error(
                    '-g/--groups is required if -c/--colors is set.')

        return args
    

cli=CLI()

cli.pca()
cli.pcangsd()
cli.polarize()
cli.flip()
cli.genomeplot()


#cli.shared_arguments()

args = cli.parse_args()

print("Name:", args.prefix)
print("Subcommand:", args.winpca)

if args.winpca == 'gt_pca':
    print("Variant file path:", args.variant_file_path)
    print("Region:", args.region)
    print("Window size:", args.w_size)
    print("Step size:", args.w_step)
    print("Min MAF:", args.min_maf)
    print("Polarize:", args.polarize)

if args.winpca == 'pcangsd':
    print("Variant file path:", args.variant_file_path)
    print("Region:", args.region)
    print("Window size:", args.w_size)
    print("Step size:", args.w_step)
    print("Min MAF:", args.min_maf)
    print("Polarize:", args.polarize)


## add option that if guide_samples is selected, they can be specified