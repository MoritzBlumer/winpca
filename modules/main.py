#!/usr/bin/env python3

"""
WinPCA. A package for windowed principal component analysis.
"""


## IMPORT CONFIG
from . import config

## IMPORT PACKAGES
import os                                                                      # pylint: disable=C0411

## IMPORT MODULES
from .cli import CLI
from .log import Log

## INSTANTIATE LOGGER
log = Log()



## MAIN

def run():
    '''
    Execute WinPCA.
    '''

    ## PARSE COMMAND LINE ARGUMENTS

    # instantiate, call subparsers & parse
    cli = CLI()
    cli.pca()
    cli.polarize()
    cli.orient()
    cli.flip()
    cli.chromplot()
    cli.genomeplot()
    cli.parse_args()
    args = cli.args

    # control thread usage
    config.N_THREADS = args['n_threads']
    os.environ["OMP_NUM_THREADS"] = str(config.N_THREADS)
    os.environ["OPENBLAS_NUM_THREADS"] = str(config.N_THREADS)
    os.environ["MKL_NUM_THREADS"] = str(config.N_THREADS)
    os.environ["VECLIB_MAXIMUM_THREADS"] = str(config.N_THREADS)

    # now import core data module
    from .data import WPCAData

    # set mode
    mode = args['mode']


    ## MODE: WINDOWED PCA

    if mode == 'pca':

        # assess variant format
        config.VAR_FMT = args['var_fmt']

        # print info
        log.newline()
        log.info('Performing windowed PCA')
        log.newline()

        # determine file format
        if args['variant_file_path'].endswith('.gz'):
            file_fmt = args['variant_file_path'].split('.')[-2].upper()
        else:
            file_fmt = args['variant_file_path'].split('.')[-1].upper()

        # load module & instantiate
        from .windowed_pca import WPCA
        wpca = WPCA(
            variant_file_path=args['variant_file_path'],
            file_fmt=file_fmt,
            var_fmt=args['var_fmt'],
            pcs = args['pcs'],
            pcangsd_em_eig = args['pcangsd_em_eig'],
            pcangsd_em_iter = args['pcangsd_em_iter'],
            sample_lst=args['sample_lst'],
            query_sample_lst=args['query_sample_lst'],
            chrom=args['chrom'],
            start=args['start'],
            stop=args['end'],
            w_size=args['w_size'],
            w_step=args['w_step'],
            gt_min_var_per_w=args['gt_min_var_per_w'],
            gl_pl_min_var_per_w=args['gl_pl_min_var_per_w'],
            skip_monomorphic=args['skip_monomorphic'],
            gt_mean_impute=args['gt_mean_impute'],
            proc_trail_trunc_w=args['proc_trail_trunc_w'],
            vcf_pass_filter=args['vcf_pass_filter'],
            min_maf=args['min_maf'],
            x_mode=args['x_mode'],
            write_win_pos=args['write_win_pos'],
            n_threads=args['n_threads'],
        )

        # run, parse and write output
        wpca.window_parser()
        data = WPCAData(
            args['prefix'],
            args['pcs'],
            wpca,
        )


    ## ELSE READ IN EXISTING DATA:

    if mode in ['polarize', 'flip', 'chromplot']:

        # read data
        data = WPCAData(
            args['prefix'],
            args['pcs'],
        )


    ## MODE: POLARIZE

    if mode in ['pca', 'polarize'] \
        and not args['polarize'] == 'skip':

        # print info
        log.newline()
        log.info('Polarizing PC data')

        # load module & instantiate
        from .transform_data import Polarize
        polarize = Polarize()

        # adaptive (also used with 'local' and 'global')
        if args['polarize'] in ['adaptive', 'local', 'global']:
            for pc in args['pol_pcs']:
                data.modify_data(                                              # pylint: disable=E0606
                    f'pc_{pc}_df',
                    polarize.adaptive,
                    args['n_prev_windows'],
                )

        # local ML model (also used before 'global')
        if args['polarize'] in ['local', 'global']:
            for pc in args['pol_pcs']:
                data.modify_data(
                    f'pc_{pc}_df', polarize.local_ml,
                    args['n_draws_local'],
                    args['flip_threshold_local'],
                    args['weights_local_path'],
                    args['epsilon'],
                )

        # global ML model
        if args['polarize'] == 'global':
            for pc in args['pol_pcs']:
                data.modify_data(
                    f'pc_{pc}_df',
                    polarize.global_ml,
                    args['weights_global_path'],
                    args['epsilon'],
                )

        # using guide samples
        if args['polarize'] == 'guide_samples':

            # check if guide samples exist
            check_df = getattr(data, f'pc_{args["pol_pcs"][0]}_df')
            if not set(args['guide_sample_lst']).issubset(check_df.columns):
                missing_lst = \
                    set(args['guide_sample_lst']) - set(check_df.columns)
                data.to_files()
                log.error_nl(
                    f'-g/--guide_samples: {", ".join(missing_lst)} not found,'
                     ' writing unpolarized output instead'
                )

            for pc in args['pol_pcs']:
                data.modify_data(
                    f'pc_{pc}_df',
                    polarize.guide_samples,
                    args['guide_sample_lst'],
                )


    ## MODE: ORIENT

    if mode == 'orient':

        # print info
        log.newline()
        log.info('Orienting all chromosomes consistently')

        # load module & instantiate
        from .transform_data import Orient
        orient = Orient(
            run_prefix=args['run_prefix'],
            run_id_lst=args['run_id_lst'],
            pcs=args['pcs'],
            orient_pcs=args['pol_pcs'],
        )

        # read all chromosomes, flip where needed, write
        orient.orient(args['epsilon'])
        orient.to_files()


    ## MODE: FLIP

    if mode == 'flip':

        # print info
        log.newline()
        log.info('Flipping PC data')

        # load module & instantiate
        from .transform_data import Flip
        flip = Flip()

        # reflect (entire)
        if args['reflect']:
            for pc in args['flip_pcs']:
                data.modify_data(
                    f'pc_{pc}_df', flip.flip_chrom
                )

        # flip specified windows
        if args['flip_windows']:
            for pc in args['flip_pcs']:
                data.modify_data(
                    f'pc_{pc}_df', flip.flip_windows, args['flip_window_lst']
                )


    ## WRITE DATA TO FILES

    # except plot modes
    if mode in ['pca', 'polarize', 'flip']:
        data.to_files()


    ## MODE: CHROMPLOT

    if mode  == 'chromplot':

        # print info
        log.newline()
        log.info('Creating chromosome plot')

        # infer stat_var
        if args['plot_var'] == 'hetp':
            if len(set(data.stat_df['n_var'][:-1])) == 1:
                stat_var = 'w_size'
            else:
                stat_var = 'n_var'
        else:
            stat_var = f'pc_{args["plot_var"]}_ve'

        # load module & instantiate
        from .plot import Plot
        plot = Plot(args['plot_var'],
                    stat_var=stat_var,
                    prefix=args['prefix'],
                    data=data,
                    pcs=args['pcs'],
                    chrom=args['chrom'],
                    start=args['start'],
                    end=args['end'],
                    metadata_path=args['metadata_path'],
                    locations_path=args['locations_path'],
                    color_by=args['color_by'],
                    hex_code_dct=args['hex_code_dct'],
                    interval=args['interval'],
                    chromplot_w=args['chromplot_w'],
                    chromplot_h=args['chromplot_h'],
                    plot_fmt_lst=args['plot_fmt_lst'],
                    numeric=args['numeric'],
                    reverse=args['reverse'],
        )
        plot.chromplot()
        plot.savefig()


    # MODE: GENOMEPLOT

    if mode == 'genomeplot':

        # print info
        log.newline()
        log.info('Creating genome-wide plot')

        # load module & instantiate
        from .plot import Plot
        plot = Plot(args['plot_var'],
                    pcs=args['pcs'],
                    run_prefix=args['run_prefix'],
                    run_id_lst=args['run_id_lst'],
                    metadata_path=args['metadata_path'],
                    locations_path=args['locations_path'],
                    color_by=args['color_by'],
                    hex_code_dct=args['hex_code_dct'],
                    interval=args['interval'],
                    genomeplot_w=args['genomeplot_w'],
                    genomeplot_h=args['genomeplot_h'],
                    plot_fmt_lst=args['plot_fmt_lst'],
                    numeric=args['numeric'],
                    reverse=args['reverse'],
        )
        plot.genomeplot()
        plot.savefig()


    # END

    # print info
    log.newline()
    log.info('Done')
    log.newline()



# EXECUTE
def main():
    '''
    Execute WinPCA from command line
    '''

    try:
        run()
    except KeyboardInterrupt:
        log.newline()
        log.info(
            'Exiting'
        )
        log.newline()
