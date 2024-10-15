'''
Plot PC and associated data chromosome- and genome-wide.
'''

## IMPORT CONFIG
from . import config

## IMPORT PACKAGES
import pandas as pd
import numpy as np
from plotly.subplots import make_subplots
import plotly.graph_objects as go

# MODULES
from modules.log import Log

## INSTANTIATE LOGGER
log = Log()


## CLASSES

class Plot:
    '''
    Plot windowed PC and associated data.
    '''

    def __init__(self,
                 plot_var,
                 stat_var=None,
                 prefix=None,
                 data=None,
                 chrom=None,
                 start=None,
                 end=None,
                 run_prefix=None,
                 run_id_lst=None,
                 metadata_path=None,
                 color_by=None,
                 hex_code_dct=None,
                 interval=config.PLOT_INTERVAL,
                 chromplot_w=config.CHROMPLOT_W,
                 chromplot_h=config.CHROMPLOT_H,
                 genomeplot_w=config.GENOMEPLOT_W,
                 genomeplot_h=config.GENOMEPLOT_H,
                 plot_fmt_lst=[config.PLOT_FMT],
                 ):

        # input variables
        self.plot_var = plot_var
        self.stat_var = stat_var
        self.prefix = prefix
        self.data = data
        self.chrom = chrom
        self.start = start
        self.end = end
        self.run_prefix = run_prefix
        self.run_id_lst = run_id_lst
        self.metadata_path = metadata_path
        self.color_by = color_by
        self.hex_code_dct = hex_code_dct
        self.interval = interval
        self.chromplot_w = chromplot_w
        self.chromplot_h = chromplot_h
        self.genomeplot_w = genomeplot_w
        self.genomeplot_h = genomeplot_h
        self.plot_fmt_lst = plot_fmt_lst

        # instance variables
        self.data_df = pd.DataFrame()
        self.stat_df = pd.DataFrame()
        self.metadata_df = pd.DataFrame()
        self.group_id = None
        self.group_lst = None
        self.color_dct = None
        self.fig = None

        # set plot_var display name
        if   self.plot_var == 'pc_1':
            self.plot_var_disp = 'PC 1'
        elif self.plot_var == 'pc_2':
            self.plot_var_disp = 'PC 2'
        elif self.plot_var == 'hetp':
            self.plot_var_disp = 'SNP Heterozygosity'

    @staticmethod
    def subset(df, interval):
        '''
        Subset a dataframe to the specified interval (-i).
        '''

        return df.iloc[::interval, :]


    def annotate(self):
        '''
        Annotate per-sample data with a metadata file if supplied, otherwise
        just reformat for plotting function.
        '''

        # fetch sample names and order (=data_df column names)
        sample_lst = list(self.data_df.columns)

        # transpose data_df
        self.data_df = self.data_df.T

        # initiate id_vars_lst for hover display (id, chrom are default)
        id_var_lst = ['id', 'chrom']

        # read metadata if provided, do sanity checks and use for annotation
        if self.metadata_path:

            # read metadata and print error message if there are non-unique IDs
            self.metadata_df = pd.read_csv(
                self.metadata_path, sep='\t', index_col=0, dtype=str
            )
            if len(self.metadata_df.index) != len(set(self.metadata_df.index)):
                log.newline()
                log.error('The provided metadata file contains non-unique'
                          ' sample IDs')
                log.newline()


            # subset and reorder metadata_df to match data_df individuals
            self.metadata_df = self.metadata_df.loc[
                self.metadata_df.index.intersection(sample_lst)
            ].reindex(sample_lst)

            # if individuals are missing in the metadata print error message
            if len(self.metadata_df) != len(sample_lst):
                log.newline()
                log.error('One or more sample IDs are missing in the'
                          ' provided metadata file')
                log.newline()


            # add metadata columns to data_df
            for column_name in self.metadata_df.columns:
                self.data_df[column_name] = list(self.metadata_df[column_name])

            # add metadata column names to id_var_lst
            id_var_lst += list(self.metadata_df.columns)

        # add id (from index) and chrom columns
        self.data_df['id'] = list(self.data_df.index)
        self.data_df['chrom'] = self.chrom

        # replace numpy NaN with 'NA' for plotting (hover_data display)
        self.data_df = self.data_df.replace(np.nan, 'NA')

        # convert to long format for plotting
        self.data_df = pd.melt(
            self.data_df,
            id_vars=id_var_lst,
            var_name='pos',
            value_name=self.plot_var,
        )


    def set_colors(self):
        '''
        Parse per-sample plot color specifications and compile color_dct.
        '''

        # fetch group_id (=color_by) if specified, else default to 'id'
        self.group_id = self.color_by if not self.color_by is None else 'id'

        # get list of groups/ids
        self.group_lst = list(set(self.data_df[self.group_id]))

        # define colors based on plotly default colors or specified HEX codes;
        # print error messages if HEX codes are missing for specified groups
        if self.hex_code_dct:
            self.color_dct = self.hex_code_dct
            if not all(x in self.color_dct.keys() for x in self.group_lst):
                log.newline()
                log.error('HEX codes missing for one or more groups')
                log.newline()

        else:
            import plotly.colors as pc
            def_col_lst = pc.DEFAULT_PLOTLY_COLORS
            self.color_dct = {
                self.group_lst[idx]: def_col_lst[idx % len(def_col_lst)] \
                    for idx in range(len(self.group_lst))
            }

    def savefig(self):
        '''
        Save figure in HTML and/or other (PDF, SVG, PNG) format(s).
        '''
        for fmt in self.plot_fmt_lst:
            if fmt == 'html':
                self.fig.write_html(f'{self.prefix}.{fmt}')
            else:
                self.fig.write_image(f'{self.prefix}.{fmt}')


    def chromplot(self):
        '''
        Plot per-sample values for one chromosome (e.g. PC 1) with small panel
        of per window values (e.g. PC 1 variance explained) on top.
        '''

        # LOAD & PREPARE DATA

        # load data
        self.data_df = getattr(self.data, self.plot_var)
        self.stat_df = getattr(self.data, 'stat')

        # subset if interval (-i) is specified
        if self.interval:
            self.data_df = self.subset(self.data_df, self.interval)
            self.stat_df = self.subset(self.stat_df, self.interval)

        # annotate per-sample data if metadata were supplied
        self.annotate()

        # set per-sample plot colors
        self.set_colors()

        # figure setup
        self.fig = make_subplots(
            rows=2, cols=1,
            row_heights=[1, 6],
            vertical_spacing=0.0,
            shared_xaxes=True,
            )


        # TOP PANEL

        # parse display name for top panel: variance explained or n of sites
        display_name = \
            '% heterozygous sites' if self.stat_var == 'hetp' else \
            '% variance explained'

        # create mask of stretches of None
        split_mask = self.stat_df[[self.stat_var]].notna().all(axis=1)

        # series defining groups
        group_srs = (split_mask != split_mask.shift()).cumsum()

        # derive sub_dfs that contain no NaN stretches
        stat_sub_dfs = [
            sub_df for _, sub_df in self.stat_df[split_mask].groupby(group_srs)]

        # only show first trace as legend item
        add_legend_item = True

        # plot each sub_df
        for sub_df in stat_sub_dfs:

            # compile per-window hover data strings
            hover_data = [
                ''.join(
                    [f'<b>pos</b>: {str(idx)}<br><b>{display_name}</b>: \
                    <b>{row[self.stat_var]}%<br>' ]
                ) for idx, row in sub_df.iterrows()
            ]

            # plot
            self.fig.add_trace(
                go.Scatter(
                    x=sub_df.index,
                    y=sub_df[self.stat_var],
                    name=display_name,
                    legendgroup=display_name,
                    showlegend=add_legend_item,
                    mode='lines',
                    text=hover_data,
                    hoverinfo='text',
                    line=dict(color='#4d61b0', width=1),
                    fill='tozeroy',
                    connectgaps=False,  # Disable connecting gaps
                ),
            row=1, col=1)

            # set to False to stop adding additional legend items after first
            add_legend_item = False


        # BOTTOM PANEL

        # plot each specified group (-g) separately (or 'id' if unspecified)
        for group in self.group_lst:

            # subset data to group
            group_df = self.data_df[self.data_df[self.group_id] == group]

            # initiate lists to hold per-sample-per-window x values, y values
            # and hover data strings
            x_val_lst = []
            y_val_lst = []
            hover_str_lst = []

            # iterate through each individual per group
            for sample in set(group_df['id']):

                # subset data
                sample_df = group_df[group_df['id'] == sample]

                # compile hover text string for each window
                hover_data = [
                    ''.join(
                        [f'<b>{c}</b>: {row[c]}<br>' for c in sample_df.columns]
                    ) for i, row in sample_df.iterrows()
                ]

                # append x, y and hover values (separated by None to isolate
                # individuals plotted as part of the same trace)
                x_val_lst += sample_df['pos'].tolist() + [None]
                y_val_lst += sample_df[self.plot_var].tolist() + [None]
                hover_str_lst += hover_data + [None]

            # plot
            self.fig.add_trace(
                go.Scatter(
                    x=x_val_lst,
                    y=y_val_lst,
                    text=hover_str_lst,
                    hoverinfo='text',
                    name=list(sample_df[self.group_id])[0],
                    legendgroup=list(sample_df[self.group_id])[0],
                    mode='lines',
                    line=dict(color=self.color_dct[group]),
                    connectgaps=False,
                ),
                row=2, col=1
            )

        # general layout
        self.fig.update_layout(
            template='simple_white',
            paper_bgcolor='rgba(0,0,0,0)',
            plot_bgcolor='rgba(0,0,0,0)',
            font_family='Arial',
            font_color='black',
            width=self.chromplot_w,
            height=self.chromplot_h,
            legend=dict(font=dict(size=10)),
            )

        # format lines
        self.fig.update_traces(
            row=1, col=1,
            line=dict(width=.7, color='lightgrey'),
        )
        self.fig.update_traces(
            row=2, col=1,
            line=dict(width=.7),
        )

        # format x axis
        self.fig.update_xaxes(
            row=1, col=1,
            range=[self.start, self.end],
            linewidth=1,
            side='top', mirror=False,
            ticks='', showticklabels=False,
        )
        self.fig.update_xaxes(
            row=2, col=1,
            range=[self.start, self.end],
            linewidth=1,
            side='bottom', mirror=True,
            ticks='outside', tickfont=dict(size=10), tickformat=',.0f',
            title_font=dict(size=12),
            title=dict(text='<b>Genomic position (bp)', standoff=10))

        # format y axis
        self.fig.update_yaxes(
            row=1, col=1,
            linewidth=1,
            side='left', mirror=True,
            ticks='outside', tickfont=dict(size=10),
        )
        self.fig.update_yaxes(
            row=2, col=1,
            linewidth=1,
            side='left', mirror=True,
            ticks='outside', tickfont=dict(size=10),
            title_font=dict(size=12),
            title=dict(text=self.plot_var_disp, standoff=0),
        )

        # save image
        self.savefig()


    def genomeplot(self):
        '''
        Plot per-sample values for multiple chromosomes (e.g. PC 1).
        '''

        # import data module
        from modules.data import WPCAData


        # LOAD & PREPARE DATA

        # initite lists to hold data_dfs from different chromosomes to be
        # concatenated, end and mid of chromosome positions for plot formatting
        data_df_lst = []
        chrom_end_lst = []
        chrom_mid_lst = []

        # initiate chromosome offset
        offset = 0

        # sequentially read in data from input chromosomes (-r)
        for run_id in self.run_id_lst:

            # load data
            data = WPCAData(self.run_prefix + run_id)
            self.data_df = getattr(data, self.plot_var)

            # subset if interval (-i) is specified
            if self.interval:
                self.data_df = self.subset(self.data_df, self.interval)

            # annotate per-sample data if metadata were supplied
            self.chrom = run_id
            self.annotate()

            # append chrom mid/end positions to lists
            chrom_end_lst.append(offset + max(self.data_df['pos']))
            chrom_mid_lst.append(offset + self.data_df['pos'].median())

            # add genome-wide plotting position as data_df column
            self.data_df['genome_pos'] = self.data_df['pos'] + offset

            # update offset
            offset += max(self.data_df['pos'])

            # append data_df to list
            data_df_lst.append(self.data_df)

        # concatenate chromosomal data_dfs
        self.data_df = pd.concat(data_df_lst)

        # set per-sample plot colors
        self.set_colors()

        # figure setup
        self.fig = go.Figure()

        # PLOT

        # plot each specified group (-g) separately (or 'id' if unspecified)
        for group in self.group_lst:

            # subset data to group
            group_df = self.data_df[self.data_df[self.group_id] == group]

            # initiate lists to hold per-sample-per-window x values, y values
            # and hover data strings
            x_val_lst = []
            y_val_lst = []
            hover_str_lst = []

            # iterate through each individual per group
            for sample in set(group_df['id']):

                # subset data
                sample_df = group_df[group_df['id'] == sample]

                # compile hover text string for each window
                hover_data = [
                    ''.join(
                        [f'<b>{c}</b>: {row[c]}<br>' for c in sample_df.columns]
                    ) for i, row in sample_df.iterrows()
                ]

                # append x, y and hover values (separated by None to isolate
                # individuals plotted as part of the same trace)
                x_val_lst += sample_df['genome_pos'].tolist() + [None]
                y_val_lst += sample_df[self.plot_var].tolist() + [None]
                hover_str_lst += hover_data + [None]

            # plot
            self.fig.add_trace(
                go.Scatter(
                    x=x_val_lst,
                    y=y_val_lst,
                    text=hover_str_lst,
                    hoverinfo='text',
                    name=list(sample_df[self.group_id])[0],
                    legendgroup=list(sample_df[self.group_id])[0],
                    mode='lines',
                    line=dict(color=self.color_dct[group]),
                    connectgaps=False,
                ),
            )

        # general layout
        self.fig.update_layout(
            template='simple_white',
            paper_bgcolor='rgba(0,0,0,0)',
            plot_bgcolor='rgba(0,0,0,0)',
            font_family='Arial',
            font_color='black',
            width=self.genomeplot_w,
            height=self.genomeplot_h,
            legend=dict(font=dict(size=10)),
            )

        # format lines
        self.fig.update_traces(
            line=dict(width=.7,),
        )

        # format x axis
        self.fig.update_xaxes(
            range=[
                min(self.data_df['genome_pos']),
                max(self.data_df['genome_pos']),
            ],
            linewidth=1,
            side='bottom', mirror=True,
            ticks='outside', tickfont=dict(size=10), tickformat=',.0f',
            tickvals=chrom_mid_lst,
            ticktext=['<b>' + x for x in self.run_id_lst],
            )

        # format y axis
        self.fig.update_yaxes(
            linewidth=1,
            side='left', mirror=True,
            ticks='outside', tickfont=dict(size=10),
            title_font=dict(size=12),
            title=dict(text=self.plot_var_disp, standoff=0),
        )

        # add vertical lines to separate chromosomes
        for genome_pos in chrom_end_lst[:-1]:

            self.fig.add_shape(
                dict(
                    type='line',
                    x0=genome_pos, x1=genome_pos,
                    y0=0, y1=1, yref='paper',
                    line=dict(color='#000000', width=.8, dash='1px, 1px',),
                )
            )

        # save image (using self.prefix to determine output prefix)
        self.prefix = f'{self.run_prefix}{self.plot_var}_genomeplot'
        self.savefig()
