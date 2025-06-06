'''
Plot PC and associated data chromosome- and genome-wide.
'''


## IMPORT PACKAGES
import pandas as pd
import numpy as np
from plotly.subplots import make_subplots
import plotly.graph_objects as go
import plotly.colors as p_colors

## IMPORT CONFIG
from . import config

# MODULES
from modules.log import Log

## INSTANTIATE LOGGER
log = Log()


## CLASSES

class Plot:
    '''
    Plot windowed PC and associated data
    '''

    def __init__(self,                                                         # pylint: disable=W0102
                 plot_var,
                 stat_var=None,
                 prefix=None,
                 data=None,
                 pcs=None,
                 chrom=None,
                 start=None,
                 end=None,
                 run_prefix=None,
                 run_id_lst=None,
                 metadata_path=None,
                 color_by=None,
                 hex_code_dct=None,
                 numeric=None,
                 reverse=None,
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
        self.pcs = pcs
        self.chrom = chrom
        self.start = start
        self.end = end
        self.run_prefix = run_prefix
        self.run_id_lst = run_id_lst
        self.metadata_path = metadata_path
        self.color_by = color_by
        self.hex_code_dct = hex_code_dct
        self.numeric = numeric
        self.reverse = reverse
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

        # reformat plot_var and set display name
        if self.plot_var == 'hetp':
            self.plot_var_disp = 'SNP Heterozygosity'
        elif self.plot_var == 'miss':
            self.plot_var_disp = 'Missingness'
        else:
            self.plot_var_disp = f'PC {self.plot_var}'
            self.plot_var =  f'pc_{self.plot_var}'
        stat_var = f'{stat_var}_df'

        # define custom color scale
        # colors = ['blue', 'blue', 'purple', 'red', 'yellow']
        colors = ['#000080', '#0000FF', '#800080', '#FF0000', '#FFFF00']
        thresholds = [0, 0.25, 0.5, 0.75, 1]

        # create the colorscale
        colorscale = [
            [threshold, color] for threshold, color in zip(thresholds, colors)
        ]

        # set color scheme
        self.color_scale = colorscale # or set inbuilt schemes like 'Plasma'

        # set allowed NA strings
        self.na_lst = [None, 'NA', 'na', 'NaN']

    @staticmethod
    def subset(df, interval):
        '''
        Subset a dataframe to the specified interval (-i)
        '''

        return df.iloc[::interval, :]


    @staticmethod
    def is_hex(code):
        '''
        Return True if strings are in ok
        '''

        ok = 'abcdefABCDEF0123456789'
        return all(x in ok for x in code)


    def annotate(self):
        '''
        Annotate per-sample data with a metadata file if supplied, otherwise
        just reformat for plotting function
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
                log.error_nl('-m/--metadata: metadata file contains'
                          ' non-unique sample IDs')
            # subset and reorder metadata_df to match data_df individuals
            self.metadata_df = self.metadata_df.loc[
                self.metadata_df.index.intersection(sample_lst)
            ].reindex(sample_lst)

            # if individuals are missing in the metadata print error message
            if len(self.metadata_df) != len(sample_lst):
                log.newline()
                log.error_nl('-m/--metadata: one or more sample IDs are missing'
                          ' in the provided metadata file')

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
            value_name=self.plot_var_disp,
        )


    def set_colors(self):
        '''
        Parse per-sample plot color specifications and compile color_dct
        '''

        # fetch group_id (=color_by) if specified, else default to 'id'
        self.group_id = self.color_by if not self.color_by is None else 'id'

        # check if group exist in metadata column names
        if self.group_id != 'id' \
            and self.group_id not in self.metadata_df.columns:
            log.error_nl(
                f'-g/--groups: "{self.group_id}" is not the name of a column'
                 ' in the specified metadata file (-m/--metadata)'
            )

        # get list of groups/ids
        self.group_lst = sorted(set(self.data_df[self.group_id]))

        # set color_dct for continuous colors
        if self.numeric:

            # replace None (pandas) with np.nan (numpy), else through an error
            try:
                val_lst = np.array(
                    [
                        np.nan if x in self.na_lst else x
                        for x in self.group_lst
                    ],
                    dtype=np.float64
                )
            except:                                                            # pylint: disable=W0702
                log.error_nl('--n/--numeric: provided column (-g/--groups)'
                          ' contains non-numerical values')
            # scale to 0-1
            norm_val_lst = (
                val_lst-np.nanmin(val_lst))/\
                (np.nanmax(val_lst) - np.nanmin(val_lst)
            )
            # compile color_dct with all values as keys, using again None
            # instead of np.nan
            self.color_dct = {}
            for val, norm_val in zip(self.group_lst, norm_val_lst):
                if val in self.na_lst:
                    self.color_dct[val] = 'lightgrey'
                else:
                    self.color_dct[val] = p_colors.sample_colorscale(
                        self.color_scale, [norm_val],
                    )[0]
            # adjust plotting order
            self.group_lst = \
                [
                    x for x in self.group_lst \
                        if self.na_lst
                ] \
                + sorted(
                [
                        x for x in self.group_lst \
                            if x not in self.na_lst
                ])

        # define colors based on plotly default colors or specified HEX codes;
        # print error messages if HEX codes are missing for specified groups or
        # malformatted
        elif self.hex_code_dct:
            self.color_dct = self.hex_code_dct
            if not all(x in self.color_dct.keys() for x in self.group_lst):
                log.error_nl(
                    '-c/--colors: HEX codes missing for one or more groups'
                )
            elif not all(self.is_hex(x[1:]) for x in self.color_dct.values()) \
                or not all(len(x) == 7 for x in self.color_dct.values()):
                log.error_nl(
                    '-c/--colors: HEX codes not formatted correctly, please'
                    ' refer to the documentation')
            else:
                # set color_dct keys as group_lst to set plotting order
                self.group_lst = list(self.color_dct.keys())

        # use default colors
        else:

            def_col_lst = p_colors.DEFAULT_PLOTLY_COLORS
            self.color_dct = {
                self.group_lst[idx]: def_col_lst[idx % len(def_col_lst)] \
                    for idx in range(len(self.group_lst))
            }

        # reverse plotting order if specified
        if self.reverse:
            self.group_lst.reverse()


    def savefig(self):
        '''
        Save figure in HTML and/or other (PDF, SVG, PNG) format(s)
        '''
        for fmt in self.plot_fmt_lst:
            if fmt == 'html':
                self.fig.write_html(f'{self.prefix}.{self.plot_var}.{fmt}')
            else:
                self.fig.write_image(f'{self.prefix}.{self.plot_var}.{fmt}')


    def chromplot(self):
        '''
        Plot per-sample values for one chromosome (e.g. PC 1) with small panel
        of per window values (e.g. PC 1 variance explained) on top
        '''

        # LOAD & PREPARE DATA

        # load data
        self.data_df = getattr(self.data, f'{self.plot_var}_df')
        self.stat_df = getattr(self.data, 'stat_df')

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
        if self.stat_var == 'n_var':
            display_name = 'Variants per window'
            unit = ''
        elif self.stat_var == 'w_size':
            display_name = 'Window size'
            unit = ' bp'
        else:
            display_name = '% variance explained'
            unit = '%'

        # create mask of stretches of None
        split_mask = self.stat_df[[self.stat_var]].notna().all(axis=1)

        # series defining groups
        group_srs = (split_mask != split_mask.shift()).cumsum()

        # derive sub_dfs that contain no NaN stretches
        stat_sub_dfs = [
            sub_df for _, sub_df in self.stat_df[split_mask].groupby(group_srs)
        ]

        # only show first trace as legend item
        add_legend_item = True

        # plot each sub_df
        for sub_df in stat_sub_dfs:

            # compile per-window hover data strings

            hover_data = [
                ''.join(
                    [f'<b>pos</b>: {str(idx)}<br><b>{display_name}</b>: \
                    <b>{row[self.stat_var]}{unit}<br>' ]
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
                    hoverlabel={'font_size': 8},
                    line={'color':'#4d61b0', 'width': 1},
                    fill='tozeroy',
                    connectgaps=False,
                ),
            row=1, col=1)

            # set to False to stop adding additional legend items after first
            add_legend_item = False


        # BOTTOM PANEL

        # set show_legend to false if using a color scale
        show_legend = not self.numeric

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
                        [
                            f'<b>{c}</b>: {row[c]}<br>'
                            for c in sample_df.columns
                        ]
                    ) for i, row in sample_df.iterrows()
                ]

                # append x, y and hover values (separated by None to isolate
                # individuals plotted as part of the same trace)
                x_val_lst += sample_df['pos'].tolist() + [None]
                y_val_lst += sample_df[self.plot_var_disp].tolist() + [None]
                hover_str_lst += hover_data + [None]

            # plot
            self.fig.add_trace(
                go.Scatter(
                    x=x_val_lst,
                    y=y_val_lst,
                    text=hover_str_lst,
                    hoverinfo='text',
                    hoverlabel={'font_size': 8},
                    name=list(sample_df[self.group_id])[0],
                    legendgroup=list(sample_df[self.group_id])[0],
                    mode='lines',
                    line={'color': self.color_dct[group]},
                    connectgaps=False,
                    showlegend=show_legend,
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
            legend={'font': {'size': 10}},
            )

        # format lines
        self.fig.update_traces(
            row=1, col=1,
            line={'width': .7, 'color': 'lightgrey'},
        )
        self.fig.update_traces(
            row=2, col=1,
            line={'width': .7},
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
            ticks='outside', tickfont={'size': 10}, tickformat=',.0f',
            title_font={'size': 12},
            title={'text': '<b>Genomic position (bp)', 'standoff': 10}
        )

        # format y axis
        self.fig.update_yaxes(
            row=1, col=1,
            linewidth=1,
            side='left', mirror=True,
            ticks='outside', tickfont={'size': 10},
        )
        self.fig.update_yaxes(
            row=2, col=1,
            linewidth=1,
            side='left', mirror=True,
            ticks='outside', tickfont={'size': 10},
            title_font={'size': 12},
            title={'text': f'<b>{self.plot_var_disp}', 'standoff': 0},
        )

        # plot colorscale instead of per-sample legend for numeric metadata
        if self.numeric:
            val_lst = [
                float(x) for x in self.data_df[self.color_by] \
                    if not x in self.na_lst
            ]
            min_val = min(val_lst)
            max_val = max(val_lst)
            self.fig.update_layout(
                coloraxis={
                    'colorscale': self.color_scale,
                    'cmin': min_val,
                    'cmax': max_val,
                    'colorbar': {
                        'len': 0.7,
                        'thickness': 10,
                        'title': {
                            'text': self.color_by,
                            'font': {'size': 10},
                            'side': 'right'
                        },
                        'tickvals': [min_val, max_val],
                        'tickfont': {'size': 10},
                    },
                },
            )
            self.fig.add_trace(go.Scatter(
                x=[None],
                y=[None],
                mode='markers',
                marker={
                    'color': [0],
                    'coloraxis': 'coloraxis'
                },
                showlegend=False,
            ),
            row=2, col=1)


    def genomeplot(self):
        '''
        Plot per-sample values for multiple chromosomes (e.g. PC 1)
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
            self.data = WPCAData(
                f'{self.run_prefix}{run_id}',
                self.pcs,
            )
            self.data_df = getattr(self.data, f'{self.plot_var}_df')

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

        # set show_legend to false if using a color scale
        show_legend = not self.numeric

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
                        [
                            f'<b>{c}</b>: {row[c]}<br>'
                            for c in sample_df.columns
                        ]
                    ) for i, row in sample_df.iterrows()
                ]

                # append x, y and hover values (separated by None to isolate
                # individuals plotted as part of the same trace)
                x_val_lst += sample_df['genome_pos'].tolist() + [None]
                y_val_lst += sample_df[self.plot_var_disp].tolist() + [None]
                hover_str_lst += hover_data + [None]

            # plot
            self.fig.add_trace(
                go.Scatter(
                    x=x_val_lst,
                    y=y_val_lst,
                    text=hover_str_lst,
                    hoverinfo='text',
                    hoverlabel={'font_size': 8},
                    name=list(sample_df[self.group_id])[0],
                    legendgroup=list(sample_df[self.group_id])[0],
                    mode='lines',
                    line={'color': self.color_dct[group]},
                    connectgaps=False,
                    showlegend=show_legend,
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
            legend={'font': {'size': 10}},
            )

        # format lines
        self.fig.update_traces(
            line={'width': .7,},
        )

        # format x axis
        self.fig.update_xaxes(
            range=[
                min(self.data_df['genome_pos']),
                max(self.data_df['genome_pos']),
            ],
            linewidth=1,
            side='bottom', mirror=True,
            ticks='outside', tickfont={'size': 10}, tickformat=',.0f',
            tickvals=chrom_mid_lst,
            ticktext=['<b>' + x for x in self.run_id_lst],
            )

        # format y axis
        self.fig.update_yaxes(
            linewidth=1,
            side='left', mirror=True,
            ticks='outside', tickfont={'size': 10},
            title_font={'size': 12},
            title={'text': f'<b>{self.plot_var_disp}', 'standoff': 0},
        )

        # add vertical lines to separate chromosomes
        for genome_pos in chrom_end_lst[:-1]:

            self.fig.add_shape(
                {
                    'type': 'line',
                    'x0': genome_pos,
                    'x1': genome_pos,
                    'y0': 0,
                    'y1': 1,
                    'yref': 'paper',
                    'line': {
                        'color': '#000000',
                        'width': .8,
                        'dash': '1px, 1px',
                    },
                }
            )

        # plot colorscale instead of per-sample legend for numeric metadata
        if self.numeric:
            val_lst = [
                float(x) for x in self.data_df[self.color_by] \
                    if not x in self.na_lst
            ]
            min_val = min(val_lst)
            max_val = max(val_lst)
            self.fig.update_layout(
                coloraxis={
                    'colorscale': self.color_scale,
                    'cmin': min_val,
                    'cmax': max_val,
                    'colorbar': {
                        'len': 0.7,
                        'thickness': 10,
                        'title': {
                            'text': self.color_by,
                            'font': {'size': 10},
                            'side': 'right'
                        },
                        'tickvals': [min_val, max_val],
                        'tickfont': {'size': 10},
                    },
                },
            )
            self.fig.add_trace(go.Scatter(
                x=[None],
                y=[None],
                mode='markers',
                marker={
                    'color': [0],
                    'coloraxis': 'coloraxis'
                },
                showlegend=False,
            ),
            )

        # save image (using self.prefix to determine output prefix)
        if self.run_prefix.endswith('/'):
            self.prefix = f'{self.run_prefix}genomeplot'
        else:
            self.prefix = f'{self.run_prefix}.genomeplot'
