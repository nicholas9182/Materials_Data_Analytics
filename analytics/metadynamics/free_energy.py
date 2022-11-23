import pandas as pd
import plotly.graph_objects as go
import plotly.express as px
import plotly.io as pio
pio.templates.default = "plotly_dark"


class FreeEnergyLandscape:

    def __init__(self, hills: pd.DataFrame = None):

        if {'time', 'height'}.issubset(hills.columns) is False:
            raise ValueError("Make sure time and height is present")

        self.hills = (pd.DataFrame(hills)
                      .loc[:, ~hills.columns.str.startswith('sigma')]
                      .drop(columns=['biasf'])
                      .assign(time=lambda x: x['time']/1000)
                      .groupby('time', group_keys=False)
                      .apply(lambda x: x.assign(walker=lambda y: range(0, y.shape[0])))
                      )

        self.cvs = (self.hills
                    .drop(columns=['time', 'height', 'walker'])
                    .columns
                    .to_list()
                    )

    def get_hills_figures(self, round_time: int = None) -> dict[str, go.Figure]:
        """
        Function to get a dictionary of plotly figure objects summarising the dynamics and hills for each walker in a metadynamics simulation.
        :param round_time: bik the data into time bins with this number of decimal places. Useful for producing smaller figures
        :return:
        """
        if round_time:
            long_hills = (self.hills
                          .melt(value_vars=self.cvs+['height'], id_vars=['time', 'walker'])
                          .assign(time=lambda x: x['time'].round(round_time))
                          .groupby(['time', 'walker', 'variable'], group_keys=False)
                          .apply(lambda x: x.assign(value=lambda y: y['value'].mean()))
                          .drop_duplicates()
                          .groupby('walker')
                          )
        else:
            long_hills = (self.hills
                          .melt(value_vars=self.cvs+['height'], id_vars=['time', 'walker'])
                          .groupby('walker')
                          )

        figs = {}
        for name, df in long_hills:
            figure = px.line(df, x='time', y='value', facet_row='variable', height=800, width=1200, facet_row_spacing=0.04)
            figure.update_yaxes(matches=None, showgrid=False, title_text='')
            figure.update_xaxes(showgrid=False)
            figure.update_traces(line_color='#2596be', line_width=2)
            figs[name] = figure

        return figs
