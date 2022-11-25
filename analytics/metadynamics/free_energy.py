import pandas as pd
import plotly.graph_objects as go
import plotly.express as px
import plotly.io as pio
pio.templates.default = "plotly_dark"
pd.set_option('mode.chained_assignment', None)


class FreeEnergyLandscape:

    def __init__(self, hills: pd.DataFrame = None):

        if {'time', 'height'}.issubset(hills.columns) is False:
            raise ValueError("Make sure time and height is present")

        hills = (pd.DataFrame(hills)
                 .loc[:, ~hills.columns.str.startswith('sigma')]
                 .drop(columns=['biasf'])
                 .assign(time=lambda x: x['time']/1000)
                 )

        self.hills = hills[hills['time'] < max(hills['time'])]
        self.n_walker = self.hills[self.hills['time'] == min(self.hills['time'])].shape[0]
        self.n_timesteps = self.hills[['time']].drop_duplicates().shape[0]
        self.max_time = self.hills['time'].max()
        self.dt = self.max_time/self.n_timesteps
        self.hills['walker'] = [i for i in range(0, self.n_walker)] * self.n_timesteps
        self.cvs = self.hills.drop(columns=['time', 'height', 'walker']).columns.to_list()

    def get_hills_figures(self, time_resolution: int = 5, height_power: float = 1) -> dict[str, go.Figure]:
        """
        Function to get a dictionary of plotly figure objects summarising the dynamics and hills for each walker in a metadynamics simulation.
        :param time_resolution: bin the data into time bins with this number of decimal places. Useful for producing smaller figures
        :param height_power: raise the height to the power of this so to see hills easier
        :return:
        """
        height_label = 'height^' + str(height_power) if height_power != 1 else 'height'
        long_hills = self.hills.rename(columns={'height': height_label})
        long_hills[height_label] = long_hills[height_label].pow(height_power)
        long_hills = (long_hills
                      .melt(value_vars=self.cvs+[height_label], id_vars=['time', 'walker'])
                      .assign(time=lambda x: x['time'].round(time_resolution))
                      .groupby(['time', 'walker', 'variable'], group_keys=False)
                      .mean()
                      .reset_index()
                      .groupby('walker')
                      )

        figs = {}
        for name, df in long_hills:
            figure = px.line(df, x='time', y='value', facet_row='variable', height=600, width=1000, facet_row_spacing=0.02)
            figure.update_yaxes(matches=None, showgrid=False, title_text='', zerolinecolor='black ')
            figure.update_xaxes(showgrid=False)
            figure.update_layout(paper_bgcolor='black', plot_bgcolor="black")
            figure.update_traces(line_color='#f2f2f2', line_width=1)
            figure.for_each_annotation(lambda a: a.update(text=a.text.split("=")[1]))
            figs[name] = figure

        return figs

