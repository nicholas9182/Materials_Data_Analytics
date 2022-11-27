import pandas as pd
import plotly.graph_objects as go
import plotly.express as px
from visualisation.themes import custom_dark_template
pd.set_option('mode.chained_assignment', None)


class metad_trajectory:

    def __init__(self, colvar: pd.DataFrame, walker_num):

        if {'time', 'metad.bias', 'metad.rct'}.issubset(colvar.columns) is False:
            raise ValueError("Make sure you have time, metad.bias and metad.rct in the colvar file")

        self.trajectory = pd.DataFrame(colvar)
        self.cvs = self.trajectory.drop(columns=['time', 'metad.bias', 'metad.rct', 'metad.rbias']).columns.to_list
        self.walker = walker_num


class FreeEnergyLine:

    def __init__(self, fes: pd.DataFrame, timestamp: int = None):

        if {'projection'}.issubset(fes.columns) is False:
            raise ValueError("Make sure the fes file has a projection column")

        self.data = pd.DataFrame(fes)
        self.cv = fes.columns.values[0]
        self.timestamp = timestamp


class FreeEnergyLandscape:

    def __init__(self, hills: pd.DataFrame = None):

        if {'time', 'height'}.issubset(hills.columns) is False:
            raise ValueError("Make sure time and height is present")

        hills = (pd.DataFrame(hills)
                 .loc[:, ~hills.columns.str.startswith('sigma')]
                 .drop(columns=['biasf'])
                 .assign(time=lambda x: x['time']/1000)
                 )

        self.hills = hills[hills['time'] < hills['time'].iloc[-1]]
        self.n_walker = self.hills[self.hills['time'] == min(self.hills['time'])].shape[0]
        self.n_timesteps = self.hills[['time']].drop_duplicates().shape[0]
        self.max_time = self.hills['time'].max()
        self.dt = self.max_time/self.n_timesteps
        self.hills['walker'] = self.hills.groupby('time').cumcount()
        self.cvs = self.hills.drop(columns=['time', 'height', 'walker']).columns.to_list()
        self.fes = []
        self.trajectories = []

    def add_metad_trajectory(self, colvar_data: pd.DataFrame, **kwargs):
        """
        function to add a metad trajectory to the landscape
        :param colvar_data: dataframe with the colvar data, read in using pl.read_as_pandas
        :param kwargs: other kwargs, mostly for the walker number
        :return:
        """
        trajectory = metad_trajectory(colvar_data, **kwargs)
        self.trajectories.append(trajectory)
        return self.trajectories

    def add_fes_line(self, fes_df: pd.DataFrame, **kwargs):
        """
        function to add a free energy line to the landscape
        :param fes_df: dataframe with the fes data, read in using pl.read_as_pandas
        :param kwargs: other kwargs, mostly for the time stamp
        :return: the fes for the landscape
        """
        fes = FreeEnergyLine(fes_df, **kwargs)
        self.fes.append(fes)
        return self.fes

    def get_long_hills(self, time_resolution: int = 5, height_power: float = 1):
        """
        Function to turn the hills into long format, and allow for time binning and height power conversion
        :param time_resolution: bin the data into time bins with this number of decimal places. Useful for producing smaller figures
        :param height_power: raise the height to the power of this so to see hills easier
        :return:
        """
        height_label = 'height^' + str(height_power) if height_power != 1 else 'height'
        long_hills = self.hills.rename(columns={'height': height_label})
        long_hills[height_label] = long_hills[height_label].pow(height_power)
        long_hills = (long_hills
                      .melt(value_vars=self.cvs + [height_label], id_vars=['time', 'walker'])
                      .assign(time=lambda x: x['time'].round(time_resolution))
                      .groupby(['time', 'walker', 'variable'], group_keys=False)
                      .mean()
                      .reset_index()
                      .groupby('walker')
                      )
        return long_hills

    def get_hills_figures(self, **kwargs) -> dict[str, go.Figure]:
        """
        Function to get a dictionary of plotly figure objects summarising the dynamics and hills for each walker in a metadynamics simulation.
        :return:
        """
        long_hills = self.get_long_hills(**kwargs)
        figs = {}
        for name, df in long_hills:
            figure = px.line(df, x='time', y='value', facet_row='variable', labels={'time': 'Time [ns]'}, template=custom_dark_template)
            figure.update_traces(line=dict(width=1))
            figure.update_yaxes(title=None, matches=None)
            figure.for_each_annotation(lambda a: a.update(text=a.text.split("=")[1]))
            figs[name] = figure
        return figs

    def get_hills_average_across_walkers(self, time_resolution: int = 5):
        """
        function to get the average hill height, averaged across the walkers.
        :return:
        """
        av_hills = (self.hills
                    .assign(time=lambda x: x['time'].round(time_resolution))
                    .groupby(['time'])
                    .mean()
                    .reset_index()
                    )

        return av_hills[['time', 'height']]

    def get_average_hills_figure(self, **kwargs):
        """
        function to get figure of average hills across walkers
        :return:
        """
        av_hills = self.get_hills_average_across_walkers(**kwargs)
        figure = px.line(av_hills, x='time', y='height', log_y=True, template=custom_dark_template,
                         labels={'time': 'Time [ns]', 'height': 'Energy [kJ/mol]'}
                         )
        figure.update_traces(line=dict(width=1))
        return figure

    def get_hills_max_across_walkers(self, time_resolution: int = 5):
        """
        function to get the average hill height, averaged across the walkers.
        :return:
        """
        max_hills = (self.hills
                     .assign(time=lambda x: x['time'].round(time_resolution))
                     .groupby(['time'])
                     .max()
                     .reset_index()
                     )

        return max_hills[['time', 'height']]

    def get_max_hills_figure(self, **kwargs):
        """
        function to get figure of average hills across walkers
        :return:
        """
        max_hills = self.get_hills_max_across_walkers(**kwargs)
        figure = px.line(max_hills, x='time', y='height', log_y=True, template=custom_dark_template,
                         labels={'time': 'Time [ns]', 'height': 'Energy [kJ/mol]'}
                         )
        figure.update_traces(line=dict(width=1))
        return figure
