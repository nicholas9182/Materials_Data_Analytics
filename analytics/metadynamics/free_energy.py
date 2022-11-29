import pandas as pd
import plotly.graph_objects as go
import plotly.express as px
import plumed as pl
from visualisation.themes import custom_dark_template
pd.set_option('mode.chained_assignment', None)


class MetaTrajectory:
    """
    Class to handle colvar files, which here are thought of as a metadynamics trajectory in CV space.
    """
    def __init__(self, colvar_file: str):

        self.colvar_file = colvar_file
        self.walker = int(colvar_file.split("/")[-1].split(".")[-1])
        self.data = pd.DataFrame(pl.read_as_pandas(colvar_file))

        if {'time', 'metad.bias', 'metad.rct'}.issubset(self.data) is False:
            raise ValueError("Make sure you have time, metad.bias and metad.rct in the colvar file")

        self.cvs = self.data.drop(columns=['time', 'metad.bias', 'metad.rct', 'metad.rbias']).columns.to_list()


class FreeEnergyLine:
    """
    Class to handle 1D fes files
    """
    def __init__(self, fes_file: str, time_data: dict[int, pd.DataFrame] = None):

        self.fes_file = fes_file
        self.data = pd.DataFrame(pl.read_as_pandas(fes_file))

        if {'projection'}.issubset(self.data.columns) is False:
            raise ValueError("Make sure the fes file has a projection column")

        self.cv = self.data.columns.values[0]
        self.time_data = time_data

    @classmethod
    def with_strides(cls, fes_files: list[str]):
        """
        function to build a fes line with time data, where the fes directories are given in a list
        :param fes_files: list of fes files to read
        :return: fes object
        """
        files = [f.split("/")[-1] for f in fes_files]
        time_stamps = [int(''.join(x for x in f if x.isdigit())) for f in files]
        fes_directories_dict = {time_stamps[i]: fes_files[i] for i in range(0, len(time_stamps))}

        time_data = {}
        for ts, fes_dir in fes_directories_dict.items():
            time_data[ts] = pd.DataFrame(pl.read_as_pandas(fes_dir))

        newest_data_dir = fes_directories_dict[max(fes_directories_dict)]
        return cls(fes_file=newest_data_dir, time_data=time_data)

    def set_datum(self, datum: float | int | tuple[float | int, float | int]):
        """
        Function to shift the fes line to set a new datum point. If a float is given, then the line will be shifted to give that x axis value an
        energy of 0.  If a tuple is given, then the fes will be shifted by the mean over that range.
        :param datum: either the point on the fes to set as the datum, or a range of the fes to set as the datum
        :param ymax: optional parameter to take off some high energy values from the line
        :return: self
        """
        if type(datum) == float or type(datum) == int:
            adjust_value = self.data.loc[(self.data[self.cv] - datum).abs().argsort()[:1], 'projection'].values[0]
            self.data['projection'] = self.data['projection'] - adjust_value
            if self.time_data:
                for key, item in self.time_data.items():
                    item['projection'] = item['projection'] - adjust_value
        elif type(datum) == tuple:
            adjust_value = self.data.loc[self.data[self.cv].between(min(datum), max(datum)), 'projection'].values.mean()
            self.data['projection'] = self.data['projection'] - adjust_value
            if self.time_data:
                for key, item in self.time_data.items():
                    item['projection'] = item['projection'] - adjust_value
        else:
            raise ValueError("Enter either a float or a tuple!")

        return self

    def get_time_difference(self, region_1: int | float, region_2: int | float) -> pd.DataFrame:

        time_data = []
        for key, df in self.time_data.items():

            if type(region_1) == int or type(region_1) == float:
                value_1 = df.loc[(df[self.cv] - region_1).abs().argsort()[:1], 'projection'].values[0]

            if type(region_2) == int or type(region_2) == float:
                value_2 = df.loc[(df[self.cv] - region_2).abs().argsort()[:1], 'projection'].values[0]

            difference = value_2 - value_1
            time_data.append(pd.DataFrame({'time_stamp': [key], 'energy_difference': [difference]}))

        time_data = pd.concat(time_data).sort_values('time_stamp')
        return time_data


class FreeEnergySpace:

    def __init__(self, hills_file: str):

        self.hills_file = hills_file
        hills = pd.DataFrame(pl.read_as_pandas(hills_file))

        if {'time', 'height'}.issubset(hills.columns) is False:
            raise ValueError("Make sure time and height is present")

        hills = (hills
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
        self.lines = {}
        self.trajectories = {}

    def add_metad_trajectory(self, meta_trajectory: MetaTrajectory, **kwargs):
        """
        function to add a metad trajectory to the landscape
        :param meta_trajectory: a metaD trajectory object to add to the landscape
        :return: the appended trajectories
        """
        if meta_trajectory not in self.trajectories:
            self.trajectories[meta_trajectory.walker] = meta_trajectory
        else:
            print(f"{meta_trajectory.colvar_file} is already in this landscape!")

        return self

    def add_line(self, fes_line: FreeEnergyLine):
        """
        function to add a free energy line to the landscape
        :param fes_line: the fes to add
        :return: the fes for the landscape
        """
        if fes_line not in self.lines:
            self.lines[fes_line.cv] = fes_line
        else:
            print(f"{fes_line.fes_file} is already in this landscape!")

        return self

    def get_long_hills(self, time_resolution: int = 6, height_power: float = 1):
        """
        Function to turn the hills into long format, and allow for time binning and height power conversion
        :param time_resolution: bin the data into time bins with this number of decimal places. Useful for producing a smaller long format hills
        dataframe.
        :param height_power: raise the height to the power of this so to see hills easier. Useful when plotting and you want to see the small hills.
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
                      )

        return long_hills

    def get_hills_figures(self, **kwargs) -> dict[int, go.Figure]:
        """
        Function to get a dictionary of plotly figure objects summarising the dynamics and hills for each walker in a metadynamics simulation.
        :return:
        """
        long_hills = self.get_long_hills(**kwargs).groupby('walker')
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
