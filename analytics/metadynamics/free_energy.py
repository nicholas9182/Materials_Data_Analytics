import pandas as pd
import numpy as np
import plotly.graph_objects as go
import plotly.express as px
from visualisation.themes import custom_dark_template
from analytics.laws_and_constants import boltzmann_energy_to_population, Kb, boltzmann_population_to_energy
pd.set_option('mode.chained_assignment', None)


class MetaTrajectory:
    """
    Class to handle colvar files, which here are thought of as a metadynamics trajectory in CV space.
    """
    def __init__(self, colvar_file: str, temperature: float = 298):

        self.data = (self._read_file(colvar_file)
                     .pipe(self._get_weights, temperature=temperature)
                     )

        self.walker = int(colvar_file.split("/")[-1].split(".")[-1])
        self.cvs = self.data.drop(columns=['time', 'bias', 'reweight_factor', 'reweight_bias']).columns.to_list()
        self.temperature = temperature

    @staticmethod
    def _read_file(file: str):
        """
        Function to read in colvar data, replacement for pl.read_as_pandas
        :param file: file to read in
        :return: data in that file in pandas format
        """
        col_names = open(file).readline().strip().split(" ")[2:]
        colvar = (pd.read_table(file, delim_whitespace=True, comment="#", names=col_names, dtype=np.float64)
                  .rename(columns={'metad.bias': 'bias', 'metad.rct': 'reweight_factor', 'metad.rbias': 'reweight_bias'})
                  )

        return colvar

    @staticmethod
    def _get_weights(data: pd.DataFrame, temperature: float = 298, y_col: str = 'reweight_bias',
                     y_col_out: str = 'weight') -> pd.DataFrame:

        new_col_args_1 = {y_col_out: lambda x: np.exp(x[y_col]/(Kb * temperature))}
        new_col_args_2 = {y_col_out: lambda x: x[y_col_out]/max(x[y_col_out])}
        data = (data
                .assign(**new_col_args_1)
                .assign(**new_col_args_2)
                )

        return data


class FreeEnergyLine:
    """
    Class to handle 1D fes files
    """
    def __init__(self, data: pd.DataFrame | dict[int | float, pd.DataFrame], temperature: float = 298):
        """
        current philosophy is to build a fes from data with an energy column and a column for the cv. Then use an alternate constructor to do it from
        a fes file. This makes it more general and allows us to calculate the fes from the hills ourselves at a later date
        :param data: data with the cv and the energy in two columns
        :param temperature: temperature at which the fes is defined
        """
        if type(data) == pd.DataFrame:
            if {'energy'}.issubset(data.columns) is False:
                raise ValueError("make sure there is an energy column in your dataframe")
            self.time_data = None
            self.data = data
        elif type(data) == dict:
            for _, value in data.items():
                if {'energy'}.issubset(value.columns) is False:
                    raise ValueError("make sure there is an energy column in each dataframe in your dict")
            self.time_data = data
            self.data = data[max(data)].copy()
        else:
            raise ValueError("fes_file must be a pd.Dataframe or a list[pd.Dataframe]")

        self.temperature = temperature
        self.cv = self.data.columns.values[0]

    @classmethod
    def from_plumed(cls, file: str | list[str], **kwargs):
        """
        alternate constructor to build the fes from a plumed file
        :param file: the file or list of files to make the plumed fes from. if list then it will make the time data
        :return: fes object
        """
        if type(file) == str:
            data = FreeEnergyLine._read_file(file, **kwargs)
        elif type(file) == list:
            individual_files = [f.split("/")[-1] for f in file]
            time_stamps = [int(''.join(x for x in f if x.isdigit())) for f in individual_files]
            data_frames = [FreeEnergyLine._read_file(f) for f in file]
            data = {time_stamps[i]: data_frames[i] for i in range(0, len(file))}
        else:
            raise ValueError("")

        return cls(data, **kwargs)

    @staticmethod
    def _read_file(file: str, temperature: float = 298):
        """
        Function to read in fes line data, replacement for pl.read_as_pandas. Does some useful other operations when reading in
        :param file: file to read in
        :param temperature: temperature of system
        :return: data in that file in pandas format
        """
        col_names = open(file).readline().strip().split(" ")[2:]
        data = (pd.read_table(file, delim_whitespace=True, comment="#", names=col_names, dtype=np.float64)
                .rename(columns={'projection': 'energy'})
                .pipe(boltzmann_energy_to_population, temperature=temperature, x_col=col_names[0])
                )

        return data

    def set_datum(self, datum: float | int | tuple[float | int, float | int]):
        """
        Function to shift the fes line to set a new datum point. If a float is given, then the line will be shifted to give that x-axis value an
        energy of 0.  If a tuple is given, then the fes will be shifted by the mean over that range.
        :param datum: either the point on the fes to set as the datum, or a range of the fes to set as the datum
        :return: self
        """
        if type(datum) == float or type(datum) == int:
            adjust_value = self.data.loc[(self.data[self.cv] - datum).abs().argsort()[:1], 'energy'].values[0]
            self.data['energy'] = self.data['energy'] - adjust_value
            if self.time_data is not None:
                for _, v in self.time_data.items():
                    v['energy'] = v['energy'] - adjust_value
        elif type(datum) == tuple:
            adjust_value = self.data.loc[self.data[self.cv].between(min(datum), max(datum)), 'energy'].values.mean()
            self.data['energy'] = self.data['energy'] - adjust_value
            if self.time_data is not None:
                for _, v in self.time_data.items():
                    v['energy'] = v['energy'] - adjust_value
        else:
            raise ValueError("Enter either a float or a tuple!")

        return self

    def get_time_difference(self, region_1: float | int | tuple[float | int, float | int],
                            region_2: float | int | tuple[float | int, float | int] = None) -> pd.DataFrame:
        """
        Function to get how the difference in energy between two points changes over time, or the energy of one point over time if region_2 is None.
        It can accept both numbers and tuples. If a tuple is given, it well take the mean of the CV over the interval given by the tupple.
        :param region_1: a point or region of the FES that you want to track as the first point
        :param region_2: a point or region of the FES that you want to track as the second point
        :return: pandas dataframe with the data
        """
        time_data = []

        for key, df in self.time_data.items():

            if type(region_1) == int or type(region_1) == float:
                value_1 = df.loc[(df[self.cv] - region_1).abs().argsort()[:1], 'energy'].values[0]
            elif type(region_1) == tuple:
                value_1 = df.loc[df[self.cv].between(min(region_1), max(region_1)), 'energy'].values.mean()
            else:
                raise ValueError("Use either a number or tuple of two numbers")

            if (type(region_2) == int or type(region_2) == float) and region_2 is not None:
                value_2 = df.loc[(df[self.cv] - region_2).abs().argsort()[:1], 'energy'].values[0]
            elif type(region_2) == tuple and region_2 is not None:
                value_2 = df.loc[df[self.cv].between(min(region_2), max(region_2)), 'energy'].values.mean()
            elif region_2 is None:
                value_2 = 0
            else:
                raise ValueError("Use either a number or tuple of two numbers")

            difference = value_2 - value_1
            time_data.append(pd.DataFrame({'time_stamp': [key], 'energy_difference': [difference]}))

        time_data = pd.concat(time_data).sort_values('time_stamp')
        return time_data


class FreeEnergySpace:

    def __init__(self, hills_file: str, temperature: float = 298):

        self.hills_file = hills_file
        self.hills = self._read_file(hills_file)
        self.n_walker = self.hills[self.hills['time'] == min(self.hills['time'])].shape[0]
        self.n_timesteps = self.hills[['time']].drop_duplicates().shape[0]
        self.max_time = self.hills['time'].max()
        self.dt = self.max_time/self.n_timesteps
        self.hills['walker'] = self.hills.groupby('time').cumcount()
        self.cvs = self.hills.drop(columns=['time', 'height', 'walker']).columns.to_list()
        self.temperature = temperature
        self.lines = {}
        self.trajectories = {}

    @staticmethod
    def _read_file(file: str):
        """
        Function to read in hills data
        :param file: file to read in
        :return: data in that file in pandas format
        """
        col_names = open(file).readline().strip().split(" ")[2:]
        data = pd.read_table(file, delim_whitespace=True, comment="#", names=col_names, dtype=np.float64)

        return (data
                .loc[:, ~data.columns.str.startswith('sigma')]
                .drop(columns=['biasf'])
                .assign(time=lambda x: x['time']/1000)
                )

    def add_metad_trajectory(self, meta_trajectory: MetaTrajectory):
        """
        function to add a metad trajectory to the landscape
        :param meta_trajectory: a metaD trajectory object to add to the landscape
        :return: the appended trajectories
        """
        if meta_trajectory not in self.trajectories.values():
            self.trajectories[meta_trajectory.walker] = meta_trajectory
        else:
            print(f"trajectory is already in this landscape!")

        return self

    def add_line(self, fes_line: FreeEnergyLine):
        """
        function to add a free energy line to the landscape
        :param fes_line: the fes to add
        :return: the fes for the landscape
        """
        self.lines[fes_line.cv] = fes_line
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
                         labels={'time': 'Time [ns]', 'height': 'Energy [kJ/mol]'})
        figure.update_traces(line=dict(width=1))

        return figure

    def get_line(self, cv, bins: int = 200):
        """
        Function to get a free energy line from a free energy space with meta trajectories in it.
        :param cv:
        :param bins:
        :return:
        """
        if cv not in self.trajectories[0].cvs:
            raise ValueError("Enter a collective variable that is in the trajectories!")

        data = [t.data for t in self.trajectories.values()]
        data = pd.concat(data)

        histogram = np.histogram(a=data[cv], bins=bins, weights=data['weight'])

        data = (pd.DataFrame(histogram, index=['population', cv])
                .transpose()
                .dropna()
                .assign(population=lambda c: c['population']/(np.trapz(y=c['population'], x=c[cv])))
                .pipe(boltzmann_population_to_energy, temperature=self.temperature)
                )

        line = FreeEnergyLine(data[[cv, 'energy', 'population']], temperature=self.temperature)
        return line
