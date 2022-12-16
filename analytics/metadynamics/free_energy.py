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
        self.cvs = self.data.drop(columns=['time', 'bias', 'reweight_factor', 'reweight_bias', 'weight']).columns.to_list()
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
                  .assign(time=lambda x: x['time'] / 1000)
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


class FreeEnergyShape:

    def __init__(self, data: pd.DataFrame | dict[int | float, pd.DataFrame], temperature: float = 298, dimension: int = None):
        """
        Current philosophy is now that there should be a super state with some general properties of free energy shapes.  Lines, surfaces and other
        shapes should inherit from this class, and then make changes depending on whether the shape has particular features
        :param data: the data needed to make the shape
        :param temperature: the temperature at which the shape is defined
        :param dimension: the dimension of the shape
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
        self.cvs = self.data.columns.values.tolist()[:dimension]
        self.dimension = dimension

    @classmethod
    def from_plumed(cls, file: str | list[str], **kwargs):
        """
        alternate constructor to build the fes from a plumed file
        :param file: the file or list of files to make the plumed fes from. if list then it will make the time data
        :return: fes object
        """
        if type(file) == str:
            data = cls._read_file(file, **kwargs)
        elif type(file) == list:
            individual_files = [f.split("/")[-1] for f in file]
            time_stamps = [int(''.join(x for x in f if x.isdigit())) for f in individual_files]
            data_frames = [cls._read_file(f) for f in file]
            data = {time_stamps[i]: data_frames[i] for i in range(0, len(file))}
        else:
            raise ValueError("")

        return cls(data, **kwargs)

    @staticmethod
    def _get_nearest_value(data: pd.DataFrame, ref_coordinate: dict[str, float | int], val_col: str) -> float:
        """
        Function to read a dataframe and get the value in val_col in the row where ref_col is closest to value using a pythagorean distance
        :param data: data
        :param ref_coordinate: dict with the column as the key and the value as the value
        :param val_col: column from which to get the value
        :return: value
        """
        new_cols = []

        for key, value in ref_coordinate.items():
            new_col = key + "_distance"
            new_cols.append(new_col)
            data[new_col] = (data[key].abs() - value)**2

        data["total_distance"] = 0
        for c in new_cols:
            data["total_distance"] = data["total_distance"] + data[c]

        sorted_data = data.sort_values('total_distance').reset_index(drop=True)
        closest_value = sorted_data.loc[0, val_col]
        return closest_value

    @staticmethod
    def _get_mean_in_range(data: pd.DataFrame, ref_col, val_col, area: tuple[int | float, int | float]):
        """
        Get the mean value in val_col over a range in ref_col, assumes ordered data
        :param data: data
        :param ref_col: the column over which you take the range
        :param val_col: the column from which you want the mean
        :param area: tuple specifying the range
        :return: value
        """
        column_filtered = data[ref_col].between(min(area), max(area))
        data_filtered = data.loc[column_filtered, val_col]
        return data_filtered.values.mean()

    def set_datum(self, datum: dict[str, float | int | tuple[float | int, float | int]]):
        """
        Function to shift the fes line to set a new datum point. If a float is given, then the line will be shifted to give that x-axis value an
        energy of 0.  If a tuple is given, then the fes will be shifted by the mean over that range.
        :param datum: either the point on the fes to set as the datum, or a range of the fes to set as the datum
        :return: self
        """
        if type(datum) != dict:
            raise TypeError("Datum must be a dictionary with the cv and value")

        for cv, d in datum.items():
            if cv not in self.cvs:
                raise ValueError("The keys for the datum dictionary need to be cvs!")

        if type(datum[self.cvs[0]]) == float or type(datum[self.cvs[0]]) == int:
            adjust_value = self._get_nearest_value(self.data, datum, 'energy')
            self.data['energy'] = self.data['energy'] - adjust_value
            if self.time_data is not None:
                for _, v in self.time_data.items():
                    adjust_value = self._get_nearest_value(v, datum, 'energy')
                    v['energy'] = v['energy'] - adjust_value
        elif type(datum[self.cvs[0]]) == tuple:
            adjust_value = self._get_mean_in_range(self.data, self.cvs[0], 'energy', datum[self.cvs[0]])
            self.data['energy'] = self.data['energy'] - adjust_value
            if self.time_data is not None:
                for _, v in self.time_data.items():
                    adjust_value = self._get_mean_in_range(v, self.cvs[0], 'energy', datum[self.cvs[0]])
                    v['energy'] = v['energy'] - adjust_value
        else:
            raise ValueError("Enter either a float or a tuple!")

        return self

    @staticmethod
    def _read_file(file, **kwargs):
        pass


class FreeEnergyLine(FreeEnergyShape):

    def __init__(self, data: pd.DataFrame | dict[int | float, pd.DataFrame], temperature: float = 298):

        super().__init__(data, temperature, dimension=1)

        cv_min = self.data[self.cvs[0]].min()
        cv_max = self.data[self.cvs[0]].max()

        if self.data.shape[0] > 2:
            lower = cv_min + (cv_max - cv_min) / 6
            upper = cv_min + (5 * (cv_max - cv_min)) / 6
            self.set_datum({self.cvs[0]: (lower, upper)})
        else:
            self.set_datum({self.cvs[0]: (cv_min, cv_max)})

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

    def get_time_difference(self, region_1: float | int | tuple[float | int, float | int],
                            region_2: float | int | tuple[float | int, float | int] = None) -> pd.DataFrame:
        """
        Function to get how the difference in energy between two points changes over time, or the energy of one point over time if region_2 is None.
        It can accept both numbers and tuples. If a tuple is given, it will take the mean of the CV over the interval given by the tuple.
        :param region_1: a point or region of the FES that you want to track as the first point
        :param region_2: a point or region of the FES that you want to track as the second point
        :return: pandas dataframe with the data
        """
        time_data = []

        for key, df in self.time_data.items():

            if type(region_1) == int or type(region_1) == float:
                value_1 = df.loc[(df[self.cvs[0]] - region_1).abs().argsort()[:1], 'energy'].values[0]
            elif type(region_1) == tuple:
                value_1 = df.loc[df[self.cvs[0]].between(min(region_1), max(region_1)), 'energy'].values.mean()
            else:
                raise ValueError("Use either a number or tuple of two numbers")

            if (type(region_2) == int or type(region_2) == float) and region_2 is not None:
                value_2 = df.loc[(df[self.cvs[0]] - region_2).abs().argsort()[:1], 'energy'].values[0]
            elif type(region_2) == tuple and region_2 is not None:
                value_2 = df.loc[df[self.cvs[0]].between(min(region_2), max(region_2)), 'energy'].values.mean()
            elif region_2 is None:
                value_2 = 0
            else:
                raise ValueError("Use either a number or tuple of two numbers")

            difference = value_2 - value_1
            time_data.append(pd.DataFrame({'time_stamp': [key], 'energy_difference': [difference]}))

        time_data = pd.concat(time_data).sort_values('time_stamp')
        return time_data

    def set_errors_from_time_dynamics(self, n_timestamps: int, bins: int = 200):
        """
        Function to get data and errors from considering the time dynamics of the FES
        :param n_timestamps: How many past FES time stamps to look at. Consider plotting the value of the minima as a function of time to see what
        an appropriate value is for this
        :param bins: Number of data points to have on your FES
        :return: dataframe with the errors.
        """
        if self.time_data is None:
            raise ValueError("You need time data to use this function")

        data = []
        min_timestamp = max(self.time_data) - n_timestamps

        for key, value in self.time_data.items():
            value['timestamp'] = key
            data.append(value)

        data = (pd.concat(data)
                .query('timestamp > @min_timestamp')
                .assign(bin=lambda x: pd.cut(x[self.cvs[0]], bins))
                )

        binned_data = pd.DataFrame({
            self.cvs[0]: data.groupby('bin').mean()[self.cvs[0]],
            'energy': data.groupby('bin').mean()['energy'],
            'energy_err': data.groupby('bin').std()['energy'],
            'population': data.groupby('bin').mean()['population'],
            'population_err': data.groupby('bin').std()['population']/np.sqrt(n_timestamps)
        }).dropna()

        self.data = binned_data

        return self


class FreeEnergySurface(FreeEnergyShape):

    def __init__(self, data: pd.DataFrame | dict[int | float, pd.DataFrame], temperature: float = 298):
        super().__init__(data, temperature, dimension=2)

    @staticmethod
    def _read_file(file: str, temperature: float = 298):
        """
        Function to read in fes surface data, replacement for pl.read_as_pandas. Does some useful other operations when reading in
        :param file: file to read in
        :param temperature: temperature of system
        :return: data in that file in pandas format
        """
        col_names = open(file).readline().strip().split(" ")[2:]
        drop_cols = [c for c in col_names if 'der_' in c]

        data = (pd.read_table(file, delim_whitespace=True, comment="#", names=col_names, dtype=np.float64)
                .drop(columns=drop_cols)
                .rename(columns={'file.free': 'energy'})
                .pipe(boltzmann_energy_to_population, temperature=temperature, x_col=col_names[0])
                )

        return data


class FreeEnergySpace:

    def __init__(self, hills_file: str, temperature: float = 298):

        self.hills_path = hills_file
        self.hills, self.sigmas = self._read_file(hills_file)
        self.n_walker = self.hills[self.hills['time'] == min(self.hills['time'])].shape[0]
        self.n_timesteps = self.hills[['time']].drop_duplicates().shape[0]
        self.max_time = self.hills['time'].max()
        self.dt = self.max_time/self.n_timesteps
        self.hills['walker'] = self.hills.groupby('time').cumcount()
        self.cvs = self.hills.drop(columns=['time', 'height', 'walker']).columns.to_list()
        self.temperature = temperature
        self.lines = {}
        self.surfaces = []
        self.trajectories = {}

    @staticmethod
    def _read_file(file: str):
        """
        Function to read in hills data
        :param file: file to read in
        :return: data in that file in pandas format
        """
        col_names = open(file).readline().strip().split(" ")[2:]
        sigmas = [col for col in col_names if col.split("_")[0] == 'sigma']
        data = pd.read_table(file, delim_whitespace=True, comment="#", names=col_names, dtype=np.float64)
        sigmas = {s.split("_")[1]: data.loc[0, s] for s in sigmas}

        data = (data
                .loc[:, ~data.columns.str.startswith('sigma')]
                .drop(columns=['biasf'])
                .assign(time=lambda x: x['time'] / 1000)
                )

        return data, sigmas

    def add_metad_trajectory(self, meta_trajectory: MetaTrajectory):
        """
        function to add a metad trajectory to the landscape
        :param meta_trajectory: a metaD trajectory object to add to the landscape
        :return: the appended trajectories
        """
        if meta_trajectory not in self.trajectories.values():
            self.trajectories[meta_trajectory.walker] = meta_trajectory
        else:
            print(f"trajectory is already in this space!")

        return self

    def add_line(self, line: FreeEnergyLine):
        """
        function to add a free energy line to the landscape
        :param line: the fes to add
        :return: the fes for the landscape
        """
        self.lines[line.cvs[0]] = line
        return self

    def add_surface(self, surface: FreeEnergySurface):
        """
        function to add a free energy surface to the landscape
        :param surface: the fes to add
        :return: the fes for the landscape
        """
        if surface not in self.surfaces:
            self.surfaces.append(surface)
        else:
            raise ValueError("This surface is already in the space")
        return self

    def get_long_hills(self, time_resolution: int = 6, height_power: float = 1):
        """
        Function to turn the hills into long format, and allow for time binning and height power conversion
        :param time_resolution: bin the data into time bins with this number of decimal places. Useful for producing a smaller long format hills
        dataframe.
        :param height_power: raise the height to the power of this so to see hills easier. Useful when plotting, and you want to see the small hills.
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

    def get_reweighted_line(self, cv, bins: int | list[int | float] = 200):
        """
        Function to get a free energy line from a free energy space with meta trajectories in it, using weighted histogram
        analysis
        :param cv: the cv in which to get the reweight
        :param bins: number of bins, or a list with the bin boundaries
        :return:
        """
        data = []

        for w, t in self.trajectories.items():
            if cv in t.cvs:
                data.append(t.data)

        if not data:
            raise ValueError("no trajectories in this space have that CV")

        data = pd.concat(data)

        histogram = np.histogram(a=data[cv], bins=bins, weights=data['weight'], density=True)

        data = (pd.DataFrame(histogram, index=['population', cv])
                .transpose()
                .dropna()
                .pipe(boltzmann_population_to_energy, temperature=self.temperature)
                )

        line = FreeEnergyLine(data[[cv, 'energy', 'population']], temperature=self.temperature)
        return line
