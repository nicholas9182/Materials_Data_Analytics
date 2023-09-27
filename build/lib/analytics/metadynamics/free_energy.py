from __future__ import annotations
import pandas as pd
import numpy as np
import os
import plotly.graph_objects as go
import plotly.express as px
from visualisation.themes import custom_dark_template
from analytics.laws_and_constants import boltzmann_energy_to_population, Kb, boltzmann_population_to_energy
pd.set_option('mode.chained_assignment', None)


class MetaTrajectory:
    """
    Class to handle colvar files, which here are thought of as a metadynamics trajectory in CV space.
    """
    def __init__(self, colvar_file: str, temperature: float = 298, metadata: dict = None):

        data, self._opes = self._read_file(colvar_file)
        self._data = data.pipe(self._get_weights, temperature=temperature)
        self.walker = int(colvar_file.split("/")[-1].split(".")[-1])
        self.cvs = (self._data.drop(columns=['time', 'bias', 'reweight_factor', 'reweight_bias', 'weight', 'zed',
                                             'neff', 'nker'], errors='ignore')
                    .columns
                    .to_list()
                    )
        self.temperature = temperature
        self._metadata = metadata

    @property
    def opes(self):
        return self._opes

    @staticmethod
    def _read_file(file: str):
        """
        Function to read in colvar _data, replacement for pl.read_as_pandas
        :param file: file to read in
        :return: _data in that file in pandas format
        """
        col_file = open(file)
        col_names = col_file.readline().strip().split(" ")[2:]
        col_file.close()
        opes = True if 'opes.bias' in col_names else False

        # TODO: Check that opes.bias is the right bias to use for reweighting!
        colvar = (pd.read_table(file, delim_whitespace=True, comment="#", names=col_names, dtype=np.float64)
                  .rename(columns={'metad.bias': 'bias', 'metad.rct': 'reweight_factor', 'metad.rbias':
                                   'reweight_bias', 'opes.bias': 'reweight_bias', 'opes.rct': 'reweight_factor',
                                   'opes.zed': 'zed', 'opes.neff': 'neff', 'opes.nker': 'nker'})
                  .assign(time=lambda x: x['time'] / 1000)
                  )

        return colvar, opes

    @staticmethod
    def _get_weights(data: pd.DataFrame, temperature: float = 298, y_col: str = 'reweight_bias',
                     y_col_out: str = 'weight') -> pd.DataFrame:
        """
        Function to get the weights for each from the data obtained in the colvar file.
        :param data:
        :param temperature:
        :param y_col:
        :param y_col_out:
        :return:
        """
        new_col_args_1 = {y_col_out: lambda x: np.exp(x[y_col]/(Kb * temperature))}
        new_col_args_2 = {y_col_out: lambda x: x[y_col_out]/max(x[y_col_out])}
        data = (data
                .assign(**new_col_args_1)
                .assign(**new_col_args_2)
                )

        return data

    def get_data(self, with_metadata: bool = False, time_resolution: int = None):
        """
        function to get the _data from a free energy shape
        :param with_metadata: print the data with the metadata?
        :param time_resolution: reduce the size of the data frame by reducing the time resolution
        :return:
        """
        data = self._data.copy()

        if with_metadata:
            data['temperature'] = self.temperature
            data['walker'] = self.walker

            if self._metadata:
                for key, value in self._metadata.items():
                    data[key] = value

        if time_resolution:
            data = (data
                    .assign(time=lambda x: x['time'].round(time_resolution))
                    .groupby(['time'])
                    .mean()
                    .reset_index()
                    .assign(walker=lambda x: x['walker'].astype(int))
                    )

        return data


class FreeEnergyShape:

    def __init__(self, data: pd.DataFrame | dict[int | float], temperature: float = 298, dimension: int = None,
                 metadata: dict = None):
        """
        Current philosophy is now that there should be a super state with some general properties of free energy shapes.
        Lines, surfaces and other shapes should inherit from this class, and then make changes depending on whether the
        shape has particular features
        :param data: the _data needed to make the shape
        :param temperature: the temperature at which the shape is defined
        :param dimension: the dimension of the shape
        """
        if type(data) == pd.DataFrame:
            if {'energy'}.issubset(data.columns) is False:
                raise ValueError("make sure there is an energy column in your dataframe")
            self._time_data = None
            self._data = data
        elif type(data) == dict:
            for _, value in data.items():
                if {'energy'}.issubset(value.columns) is False:
                    raise ValueError("make sure there is an energy column in each dataframe in your dict")
            self._time_data = data
            self._data = data[max(data)].copy()
        else:
            raise ValueError("fes_file must be a pd.Dataframe or a list[pd.Dataframe]")

        self.temperature = temperature
        self.cvs = self._data.columns.values.tolist()[:dimension]
        self.dimension = dimension
        self._metadata = metadata

    @classmethod
    def from_plumed(cls, file: str | list[str], **kwargs):
        """
        alternate constructor to build the fes from a plumed file
        :param file: the file or list of files to make the plumed fes from. if list then it will make the time _data
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
        Function to read a dataframe and get the value in val_col in the row where ref_col is closest to value using a
        pythagorean distance
        :param data: _data
        :param ref_coordinate: dict with the column as the key and the value as the value
        :param val_col: column from which to get the return
        :return: value
        """
        new_cols = []
        data = data.copy()

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
        Get the mean value in val_col over a range in ref_col, assumes ordered _data
        :param data: _data
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
        Function to shift the fes line to set a new datum point. If a float is given, then the line will be shifted to
        give that x-axis value an energy of 0.  If a tuple is given, then the fes will be shifted by the mean over that
        range.
        :param datum: either the point on the fes to set as the datum, or a range of the fes to set as the datum
        :return: self
        """
        if type(datum) != dict:
            raise TypeError("Datum must be a dictionary with the cv and value")

        for cv, d in datum.items():
            if cv not in self.cvs:
                raise ValueError("The keys for the datum dictionary need to be cvs!")

        if type(datum[self.cvs[0]]) == float or type(datum[self.cvs[0]]) == int:
            adjust_value = self._get_nearest_value(self._data, datum, 'energy')
            self._data['energy'] = self._data['energy'] - adjust_value
            if self._time_data is not None:
                for _, v in self._time_data.items():
                    adjust_value = self._get_nearest_value(v, datum, 'energy')
                    v['energy'] = v['energy'] - adjust_value
        elif type(datum[self.cvs[0]]) == tuple:
            adjust_value = self._get_mean_in_range(self._data, self.cvs[0], 'energy', datum[self.cvs[0]])
            self._data['energy'] = self._data['energy'] - adjust_value
            if self._time_data is not None:
                for _, v in self._time_data.items():
                    adjust_value = self._get_mean_in_range(v, self.cvs[0], 'energy', datum[self.cvs[0]])
                    v['energy'] = v['energy'] - adjust_value
        else:
            raise ValueError("Enter either a float or a tuple!")

        return self

    @staticmethod
    def _read_file(file, **kwargs):
        pass

    def get_data(self, with_metadata: bool = False, with_timedata: bool = False):
        """
        function to get the _data from a free energy shape
        :param with_metadata: also return the metadata as columns
        :param with_timedata: return the time data with the time stamp as a column
        :return:
        """
        if with_timedata:
            data = []
            for ts, d in self._time_data.items():
                d['timestamp'] = ts
                data.append(d)
            data = pd.concat(data)
        else:
            data = self._data.copy()

        if with_metadata:
            data['temperature'] = self.temperature

            if self._metadata:
                for key, value in self._metadata.items():
                    data[key] = value

        return data


class FreeEnergyLine(FreeEnergyShape):

    def __init__(self, data: pd.DataFrame | dict[int | float, pd.DataFrame], temperature: float = 298,
                 metadata: dict = None):

        super().__init__(data, temperature, metadata=metadata, dimension=1)

        cv_min = self._data[self.cvs[0]].min()
        cv_max = self._data[self.cvs[0]].max()

        if self._data.shape[0] > 2:
            lower = cv_min + (cv_max - cv_min) / 6
            upper = cv_min + (5 * (cv_max - cv_min)) / 6
            self.set_datum({self.cvs[0]: (lower, upper)})
        else:
            self.set_datum({self.cvs[0]: (cv_min, cv_max)})

    @staticmethod
    def _read_file(file: str, temperature: float = 298):
        """
        Function to read in fes line _data, replacement for pl.read_as_pandas. Does some useful other operations when
        reading in
        :param file: file to read in
        :param temperature: temperature of system
        :return: _data in that file in pandas format
        """
        col_file = open(file)
        col_names = col_file.readline().strip().split(" ")[2:]
        col_file.close()
        data = (pd.read_table(file, delim_whitespace=True, comment="#", names=col_names, dtype=np.float64)
                .rename(columns={'projection': 'energy'})
                .pipe(boltzmann_energy_to_population, temperature=temperature, x_col=col_names[0])
                )

        return data

    def get_time_difference(self, region_1: float | int | tuple[float | int, float | int],
                            region_2: float | int | tuple[float | int, float | int] = None,
                            with_metadata: bool = False) -> pd.DataFrame:
        """
        Function to get how the difference in energy between two points changes over time, or the energy of one point
        over time if region_2 is None. It can accept both numbers and tuples. If a tuple is given, it will take the mean
        of the CV over the interval given by the tuple.
        :param region_1: a point or region of the FES that you want to track as the first point
        :param region_2: a point or region of the FES that you want to track as the second point
        :param with_metadata: whether to return _data with the line _metadata
        :return: pandas dataframe with the _data
        """
        time_data = []

        for key, df in self._time_data.items():

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

        if with_metadata:
            time_data['temperature'] = self.temperature
            for key, value in self._metadata.items():
                time_data[key] = value

        return time_data

    def set_errors_from_time_dynamics(self, n_timestamps: int, bins: int = 200):
        """
        Function to get _data and errors from considering the time dynamics of the FES
        :param n_timestamps: How many past FES time stamps to look at. Consider plotting the value of the minima as a
        function of time to see what
        an appropriate value is for this
        :param bins: Number of _data points to have on your FES
        :return: dataframe with the errors.
        """
        if self._time_data is None:
            raise ValueError("You need time _data to use this function")

        data = []
        min_timestamp = max(self._time_data) - n_timestamps
        for key, value in self._time_data.items():
            new_data = value.copy()
            new_data['timestamp'] = key
            data.append(new_data)

        data = (pd.concat(data)
                .query('timestamp > @min_timestamp')
                .assign(bin=lambda x: pd.cut(x[self.cvs[0]], bins))
                )

        binned_data = pd.DataFrame({
            self.cvs[0]: data.groupby('bin', observed=True).mean()[self.cvs[0]],
            'energy': data.groupby('bin', observed=True).mean()['energy'],
            'energy_err': data.groupby('bin', observed=True).std()['energy'],
            'population': data.groupby('bin', observed=True).mean()['population'],
            'population_err': data.groupby('bin', observed=True).std()['population']/np.sqrt(n_timestamps)
        }).dropna()

        self._data = binned_data

        return self


class FreeEnergySurface(FreeEnergyShape):

    def __init__(self, data: pd.DataFrame | dict[int | float, pd.DataFrame], temperature: float = 298,
                 metadata: dict = None):

        super().__init__(data, temperature, dimension=2, metadata=metadata)

    @staticmethod
    def _read_file(file: str, temperature: float = 298):
        """
        Function to read in fes surface _data, replacement for pl.read_as_pandas. Does some useful other operations when
        reading in
        :param file: file to read in
        :param temperature: temperature of system
        :return: _data in that file in pandas format
        """
        col_file = open(file)
        col_names = col_file.readline().strip().split(" ")[2:]
        col_file.close()
        drop_cols = [c for c in col_names if 'der_' in c]

        data = (pd.read_table(file, delim_whitespace=True, comment="#", names=col_names, dtype=np.float64)
                .drop(columns=drop_cols)
                .rename(columns={'file.free': 'energy'})
                .pipe(boltzmann_energy_to_population, temperature=temperature, x_col=col_names[0])
                )

        return data


class FreeEnergySpace:

    def __init__(self, hills_file: str = None, temperature: float = 298, metadata: dict = None):

        if hills_file is not None:
            self._hills, self.sigmas = self._read_file(hills_file)
            self.n_walker = self._hills[self._hills['time'] == min(self._hills['time'])].shape[0]
            self.n_timesteps = self._hills[['time']].drop_duplicates().shape[0]
            self.max_time = self._hills['time'].max()
            self.dt = self.max_time/self.n_timesteps
            self.cvs = (self
                        ._hills.drop(columns=['time', 'height', 'walker', 'logweight'], errors='ignore')
                        .columns
                        .to_list()
                        )
            self._opes = False if 'logweight' not in self._hills.columns.to_list() else True
        elif hills_file is None:
            self.n_walker = 0
        self.temperature = temperature
        self.lines = {}
        self.surfaces = []
        self.trajectories = {}
        self._metadata = metadata

    @classmethod
    def from_standard_directory(cls, standard_dir, colvar_string_matcher: str = "COLVAR_REWEIGHT.", **kwargs):
        """
        alternate constructor to make a free energy space from a standard metadynamics directory. In this directory,
        the free energy lines and surfaces are held in folders called FES_* . The reweight data is held in COLVAR files
        called COLVAR_REWEIGHT.* .
        :param standard_dir: The directory with the plumed/gromacs files
        :param colvar_string_matcher: the string that matches to the colvar files names
        :return: a populated FreeEnergySpace
        """
        space = cls(**kwargs)

        for f in os.scandir(standard_dir):

            if (f.is_dir() and f.path.split("/")[-1].split("_")[0] == "FES"
                    and len(f.path.split("/")[-1].split("_")) == 2):
                path = f.path + "/"
                print(f"Adding a free energy line from files in {path}")
                files = [path + d for d in os.listdir(path)]
                files = files[0] if len(files) == 1 else files
                line = FreeEnergyLine.from_plumed(files)
                space.add_line(line)

        for f in os.scandir(standard_dir):
            if (f.is_dir() and f.path.split("/")[-1].split("_")[0] == "FES"
                    and len(f.path.split("/")[-1].split("_")) == 3):
                path = f.path + "/"
                print(f"Adding a free energy surface from files in {path}")
                files = [path + d for d in os.listdir(path)]
                files = files[0] if len(files) == 1 else files
                surface = FreeEnergySurface.from_plumed(files)
                space.add_surface(surface)

        for f in [standard_dir + "/" + f for f in os.listdir(standard_dir)
                  if colvar_string_matcher in f and 'bck' not in f]:
            file = f.split("/")[-1]
            print(f"Adding {file} as a metaD trajectory")
            traj = MetaTrajectory(f)
            space.add_metad_trajectory(traj)

        return space

    @property
    def opes(self):
        if hasattr(self, "_opes"):
            return self._opes
        else:
            return None

    @staticmethod
    def _read_file(file: str):
        """
        Function to read in _hills _data
        :param file: file to read in
        :return: _data in that file in pandas format
        """
        col_file = open(file)
        col_names = col_file.readline().strip().split(" ")[2:]
        col_file.close()
        sigmas = [col for col in col_names if col.split("_")[0] == 'sigma']
        data = pd.read_table(file, delim_whitespace=True, comment="#", names=col_names, dtype=np.float64)
        sigmas = {s.split("_")[1]: data.loc[0, s] for s in sigmas}

        data = (data
                .loc[:, ~data.columns.str.startswith('sigma')]
                .drop(columns=['biasf'], errors='ignore')
                .assign(time=lambda x: x['time'] / 1000)
                .assign(walker=lambda x: x.groupby('time').cumcount())
                )

        return data, sigmas

    def add_metad_trajectory(self, meta_trajectory: MetaTrajectory):
        """
        function to add a metad trajectory to the landscape
        :param meta_trajectory: a metaD trajectory object to add to the landscape
        :return: the appended trajectories
        """
        if not meta_trajectory.temperature == self.temperature:
            raise ValueError("Your trajectory has a different temperature to your space!")
        meta_trajectory._metadata = self._metadata
        self.trajectories[meta_trajectory.walker] = meta_trajectory
        opes_before = self._opes if hasattr(self, "_opes") else None
        self._opes = meta_trajectory.opes
        self.n_walker = self.n_walker if hasattr(self, "_hills") else self.n_walker + 1

        if opes_before is False and self._opes is True:
            raise ValueError("The opes status of this space has changed!! Check your files")
        elif opes_before is True and self._opes is False:
            raise ValueError("The opes status of this space has changed!! Check your files")

        return self

    def add_line(self, line: FreeEnergyLine):
        """
        function to add a free energy line to the landscape
        :param line: the fes to add
        :return: the fes for the landscape
        """
        cv = line.cvs[0]
        line._metadata = self._metadata
        line.temperature = self.temperature
        self.lines[cv] = line
        return self

    def add_surface(self, surface: FreeEnergySurface):
        """
        function to add a free energy surface to the landscape
        :param surface: the fes to add
        :return: the fes for the landscape
        """
        surface._metadata = self._metadata
        surface.temperature = self.temperature
        self.surfaces.append(surface)
        return self

    def get_long_hills(self, time_resolution: int = 6, height_power: float = 1):
        """
        Function to turn the _hills into long format, and allow for time binning and height power conversion. Works for
        both OPES and MetaD
        :param time_resolution: bin the _data into time bins with this number of decimal places. Useful for producing a
        smaller long format _hills
        dataframe.
        :param height_power: raise the height to the power of this so to see _hills easier. Useful when plotting, and
        you want to see the small
        _hills.
        :return:
        """

        if self._hills is None:
            raise ValueError("The space needs some hills/kernels data!")

        height_label = 'height^' + str(height_power) if height_power != 1 else 'height'
        long_hills = self._hills.rename(columns={'height': height_label})
        long_hills[height_label] = long_hills[height_label].pow(height_power)

        if not self._opes:
            long_hills = long_hills.melt(value_vars=self.cvs + [height_label], id_vars=['time', 'walker'])
        elif self._opes:
            long_hills = (long_hills
                          .melt(value_vars=self.cvs + [height_label] + ['logweight'], id_vars=['time', 'walker'])
                          )

        long_hills = (long_hills
                      .assign(time=lambda x: x['time'].round(time_resolution))
                      .groupby(['time', 'walker', 'variable'], group_keys=False)
                      .mean()
                      .reset_index()
                      )

        return long_hills

    def get_hills_figures(self, **kwargs) -> dict[int, go.Figure]:
        """
        Function to get a dictionary of plotly figure objects summarising the dynamics and _hills for each walker in a
        metadynamics simulation.
        :return:
        """

        if self._hills is None:
            raise ValueError("The space needs some hills data!")

        long_hills = self.get_long_hills(**kwargs).groupby('walker')
        figs = {}
        for name, df in long_hills:
            figure = px.line(df, x='time', y='value', facet_row='variable', labels={'time': 'Time [ns]'},
                             template=custom_dark_template)
            figure.update_traces(line=dict(width=1))
            figure.update_yaxes(title=None, matches=None)
            figure.for_each_annotation(lambda a: a.update(text=a.text.split("=")[1]))
            figs[name] = figure

        return figs

    def get_hills_average_across_walkers(self, time_resolution: int = 5, with_metadata: bool = False):
        """
        function to get the average hill height, averaged across the walkers.
        :return:
        """

        if self._hills is None:
            raise ValueError("The space needs some hills data!")

        av_hills = (self._hills
                    .assign(time=lambda x: x['time'].round(time_resolution))
                    .groupby(['time'])
                    .mean()
                    .reset_index()
                    )

        av_hills = av_hills[['time', 'height']]

        if with_metadata:
            av_hills['temperature'] = self.temperature
            for key, value in self._metadata.items():
                av_hills[key] = value

        return av_hills

    def get_average_hills_figure(self, **kwargs):
        """
        function to get figure of average _hills across walkers
        :return:
        """

        if self._hills is None:
            raise ValueError("The space needs some hills data!")

        log_y = True if not self._opes else False

        av_hills = self.get_hills_average_across_walkers(**kwargs)
        figure = px.line(av_hills, x='time', y='height', log_y=log_y, template=custom_dark_template,
                         labels={'time': 'Time [ns]', 'height': 'Energy [kJ/mol]'}
                         )
        figure.update_traces(line=dict(width=1))
        return figure

    def get_hills_max_across_walkers(self, time_resolution: int = 5, with_metadata: bool = False):
        """
        function to get the average hill height, averaged across the walkers.
        :return:
        """

        if self._hills is None:
            raise ValueError("The space needs some hills data!")

        max_hills = (self._hills
                     .assign(time=lambda x: x['time'].round(time_resolution))
                     .groupby(['time'])
                     .max()
                     .reset_index()
                     )

        max_hills = max_hills[['time', 'height']]

        if with_metadata:
            max_hills['temperature'] = self.temperature
            for key, value in self._metadata.items():
                max_hills[key] = value

        return max_hills

    def get_max_hills_figure(self, **kwargs):
        """
        function to get figure of average _hills across walkers
        :return:
        """

        if self._hills is None:
            raise ValueError("The space needs some hills data!")

        log_y = True if not self._opes else False

        max_hills = self.get_hills_max_across_walkers(**kwargs)
        figure = px.line(max_hills, x='time', y='height', log_y=log_y, template=custom_dark_template,
                         labels={'time': 'Time [ns]', 'height': 'Energy [kJ/mol]'})
        figure.update_traces(line=dict(width=1))

        return figure

    @staticmethod
    def _reweight_traj_data(data: pd.DataFrame, cv: str | list[str], bins: int | list[int | float] = 200,
                            temperature: float = 298):
        """
        Function to reweight a _data frame using weights. Can do both one dimensional binning and two-dimensional
        binning
        :param data: _data frame to reweight
        :param cv: the collective variable you are reweighting over
        :param bins: number of bins, or a list of bin boundaries
        :param temperature to get the population
        :return: reweighted dataframe
        """
        if type(cv) == str:
            histogram = np.histogram(a=data[cv], bins=bins, weights=data['weight'], density=True)
            x_points = [(histogram[1][i] + histogram[1][i + 1]) / 2 for i in range(0, len(histogram[1]) - 1)]
            if type(bins) == list:
                x_widths = [(histogram[1][i+1] - histogram[1][i]) for i in range(0, len(histogram[1]) - 1)]
                pop = [histogram[0][i] * x_widths[i] for i in range(0, len(histogram[0]))]
            else:
                pop = [p for p in histogram[0]]

            reweighted_data = pd.DataFrame({
                'population': pop,
                cv: x_points
            }).pipe(boltzmann_population_to_energy, temperature=temperature)

        elif type(cv) == list and len(cv) == 2:
            histogram = np.histogram2d(x=data[cv[0]], y=data[cv[1]], bins=bins, weights=data['weight'], density=True)
            x_points = [(histogram[1][i] + histogram[1][i + 1]) / 2 for i in range(0, len(histogram[1]) - 1)]
            y_points = [(histogram[2][i] + histogram[2][i + 1]) / 2 for i in range(0, len(histogram[2]) - 1)]
            reweighted_data = (pd.DataFrame(histogram[0], index=x_points, columns=y_points)
                               .melt(var_name=cv[1], value_name="population", ignore_index=False)
                               .reset_index(names=cv[0])
                               .pipe(boltzmann_population_to_energy, temperature=temperature)
                               )
        else:
            raise ValueError('Reweighting only supports one or two CVs at the moment')

        return reweighted_data

    def get_reweighted_surface(self, cvs: list[str, str], bins: list[int, int]):
        """
        Function to get a reweighted surface
        :param cvs: list with the two cvs. The first will go on the x-axis, the second on the y-axis
        :param bins: list with two integers for the number of bins in each CV
        :return: a free energy surface
        """
        data = []
        for w, t in self.trajectories.items():
            if cvs[0] in t.cvs and cvs[1] in t.cvs:
                data.append(t.get_data())
        if not data:
            raise ValueError("no trajectories in this space have that CV")
        data = pd.concat(data).sort_values('time')
        fes_data = self._reweight_traj_data(data, cvs, bins, self.temperature)
        surface = FreeEnergySurface(fes_data, temperature=self.temperature, metadata=self._metadata)
        return surface
    
    @staticmethod
    def _reweight_traj_list(traj_list: list, cv: str, bins: int | list[int | float] = 200, n_timestamps: int = None,
                            verbosity: bool = False, conditions: str | list[str] = None, temperature: float = 298,
                            ) -> (pd.DataFrame | dict[pd.DataFrame]):
        """
        Function to reweight a list of trajectories.
        :param traj_list: list of trajectories to reweight.
        :param cv: the cv in which to get the reweight.
        :param bins: number of bins, or a list of bin boundaries.
        :param n_timestamps: number of time stamps to have in the _time_data.
        :param verbosity: print progress?
        :param conditions: some query style conditions to put on the histogram.
        :param temperature: temperature to get the population.
        :return: reweighted trajectory data.
        """
        # grab the data from the trajectories
        data = []
        for t in traj_list:
            if cv in t.cvs:
                data.append(t.get_data())
            else:
                raise ValueError("no trajectories in this space have that CV")
        data = pd.concat(data).sort_values('time')

        # filter the data if there is a condition
        if conditions:
            if type(conditions) == str:
                data = data.query(conditions)
            elif type(conditions) == list:
                for c in conditions:
                    data = data.query(c)
        
        # reweight the data
        if n_timestamps is None:
            fes_data = (FreeEnergySpace
                        ._reweight_traj_data(data, cv, bins, temperature)
                        .filter([cv, 'energy', 'population'])
                        )
        elif type(n_timestamps) == int:
            fes_data = {}
            max_time = data['time'].max()
            for i in range(0, n_timestamps):
                time = (i + 1) * max_time / n_timestamps
                filtered_data = data.query('time <= @time')
                fes_data[i+1] = (FreeEnergySpace
                                 ._reweight_traj_data(filtered_data, cv, bins, temperature)
                                 .filter([cv, 'energy', 'population'])
                                 )
                if verbosity:
                    print(f"Made histogram for {i} timestamp")
        else:
            raise ValueError("n_timestamps needs to be None or integer!")

        return fes_data

    def get_reweighted_line(self, cv: str, bins: int | list[int | float] = 200, n_timestamps: int = None,
                            verbosity: bool = False, conditions: str | list[str] = None,
                            ) -> FreeEnergyLine:
        """
        Function to get a free energy line from a free energy space with meta trajectories in it, using weighted
        histogram analysis.
        :param cv: the cv in which to get the reweight
        :param bins: number of bins, or a list with the bin boundaries
        :param n_timestamps: number of time stamps to have in the _time_data
        :param verbosity: print progress?
        :param conditions: some query style conditions to put on the histogram
        :return:
        """
        # grab the trajectories and put them in a list
        traj_list = []
        for w, t in self.trajectories.items():
            traj_list.append(t)

        # reweight the trajectories
        fes_data = self._reweight_traj_list(traj_list, cv, bins, n_timestamps, verbosity, conditions, self.temperature)

        line = FreeEnergyLine(fes_data, temperature=self.temperature, metadata=self._metadata)
        return line
    
    @staticmethod
    def _calc_errors_from_walker_std(data_list: list[pd.DataFrame], cv: str, n_walkers: int,
                                     bins: int | list[int | float] = 200) -> pd.DataFrame:
        """
        Function to calculate the errors from the standard deviation of the multiple walkers.
        :param data_list: list of dataframes with data from each walker.
        :param cv: the collective variable to bin over.
        :param n_walkers: number of walkers.
        :param bins: number of bins, or a list of bin boundaries.
        :return: dataframe with errors.
        """
        if n_walkers != len(data_list):
            raise ValueError("n_walkers does not match the number of walkers in the data!")

        # bin the data
        data = (pd.concat(data_list)
                .assign(bin=lambda x: pd.cut(x[cv], bins))
                )

        # calculate errors
        binned_data = pd.DataFrame({
            cv: data.groupby('bin', observed=True).mean()[cv],
            'energy': data.groupby('bin', observed=True).mean()['energy'],
            'energy_err': data.groupby('bin', observed=True).std()['energy']/np.sqrt(n_walkers),
            'population': data.groupby('bin', observed=True).mean()['population'],
            'population_err': data.groupby('bin', observed=True).std()['population']/np.sqrt(n_walkers),
        }).dropna()

        return binned_data

    def get_reweighted_line_with_walker_error(self, cv: str, bins: int | list[int | float] = 200,
                                              n_timestamps: int = None, verbosity: bool = False,
                                              conditions: str | list[str] = None) -> FreeEnergyLine:
        """
        Function to get a free energy line from a free energy space with meta trajectories in it, using weighted
        histogram
        analysis. Errors are calculated from the standard deviation of the multiple walkers.
        :param cv: the cv in which to get the reweight
        :param bins: number of bins, or a list with the bin boundaries
        :param n_timestamps: number of time stamps to have in the _time_data
        :param verbosity: print progress?
        :param conditions: some query style conditions to put on the histogram
        :return:
        """
        if self.n_walker == 1:
            raise ValueError("there is only data from one walker in this space!")
        
        # reweight each trajectory individually
        fes_data = []
        for w, t in self.trajectories.items():
            if verbosity:
                print(f"Getting reweighted data for walker {w}")
            fes_data.append(self._reweight_traj_list([t], cv, bins, n_timestamps=n_timestamps,
                                                     verbosity=verbosity, conditions=conditions,
                                                     temperature=self.temperature)
                            )
        
        if type(fes_data[0]) == dict:
            binned_fes = {}
            for time_idx in range(0, n_timestamps):
                time_fes_data = [fes_data[walker_idx][time_idx+1] for walker_idx in range(0, len(fes_data))]
                binned_fes[time_idx+1] = self._calc_errors_from_walker_std(time_fes_data, cv, self.n_walker, bins,)
        else:
            binned_fes = self._calc_errors_from_walker_std(fes_data, cv, self.n_walker, bins,)

        line = FreeEnergyLine(binned_fes, temperature=self.temperature, metadata=self._metadata)

        return line

    def get_data(self, with_metadata: bool = False):
        """
        function to get the _data from a free energy shape
        :param with_metadata:
        :return:
        """
        data = self._hills.copy()
        if with_metadata:
            data['temperature'] = self.temperature

            if self._metadata:
                for key, value in self._metadata.items():
                    data[key] = value

        return data
