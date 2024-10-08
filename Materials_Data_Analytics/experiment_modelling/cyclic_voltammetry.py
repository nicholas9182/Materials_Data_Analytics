from Materials_Data_Analytics.experiment_modelling.core import ElectrochemicalMeasurement
from Materials_Data_Analytics.materials.electrolytes import Electrolyte
from Materials_Data_Analytics.materials.ions import Cation, Anion
import pandas as pd
import numpy as np
from typing import Union
import plotly.express as px
import plotly.graph_objects as go
import scipy.integrate as integrate
import base64
import io


class CyclicVoltammogram(ElectrochemicalMeasurement):
    """
    A general class for the analysis of cyclic voltammograms.
    Main contributors:
    Nicholas Siemons
    Contributors:
    """
    def __init__(self,  
                 potential: Union[list, pd.Series, np.array] = None,
                 current: Union[list, pd.Series, np.array] = None,
                 time: Union[list, pd.Series, np.array] = None,
                 electrolyte: Electrolyte = None,
                 metadata: dict = None
                 ) -> None:
        
        super().__init__(electrolyte, metadata=metadata)

        self._data = pd.DataFrame()

        if len(potential) and len(current) and len(time) != 0:
            self._data = (pd
                          .DataFrame({'potential': potential, 'current': current, 'time': time})
                          .pipe(self._wrangle_data)
                          )
        
        self._max_cycle = self._data['cycle'].max()
        self._max_segment = self._data['segment'].max()
        self._steps_per_cycle = (self._data
                                 .query('cycle != cycle.max() and cycle != cycle.min()')
                                 .astype({'cycle': 'int'})
                                 .groupby(['cycle']).count()['potential'].mean().round(0)
                                 )

    def _wrangle_data(self, data, first_index = 5) -> pd.DataFrame:
        """
        Function to wrangle the data
        :param data: pd.DataFrame with columns potential, current, cycle, time
        """
        data = (data
                .query('index > @first_index')
                .dropna()
                .reset_index(drop=True)
                .sort_values(by=['time'])
                .assign(time = lambda x: x['time'] - x['time'].min())
                .pipe(self._find_current_roots)
                .pipe(self._determine_direction)
                .pipe(self._make_segments)
                .pipe(self._add_endpoints)
                .pipe(self._make_cycles)
                .pipe(self._check_types)
                .sort_values(by=['time', 'segment'])
                .reset_index(drop=True)
                )
        
        return data
    
    def _find_current_roots(self, data: pd.DataFrame) -> pd.DataFrame:
        """
        Function to find the time and voltage points where the current passes through 0
        """
        raw_roots_indexes = (data
                             .copy()
                             .assign(sign_diff = lambda x: np.sign(x['current']).diff())
                             .query('sign_diff != 0 and sign_diff.notna()')
                             .index
                             )
        
        for i in raw_roots_indexes:
            data_interpolate = data.iloc[i-1:i+1]
            time_root = self.get_root_linear_interpolation(data_interpolate['time'], data_interpolate['current'])
            potential_root = self.get_root_linear_interpolation(data_interpolate['potential'], data_interpolate['current'])

            new_row = (data_interpolate
                       .query('index == index.min()')
                       .assign(potential = potential_root)
                       .assign(time = time_root)
                       .assign(current = 0)
                       )

            data = pd.concat([data, new_row], ignore_index=True) 

        return data.sort_values(by=['time']).reset_index(drop=True)
    
    def _determine_direction(self, data: pd.DataFrame) -> pd.DataFrame:
        """
        Function to determine the direction of the cycle (either oxidation or reduction)
        :param data: pd.DataFrame with potential and time
        """
        potential = data['potential']
        potential_shifted = potential.shift(1)
        dv = potential - potential_shifted
        directions = ['oxidation'] if dv[1] > 0 else ['reduction']

        for i in range(1, len(dv)):
            if dv[i] > 0:
                directions.append('oxidation')
            elif dv[i] < 0:
                directions.append('reduction')
            elif dv[i] == 0:
                directions.append(directions[i-1])
        
        for i in range(1, len(directions)-1):
            if directions[i] != directions[i-1] and directions[i] != directions[i+1]:
                directions[i] = directions[i-1]

        data = data.assign(direction = directions)
            
        return data
    
    def _add_endpoints(self, data: pd.DataFrame) -> pd.DataFrame:
        """
        Function to double up data points at the end of the cycles so that each cycle is complete when filtered by direction
        """
        for s in range(data['segment'].min()+1, data['segment'].max()+1):
            prev_seg_num = s - 1
            prev_time = data.query('segment == @prev_seg_num')['time'].iloc[-1]
            prev_current = data.query('segment == @prev_seg_num')['current'].iloc[-1]
            prev_potential = data.query('segment == @prev_seg_num')['potential'].iloc[-1]

            new_row = (data
                    .query('segment == @s')
                    .query('time == time.min()')
                    .assign(time = prev_time, current = prev_current, potential = prev_potential)
                    )
            
            data = pd.concat([data, new_row], ignore_index=True).sort_values(by=['time'])

        return data
    
    def _make_cycles(self, data: pd.DataFrame) -> pd.DataFrame:
        """
        Function to find the cycles of the cyclic voltammogram
        """
        cycle = 0
        new_data = []
        for group_name, group_df in data.groupby('segment'):
            if group_name % 2 == 1:
                cycle += 1
            group_df['cycle'] = cycle
            new_data.append(group_df)

        data = pd.concat(new_data, ignore_index=True)

        return data
    
    def _check_types(self, data) -> pd.DataFrame:
        """
        Function to check the data types in the data columns
        :param data: pd.DataFrame with columns potential, current, cycle, time
        """
        return (data
                .assign(cycle = lambda x: x['cycle'].astype(int))
                .assign(time = lambda x: x['time'].astype(float))
                .assign(potential = lambda x: x['potential'].astype(float))
                .assign(current = lambda x: x['current'].astype(float))
                .assign(segment = lambda x: x['segment'].astype(int))
                .assign(direction = lambda x: x['direction'].astype(str))
                )
    
    @staticmethod
    def _make_segments(data) -> pd.DataFrame:
        """
        Function to find the segments of the cyclic voltammogram
        """
        direction_change = data['direction'] != data['direction'].shift(1)
        segments = direction_change.cumsum() - 1
        data = data.assign(segment=segments)
        return data

    @property
    def data(self) -> pd.DataFrame:
        
        data = self._data.copy()
        metadata = self.metadata

        for k in metadata.keys():
            data[k] = self.metadata[k]

        return data

    @classmethod
    def from_html_base64(cls, file_contents, source, scan_rate = None, **kwargs):
        """
        Function to make a CyclicVoltammogram object from an html file
        """
        content_type, content_string = file_contents.split(',')
        decoded = base64.b64decode(content_string)
        file_data = io.StringIO(decoded.decode('utf-8'))

        if source == 'biologic':
            data = pd.read_table(file_data, sep='\t')
            cv = cls.from_biologic(data=data, **kwargs)
        elif source == 'aftermath':
            data = pd.read_table(file_data, sep=",")
            cv = cls.from_aftermath(data=data, scan_rate=scan_rate, **kwargs)
        else:
            raise ValueError('The source must be either biologic or aftermath')
        
        return cv

    @classmethod
    def from_biologic(cls, path: str = None, data: pd.DataFrame = None, **kwargs):
        """
        Function to make a CyclicVoltammogram object from a biologic file
        """

        if path is None and data is not None:
            data = data
        elif path is not None and data is None:
            data = pd.read_table(path, sep='\t')

        data = (data
                .rename({'Ewe/V': 'potential', '<I>/mA': 'current', 'time/s': 'time'}, axis=1)
                .filter(['potential', 'current', 'time'])
                )

        cv = cls(potential=data['potential'], current=data['current'], time=data['time'], **kwargs)
        
        return cv
    
    @classmethod
    def from_aftermath(cls, path: str = None, scan_rate: float = None, data: pd.DataFrame = None, **kwargs):
        """
        Function to make a CyclicVoltammogram object from an AfterMath file
        """

        if path is None and data is not None:
            data = data
        elif path is not None and data is None:
            data = pd.read_table(path, sep=",")

        if type(scan_rate) != float:
            scan_rate = float(scan_rate)

        data = (data
                .rename({'Potential (V)': 'potential', 'Current (A)': 'current'}, axis=1)
                .filter(['potential', 'current'])
                .assign(current = lambda x: x['current']/1000)
                )
        
        dv = (pd
              .DataFrame({'dv': data['potential'] - data['potential'].shift(1)})
              .query('index > 0')
              .abs()
              .mean()
              .iloc[0]
              )
        
        time = [(i*dv)/(scan_rate/1000) for i in range(0, len(data))]
        
        cv = cls(potential=data['potential'], current=data['current'], time=time, **kwargs)

        return cv
    
    def drop_cycles(self, drop: list[int] | int = None, keep: list[int] | int = None) -> pd.DataFrame:
        """
        Function to edit which cycles are being considered
        """
        if type(drop) == int:
            drop = [int]

        if type(keep) == int:
            keep = [keep]

        if drop is not None:
            self._data = self._data.query('cycle not in @drop')

        if keep is not None:
            self._data = self._data.query('cycle in @keep')

        return self
    
    def get_current_potential_plot(self, **kwargs):
        """
        Function to plot the cyclic voltammogram
        """
        data = self.data.assign(cycle_direction = lambda x: x['cycle'].astype('str') + ', ' + x['direction'])

        figure = px.line(data, x='potential', y='current', color='cycle_direction', markers=True, 
                         labels={'potential': 'Potential [V]', 'current': 'Current [mA]'}, **kwargs)
        
        return figure
    
    def show_current_potential(self, **kwargs):
        """
        Function to show the cyclic voltammogram
        """
        figure = self.get_current_potential_plot(**kwargs)
        figure.show()
        return self
    
    def get_current_time_plot(self, **kwargs):
        """
        Function to plot the current vs time
        """
        data = self.data.assign(cycle_direction = lambda x: x['cycle'].astype('str') + ', ' + x['direction'])

        figure = px.line(data, x='time', y='current', color='cycle_direction', markers=True, 
                         labels={'time': 'Time [s]', 'current': 'Current [mA]', 'cycle_direction': 'Cycle, Direction'}, **kwargs)
        
        return figure
    
    def show_current_time(self, **kwargs):
        """
        Function to show the current vs time plot
        """
        figure = self.get_current_time_plot(**kwargs)
        figure.show()
        return self
    
    def get_potential_time_plot(self, **kwargs):
        """
        Function to plot the potential vs time
        """
        data = self.data.assign(cycle_direction = lambda x: x['cycle'].astype('str') + ', ' + x['direction'])
        
        figure = px.line(data, x='time', y='potential', color='cycle_direction', markers=True, 
                         labels={'time': 'Time [s]', 'potential': 'Potential [V]'}, **kwargs)
        
        return figure
    
    def show_potential_time(self, **kwargs):
        """
        Function to show the potential vs time plot
        """
        figure = self.get_potential_time_plot(**kwargs)
        figure.show()
        return self

    @property
    def pH(self) -> float:
        return self.electrolyte.pH
    
    @property
    def temperature(self) -> Electrolyte:
        return self.electrolyte.temperature
    
    @property
    def cation(self) -> Cation:
        return self.electrolyte.cation
    
    @property
    def anion(self) -> Anion:
        return self.electrolyte.anion
    
    @property
    def max_cycle(self) -> int:
        return self._max_cycle
    
    @property
    def steps_per_cycle(self) -> int:
        return self._steps_per_cycle
    
    #TODO: I dont like this function, it should be more general
    def _integrate_curves(self, data: pd.DataFrame, direction: str, valence: str) -> float:
        """
        Function to integrate data
        """
        if valence == 'positive':
            current_data = data.query('current >= 0')
        elif valence == 'negative':
            current_data = data.query('current <= 0')
        else:
            raise ValueError('Valence must be either positive or negative')

        int_current = current_data.query('direction == @direction')['current'].to_numpy()
        int_time = current_data.query('direction == @direction')['time'].to_numpy()
        integral = abs(integrate.simpson(int_current, x=int_time))

        return integral

    def get_charge_passed(self, average_segments = False) -> pd.DataFrame:
        """
        Function to get the integrals of the current
        """ 
        data = self._data.query('segment != 0 and segment != @self._max_segment').copy()

        integrals = (data
                     .groupby(['segment'], group_keys=False)
                     .apply(lambda df: (df
                                        .assign(anodic_charge = lambda x: self._integrate_curves(x, direction=x['direction'].iloc[0], valence='positive'))
                                        .assign(cathodic_charge = lambda x: self._integrate_curves(x, direction=x['direction'].iloc[0], valence='negative'))
                                        ))
                     .drop(columns=['potential', 'time', 'current'])
                     .drop_duplicates()
                     .reset_index(drop=True)
                     .assign(valence = lambda x: [-1 if i == 'reduction' else 1 for i in x.direction])
                     .assign(total_charge = lambda x: (x.anodic_charge - x.cathodic_charge)*x.valence)
                     .drop(columns=['valence'])
                     )
        
        if average_segments is True:
            integrals = (integrals
                         .groupby(['direction'], group_keys = False)
                         .apply(lambda df: (df
                                            .assign(anodic_charge_err = lambda x: x['anodic_charge'].std()/np.sqrt(x['anodic_charge'].count()))
                                            .assign(cathodic_charge_err = lambda x: x['cathodic_charge'].std()/np.sqrt(x['cathodic_charge'].count()))
                                            .assign(anodic_charge = lambda x: x['anodic_charge'].mean())
                                            .assign(cathodic_charge = lambda x: x['cathodic_charge'].mean())
                                            ))
                         .filter(['direction', 'anodic_charge', 'anodic_charge_err', 'cathodic_charge', 'cathodic_charge_err'])
                         .drop_duplicates()
                         )
        
        return integrals
    
    def get_maximum_charges_passed(self, average_sections = False) -> pd.DataFrame:
        """
        Function to get the maximum charges passed in each direction
        """ 
        data = (self
                ._data.query('segment != 0 and segment != @self._max_segment')
                .copy()
                .assign(section = lambda x: (x['current'] == 0).cumsum())
                )
        
        additional_points = (data
                             .assign(sos = lambda x: (x['current'] == 0))
                             .query('sos == True')
                             .assign(section = lambda x: x['section'] - 1)  
                             .drop(columns=['sos'])
                             )

        max_charges_passed = (pd
                              .concat([data, additional_points], axis=0)
                              .sort_values(by=['time'])
                              .reset_index(drop=True)
                              .groupby(['section'], group_keys=False)
                              .apply(lambda df: (df
                                                  .assign(total_charge = integrate.simpson(df.current, x=df.time))
                                                  .assign(t_min = df.time.min())
                                                  .assign(t_max = df.time.max())
                                                  .filter(['total_charge', 'section', 't_min', 't_max'])
                                                  .drop_duplicates()
                                                  ))
                              .assign(type = lambda x: ['anodic_charge' if i > 0 else 'cathodic_charge' for i in x['total_charge']])
                              .assign(total_charge = lambda x: x['total_charge'].abs())
                              .query('section != section.max() and section != section.min()')
                              )
        
        if average_sections is True:
            max_charges_passed = (max_charges_passed
                                  .groupby(['type'], group_keys=False)
                                  .apply(lambda df: (df
                                                      .assign(total_charge_err = lambda x: x.total_charge.std() / x.total_charge.count())
                                                      .assign(total_charge = lambda x: x.total_charge.mean())
                                                      .filter(['type', 'total_charge', 'total_charge_err'])
                                                      .drop_duplicates()
                                                      ))
                                  )
        
        return max_charges_passed
    
    def get_maximum_charge_passed_plot(self, **kwargs):
        """
        Function to plot the maximum charges passed
        """
        max_charges_passed = self.get_maximum_charges_passed()

        figure = px.line(max_charges_passed, x='section', y='total_charge', color='type', markers=True, labels={'total_charge': 'Charge [mC]', 'type': ''} , custom_data=['type'], **kwargs)
        
        return figure
    
    def get_charge_passed_plot(self, **kwargs):
        """
        Function to plot the charge passed
        """
        passed_charge = (self
                         .get_charge_passed()
                         .melt(value_vars=['anodic_charge', 'cathodic_charge', 'total_charge'], id_vars=['cycle', 'direction'], var_name='type', value_name='charge')
                         )

        figure = px.line(passed_charge, x='cycle', y='charge', color='type', facet_row='direction', markers=True,
                        labels={'charge': 'Charge [mC]', 'type': ''} , custom_data=['direction', 'type'], facet_row_spacing=0.12, **kwargs)
        figure.update_xaxes(dtick=1)

        return figure
    
    def show_charge_passed(self, **kwargs):
        """
        Function to show the charge passed
        """
        figure = self.get_charge_passed_plot(**kwargs)
        figure.show()
        return self
    
    def get_maximum_charge_integration_plot(self, section: int, **kwargs):
        """
        Function to return a plot showing the area integrated to get the maximum charges passed 
        """
        data = self.data.assign(cycle_direction = lambda x: x['cycle'].astype('str') + ', ' + x['direction'])
        charge_data = self.get_maximum_charges_passed()
        charge_valence = charge_data.query('section == @section')['type'].iloc[0]

        t_min = charge_data.query('section == @section')['t_min'].iloc[0]
        t_max = charge_data.query('section == @section')['t_max'].iloc[0]
        c_min = data.query('time >= @t_min and time <= @t_max')['current'].min()
        c_max = data.query('time >= @t_min and time <= @t_max')['current'].max()

        data_area = data.query('time >= @t_min and time <= @t_max')

        t_plot_min = t_min - (t_max - t_min) * 0.1
        t_plot_max = t_max + (t_max - t_min) * 0.1
        c_plot_min = c_min - (c_max - c_min) * 0.1
        c_plot_max = c_max + (c_max - c_min) * 0.1

        figure = px.line(data, x='time', y='current', color='cycle_direction', labels={'time': 'Time [s]', 'current': 'Current [mA]', 'cycle_direction': 'Cycle, Direction'}, **kwargs)

        if charge_valence == 'anodic_charge':
            figure.add_trace(go.Scatter(x=data_area['time'], y=data_area['current'], mode='lines', name='anodic charge', fill='tozeroy', fillcolor='#ADD8E6', line=dict(color='rgba(0, 0, 0, 0)')))
        if charge_valence == 'cathodic_charge':
            figure.add_trace(go.Scatter(x=data_area['time'], y=data_area['current'], mode='lines', name='cathodic charge', fill='tozeroy', fillcolor='#FFA07A', line=dict(color='rgba(0, 0, 0, 0)')))

        figure.add_shape(x0=t_plot_min, x1=t_plot_max, y0=0, y1=0, line=dict(color='black', width=2))
        figure.update_xaxes(range=[t_plot_min, t_plot_max])
        figure.update_yaxes(range=[c_plot_min, c_plot_max])

        return figure

    def get_charge_integration_plot(self, cycle: int, direction: str, **kwargs):
        """
        Function to return a plot which shows the area under a curve which is used for calculating the charges passed in the CV
        """
        if direction.lower() != "oxidation" and direction.lower() != "reduction":
            raise ValueError("Direction must be either oxidation or reduction")

        direction = direction.lower()

        data = self.data.copy().assign(cycle_direction = lambda x: x['cycle'].astype('str') + ', ' + x['direction'])
        segment = data.query('cycle == @cycle and direction == @direction')['segment'].values[0]
        data_area = data.query('segment == @segment')
        data_area_positive = data_area.query('current >= 0')
        data_area_negative = data_area.query('current <= 0')

        t_min = data_area['time'].min()
        t_max = data_area['time'].max()
        c_min = data_area['current'].min()
        c_max = data_area['current'].max()
        t_min = t_min - (t_max - t_min) * 0.1
        t_max = t_max + (t_max - t_min) * 0.1
        c_min = c_min - (c_max - c_min) * 0.1
        c_max = c_max + (c_max - c_min) * 0.1

        figure = px.line(data, x='time', y='current', color='cycle_direction', 
                         labels={'time': 'Time [s]', 'current': 'Current [mA]', 'cycle_direction': 'Cycle, Direction'}, **kwargs)
        
        figure.add_trace(go.Scatter(x=data_area_positive['time'], y=data_area_positive['current'], mode='lines', name='anodic charge', 
                                    fill='tozeroy', fillcolor='#ADD8E6', line=dict(color='rgba(0, 0, 0, 0)')))
        figure.add_trace(go.Scatter(x=data_area_negative['time'], y=data_area_negative['current'], mode='lines', name='cathodic charge', 
                                    fill='tozeroy', fillcolor='#FFA07A', line=dict(color='rgba(0, 0, 0, 0)')))
        
        figure.add_shape(x0=t_min, x1=t_max, y0=0, y1=0, line=dict(color='black', width=2))
        figure.update_xaxes(range=[t_min, t_max])
        figure.update_yaxes(range=[c_min, c_max])

        return figure

    def downsample(self, n: int | list[float] = 400) -> pd.DataFrame:
        """
        Function to downsample the data
        """
        # from a range and number of intervals, get the bin edges
        def get_bins(df, n):
            t_min = df['time'].min()
            t_max = df['time'].max()
            t_mids = np.linspace(t_min, t_max, n)
            dt = (t_max - t_min) / n
            bins = np.concatenate(([t_min - dt / 2], t_mids + dt / 2))
            return bins

        # get the original data and remove the additional segment points and current roots
        data = (self
                ._data
                .query('current != 0')
                .groupby(['segment'], group_keys=False)
                .apply(lambda df: df.query('index != index.max()'))
                )

        # downsample the data. Attention to use the .agg as it is much faster than .apply
        down_sampled_data = (data
                             .groupby(['cycle', 'direction', 'segment'], group_keys=False)
                             .apply(lambda df: df.assign(time_bin = pd.cut(df['time'], bins=get_bins(df, n), labels=False)))
                             .groupby(['cycle', 'direction', 'segment', 'time_bin'], group_keys=False)
                             .agg(
                                 potential = ('potential', 'mean'),
                                 current = ('current', 'mean'),
                                 time = ('time', 'mean')
                                 )
                             .reset_index()
                             .drop(columns=['time_bin'])
                             .pipe(self._add_endpoints)
                             .pipe(self._find_current_roots)
                             )
        
        self._data = down_sampled_data
        return self
    
