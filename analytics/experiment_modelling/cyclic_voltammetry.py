from analytics.experiment_modelling.core import ElectrochemicalMeasurement
from analytics.materials.electrolytes import Electrolyte
from analytics.materials.ions import Cation, Anion
import pandas as pd
import numpy as np
from typing import Union
import plotly.express as px
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

    def _wrangle_data(self, data) -> pd.DataFrame:
        """
        Function to wrangle the data
        :param data: pd.DataFrame with columns potential, current, cycle, time
        """
        data = (data
                .query('index > 5')
                .dropna()
                .reset_index(drop=True)
                .sort_values(by=['time'])
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
            new_row = pd.DataFrame({'potential': [potential_root], 'time': [time_root], 'current': [0]})     
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
        for s in range(1, data['segment'].max()+1):
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
        segments = [0]
        for index, row in data.iterrows():
            if index == 0:
                continue
            elif row['direction'] != data['direction'].iloc[index-1]:
                segments.append(segments[-1] + 1)
            elif row['direction'] == data['direction'].iloc[index-1]:
                segments.append(segments[-1])

        data = data.assign(segment = segments)

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

        figure = px.line(data, x='potential', y='current', color='cycle_direction', markers=True, title='Current vs Potential',
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

        figure = px.line(data, x='time', y='current', color='cycle_direction', markers=True, title='Current vs Time',
                         labels={'time': 'Time [s]', 'current': 'Current [mA]'}, **kwargs)
        
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
        
        figure = px.line(data, x='time', y='potential', color='cycle_direction', markers=True, title='Potential vs Time',
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
    
    def _integrate_curves(self, data: pd.DataFrame, direction: str, valence: str) -> float:
        """
        Function to add a point to the current and time arrays to calculate the integral
        """

        if valence == 'positive':
            current_data = data.query('current >= 0')
        elif valence == 'negative':
            current_data = data.query('current <= 0')
        else:
            raise ValueError('Valence must be either positive or negative')
        
        int_current = current_data.query('direction == @direction')['current'].to_numpy()
        int_time = current_data.query('direction == @direction')['time'].to_numpy()

        integral = abs(integrate.simpson(int_current, int_time))

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
    
    def get_charge_passed_plot(self, **kwargs):
        """
        Function to plot the charge passed
        """
        passed_charge = (self
                         .get_charge_passed()
                         .melt(value_vars=['anodic_charge', 'cathodic_charge'], id_vars=['cycle', 'direction'], var_name='type', value_name='charge')
                         )

        figure = px.bar(passed_charge, x='cycle', y='charge', color='type', barmode='group', facet_row='direction', title='Charge Passed for Each Cycle',
                        labels={'charge': 'Charge [mC]'} ,**kwargs)
        
        return figure
    
    def show_charge_passed(self, **kwargs):
        """
        Function to show the charge passed
        """
        figure = self.get_charge_passed_plot(**kwargs)
        figure.show()
        return self

