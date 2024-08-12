from analytics.experiment_modelling.core import ElectrochemicalExperiment
from analytics.materials.electrolytes import Electrolyte
from analytics.materials.ions import Cation, Anion
import pandas as pd
import numpy as np
from typing import Union
import plotly.express as px
import scipy.integrate as integrate


class CyclicVoltammogram(ElectrochemicalExperiment):

    def __init__(self,  
                 potential: Union[list, pd.Series, np.array] = None,
                 current: Union[list, pd.Series, np.array] = None,
                 cycle: Union[list, pd.Series, np.array] = None,
                 time: Union[list, pd.Series, np.array] = None,
                 electrolyte: Electrolyte = None,
                 metadata: dict = None
                 ) -> None:
        
        super().__init__(electrolyte, metadata=metadata)

        self._data = pd.DataFrame()
        self._max_cycle = max(cycle)

        if len(potential) and len(current) and len(cycle) and len(time) != 0:
            self._data = (pd
                          .DataFrame({'potential': potential, 'current': current, 'cycle': cycle, 'time': time})
                          .pipe(self._wrangle_data)
                          )

    def _wrangle_data(self, data):
        """
        Function to clean the data
        """
        data = (data
                .dropna()
                .pipe(self._remake_cycle_numbers)
                .assign(cycle = lambda x: x['cycle'].astype(int))
                .assign(time = lambda x: x['time'].astype(float))
                .assign(potential = lambda x: x['potential'].astype(float))
                .assign(current = lambda x: x['current'].astype(float))
                .sort_values(by=['cycle', 'time'])
                .groupby(['cycle'], group_keys=False)
                .apply(lambda df: df.pipe(self._add_endpoints, df['cycle'].iloc[0]))
                .reset_index(drop=True)
                )
        
        return data
    
    def _add_endpoints(self, df, cycle):
        """
        Function to double up data points at the end of the cycles so that each cycle is complete when filtered by direction
        """
        if cycle != self._max_cycle:
            first_redox = df['direction'].iloc[0]
            last_redox = df['direction'].iloc[-1]
            last_index = df.query('direction == @first_redox').index[-1]
            first_chunk = df.query('direction == @first_redox')
            second_chunk = df.query('direction == @first_redox and index == @last_index').assign(direction = last_redox)
            third_chunk = df.query('direction == @last_redox') 
            return pd.concat([first_chunk, second_chunk, third_chunk])
        else:
            return df
    
    @staticmethod
    def _remake_cycle_numbers(data):
        """
        Function to remake the cycle numbers so that one cycle is from peak to peak
        """
        data = (data
                .groupby(['cycle'], group_keys=False)
                .apply(lambda df: df.assign(
                    _peak_potential = lambda x: x['potential'].max(),
                    _trough_potential = lambda x: x['potential'].min()
                    ))
                )
        
        new_cycle_num_list = []
        cycle = -0.1

        for index, row in data.iterrows():
            if row['potential'] == row['_peak_potential'] or row['potential'] == row['_trough_potential']:
                cycle += 1
            new_cycle_num_list.append(cycle)

        data = (data
                .assign(cycle = new_cycle_num_list)
                .assign(cycle = lambda x: x['cycle'].div(2).round(0).astype(int))
                .drop(columns=['_peak_potential', '_trough_potential'])
                .groupby(['cycle'], group_keys=False)
                .apply(lambda df: df.assign(direction = lambda x: ['reduction' if i < 0 else 'oxidation' for i in x['potential'].diff()]))
                .query('index > index.min()+5')
                )

        return data

    @property
    def data(self, with_metadata = True) -> pd.DataFrame:
        
        data = self._data

        if with_metadata:
            for k in self.metadata.keys():
                data = data.assign(k = self.metadata[k])

        return data

    @classmethod
    def from_biologic(cls, path: str, electrolyte: Electrolyte = None, **kwargs):
        """
        Function to make a CyclicVoltammogram object from a biologic file
        """
        data = pd.read_table(path, sep='\s+', names=['potential', 'current', 'cycle', 'time'], skiprows=1)
        cv = cls(electrolyte=electrolyte, potential=data['potential'], current=data['current'], cycle=data['cycle'], time=data['time'], **kwargs)
        return cv
    
    def drop_cycles(self, cycles: list[int]):
        """
        Function to remove cycles from the data
        """
        self._data = self._data.query('cycle not in @cycles')
        return self
    
    def show_current_potential(self, **kwargs):
        """
        Function to plot the cyclic voltammogram
        """
        data = self.data.assign(cycle_direction = lambda x: x['cycle'].astype('str') + ', ' + x['direction'])
        px.line(data, x='potential', y='current', color='cycle_direction', markers=True, 
                labels={'potential': 'Potential [V]', 'current': 'Current [A]'}, **kwargs).show()
        
        return self
    
    def show_current_time(self, **kwargs):
        """
        Function to plot the current vs time
        """
        data = self.data.assign(cycle_direction = lambda x: x['cycle'].astype('str') + ', ' + x['direction'])
        px.line(data, x='time', y='current', color='cycle_direction', markers=True, 
                labels={'time': 'Time [s]', 'current': 'Current [A]'}, **kwargs).show()
        
        return self
    
    def show_potential_time(self, **kwargs):
        """
        Function to plot the potential vs time
        """
        data = self.data.assign(cycle_direction = lambda x: x['cycle'].astype('str') + ', ' + x['direction'])
        px.line(data, x='time', y='potential', color='cycle_direction', markers=True, 
                labels={'time': 'Time [s]', 'potential': 'Potential [V]'}, **kwargs).show()
        
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
    
    def _integrate_curves(self, data: pd.DataFrame, direction: str, valence: str):
        """
        Function to add a point to the current and time arrays to calculate the integral
        """

        if valence == 'positive':
            current_data = data.query('current > 0')
        elif valence == 'negative':
            current_data = data.query('current < 0')
        else:
            raise ValueError('Valence must be either positive or negative')
        
        if direction == 'oxidation':
            int_current = current_data.query('direction == "oxidation"')['current'].to_numpy()
            int_time = current_data.query('direction == "oxidation"')['time'].to_numpy()
        elif direction == 'reduction':
            int_current = current_data.query('direction == "reduction"')['current'].to_numpy()
            int_time = current_data.query('direction == "reduction"')['time'].to_numpy()
        else:
            raise ValueError('direction must be either oxidation or reduction')
        
        if direction == 'oxidation' and valence == 'negative':
            max_index = data.query('direction == "oxidation" and current < 0').index.max() + 1
            y1 = int_current[-1]
            x1 = int_time[-1]
            y2 = data.query('index == @max_index')['current'].iloc[0]
            x2 = data.query('index == @max_index')['time'].iloc[0]
            root = self.get_root(x1, x2, y1, y2)
            int_current = np.append(int_current, 0)
            int_time = np.append(int_time, root)

        elif direction == 'oxidation' and valence == 'positive':
            min_index = data.query('direction == "oxidation" and current > 0').index.min() - 1
            y1 = int_current[0]
            x1 = int_time[0]
            y2 = data.query('index == @min_index')['current'].iloc[0]
            x2 = data.query('index == @min_index')['time'].iloc[0]
            root = self.get_root(x1, x2, y1, y2)
            int_current = np.append(0, int_current)
            int_time = np.append(root, int_time)

        elif direction == 'reduction' and valence == 'negative':
            min_index = data.query('direction == "reduction" and current < 0').index.min() - 1
            y1 = int_current[0]
            x1 = int_time[0]
            y2 = data.query('index == @min_index')['current'].iloc[0]
            x2 = data.query('index == @min_index')['time'].iloc[0]
            root = self.get_root(x1, x2, y1, y2)
            int_current = np.append(0, int_current)
            int_time = np.append(root, int_time)

        elif direction == 'reduction' and valence == 'positive':
            max_index = data.query('direction == "reduction" and current > 0').index.max() + 1
            y1 = int_current[-1]
            x1 = int_time[-1]
            y2 = data.query('index == @max_index')['current'].iloc[0]
            x2 = data.query('index == @max_index')['time'].iloc[0] 
            root = self.get_root(x1, x2, y1, y2)
            int_current = np.append(int_current, 0)
            int_time = np.append(int_time, root)

        integral = abs(integrate.simpson(int_current, int_time))

        return integral


    def get_charge_passed(self):
        """
        Function to get the integrals of the current
        """ 
        data = self._data.query('cycle != 0 and cycle != @self._max_cycle').copy()

        integrals = (data
                     .groupby(['cycle', 'direction'], group_keys=False)
                     .apply(lambda df: (df
                                        .assign(anodic_charge = lambda x: self._integrate_curves(x, direction=x['direction'].iloc[0], valence='positive'))
                                        .assign(cathodic_charge = lambda x: self._integrate_curves(x, direction=x['direction'].iloc[0], valence='negative'))
                                        ))
                     .drop(columns=['potential', 'time', 'current'])
                     .drop_duplicates()
                     .reset_index(drop=True)
                     )
        
        return integrals
    
    def show_charge_passed(self, **kwargs):
        """
        Function to plot the charge passed
        """
        passed_charge = (self
                         .get_charge_passed()
                         .melt(value_vars=['anodic_charge', 'cathodic_charge'], id_vars=['cycle', 'direction'], var_name='type', value_name='charge')
                         )

        px.bar(passed_charge, x='cycle', y='charge', color='type', barmode='group', facet_row='direction', **kwargs).show()
        return self

