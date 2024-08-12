from analytics.experiment_modelling.core import ElectrochemicalExperiment
from analytics.materials.electrolytes import Electrolyte
from analytics.materials.ions import Cation, Anion
import pandas as pd
import numpy as np
from typing import Union
import plotly.express as px


class CyclicVoltammogram(ElectrochemicalExperiment):

    def __init__(self,  
                 potential: Union[list, pd.Series, np.array] = None,
                 current: Union[list, pd.Series, np.array] = None,
                 cycle: Union[list, pd.Series, np.array] = None,
                 time: Union[list, pd.Series, np.array] = None,
                 electrolyte: Electrolyte = None,
                 metadata: dict = None,
                 keep_cycle_1: bool = False
                 ) -> None:
        
        super().__init__(electrolyte, metadata=metadata)

        self._keep_cycle_1 = keep_cycle_1

        self._data = pd.DataFrame()

        if len(potential) and len(current) and len(cycle) and len(time) != 0:
            self._data = pd.DataFrame({'potential': potential, 'current': current, 'cycle': cycle, 'time': time})
            self._wrangle_data()

    def _wrangle_data(self):
        """
        Function to clean the data
        """
        self._data = self._data.dropna()

        if self._keep_cycle_1 is False:
            self._data = self._data.query('cycle != 1')

        self._remake_cycle_numbers()

        if 'cycle' in self._data.columns:
            self._data['cycle'] = self._data['cycle'].astype(int).sort_values()

        if 'time' in self._data.columns:
            self._data['time'] = self._data['time'].astype(float).sort_values()

        self._data = self._data.groupby(['cycle'], group_keys=False).apply(lambda df: df.pipe(self._add_endpoints))

        return self
    
    @staticmethod
    def _add_endpoints(df):
        first_redox = df['redox'].iloc[0]
        last_redox = df['redox'].iloc[-1]
        last_index = df.query('redox == @first_redox').index[-1]
        first_chunk = df.query('redox == @first_redox')
        second_chunk = df.query('redox == @first_redox and index == @last_index').assign(redox = last_redox)
        third_chunk = df.query('redox == @last_redox') 
        return pd.concat([first_chunk, second_chunk, third_chunk])
    
    def _remake_cycle_numbers(self):
        """
        Function to remake the cycle numbers so that one cycle is from peak to peak
        """

        self._data = (self._data
                      .groupby(['cycle'], group_keys=False)
                      .apply(lambda df: df.assign(
                          _peak_potential = lambda x: x['potential'].max(),
                          _trough_potential = lambda x: x['potential'].min()
                          ))
                      )
        
        new_cycle_num_list = []
        cycle = -0.1
        for index, row in self._data.iterrows():
            if row['potential'] == row['_peak_potential'] or row['potential'] == row['_trough_potential']:
                cycle += 1
            new_cycle_num_list.append(cycle)

        self._data = (self._data
                      .assign(cycle = new_cycle_num_list)
                      .assign(cycle = lambda x: x['cycle'].div(2).round(0).astype(int))
                      .drop(columns=['_peak_potential', '_trough_potential'])
                      .groupby(['cycle'], group_keys=False)
                      .apply(lambda df: df.assign(redox = lambda x: ['reduction' if i < 0 else 'oxidation' for i in x['potential'].diff()]))
                      .query('index > index.min()+5')
                      )

        return self

    @property
    def data(self, with_metadata = True) -> pd.DataFrame:
        
        data = self._data

        if with_metadata:
            for k in self.metadata.keys():
                data = data.assign(k = self.metadata[k])

        return data

    @classmethod
    def from_benelogic(cls, path: str, electrolyte: Electrolyte = None, **kwargs):
        """
        Function to make a CyclicVoltammogram object from a Benelogic file
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
        data = self.data.assign(cycle_direction = lambda x: x['cycle'].astype('str') + ', ' + x['redox'])
        px.line(data, x='potential', y='current', color='cycle_direction', markers=True, 
                labels={'potential': 'Potential [V]', 'current': 'Current [A]'}, **kwargs).show()
        
        return self
    
    def show_current_time(self, **kwargs):
        """
        Function to plot the current vs time
        """
        data = self.data.assign(cycle_direction = lambda x: x['cycle'].astype('str') + ', ' + x['redox'])
        px.line(data, x='time', y='current', color='cycle_direction', markers=True, 
                labels={'time': 'Time [s]', 'current': 'Current [A]'}, **kwargs).show()
        
        return self
    
    def show_potential_time(self, **kwargs):
        """
        Function to plot the potential vs time
        """
        data = self.data.assign(cycle_direction = lambda x: x['cycle'].astype('str') + ', ' + x['redox'])
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
