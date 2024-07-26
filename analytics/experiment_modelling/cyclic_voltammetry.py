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
                 metadata: dict = None
                 ) -> None:
        
        super().__init__(electrolyte, metadata=metadata)

        if potential is not None and current is not None and cycle is not None and time is not None:
            if len(potential) != len(current) or len(potential) != len(cycle) or len(potential) != len(time):
                raise ValueError('The length of the potential, current, cycle and time arrays must be the same')
            self._data = pd.DataFrame({'potential': potential, 'current': current, 'cycle': cycle, 'time': time})
        else:
            self._data = pd.DataFrame()

    @property
    def data(self, with_metadata = True) -> pd.DataFrame:
        
        data = self._data

        if with_metadata:
            for k in self.metadata.keys():
                data = data.assign(k = self.metadata[k])

        return data

    @classmethod
    def from_benelogic(cls, path: str, electrolyte: Electrolyte = None):
        """
        Function to make a CyclicVoltammogram object from a Benelogic file
        """
        cv = cls(electrolyte=electrolyte)
        data = pd.read_table(path, sep='\s+', names=['potential', 'current', 'cycle', 'time'], skiprows=1)
        cv._data = data
        return cv
    
    def drop_cycles(self, cycles: list[int]):
        """
        Function to remove cycles from the data
        """
        self._data = self._data.query('cycle not in @cycles')
        return self
    
    def make_plot(self, **kwargs):
        """
        Function to plot the cyclic voltammogram
        """

        px.line(self.data, x='potential', y='current', color='cycle', markers=True, 
                labels={'potential': 'Potential [V]', 'current': 'Current [A]'}, **kwargs).show()
        
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
