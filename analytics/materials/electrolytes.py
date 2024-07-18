from analytics.materials.solvents import Solvent
from analytics.materials.ions import Ion


class Electrolyte():
    """
    Master class for an electrolyte
    """
    def __init__(self, solvent: Solvent, cation: Ion, anion: Ion, pH: float = None, temperature: float = 298) -> None:
        self._pH = pH
        self._temperature = temperature

        if cation.charge <= 0:
            raise ValueError('Your cation needs to have a positive charge!')
        if anion.charge >= 0:
            raise ValueError('Your anion needs to have a negative charge!')
        
        self._cation = cation
        self._anion = anion
        self._solvent = solvent

    @property
    def pH(self):
        return self._pH
    
    @property
    def temperature(self):
        return self._temperature

    @property
    def solvent(self):
        return self._solvent
    
    @property
    def cation(self):
        return self._cation
    
    @property
    def anion(self):
        return self._anion

