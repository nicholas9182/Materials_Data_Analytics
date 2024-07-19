from analytics.materials.solvents import Solvent
from analytics.materials.ions import Ion, Cation, Anion
from analytics.materials.solutes import Solute


class Electrolyte():
    """
    Master class for an electrolyte
    """
    def __init__(self, solvent: Solvent, cation: Cation, anion: Anion, solutes: list[Solute] = None, pH: float = None, temperature: float = 298) -> None:
        self._pH = pH
        self._temperature = temperature        
        self._cation = cation
        self._anion = anion
        self._solvent = solvent
        self._solutes = solutes

    @property
    def solutes(self):
        return self._solutes

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

