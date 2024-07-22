from analytics.materials.solvents import Solvent
from analytics.materials.ions import Ion, Cation, Anion
from analytics.materials.solutes import Solute


class Electrolyte():
    """
    Master class for an electrolyte
    :param solvent: the solvent in the electrolyte
    :param cation: the cation in the electrolyte
    :param anion: the anion in the electrolyte
    :param concentrations: the concentrations of the ions and solutes in the electrolyte, given as a dictionary. concentrations in mol/dm^3
    :param solutes: the solutes in the electrolyte
    :return: an electrolyte object
    """
    def __init__(self, solvent: Solvent, cation: Cation | list[Cation], anion: Anion | list[Anion], concentrations: dict[Cation|Anion|Solute, float], 
                 solute: Solute | list[Solute] = None, pH: float = None, temperature: float = 298) -> None:
        self._pH = pH
        self._temperature = temperature        
        
        if type(cation) == Cation: 
            self._cation = [cation]
        elif type(cation) == list:
            for c in cation:
                if type(c) != Cation:
                    raise ValueError('Cation must be a Cation object')
            self._cation = cation
        else:
            raise ValueError('Cation must be a Cation object')
        
        if type(anion) == Anion:
            self._anion = [anion]
        elif type(anion) == list:
            for a in anion:
                if type(a) != Anion:
                    raise ValueError('Anion must be an Anion object')
            self._anion = anion
        else:
            raise ValueError('Anion must be an Anion object')
        
        if type(solvent) == Solvent:
            self._solvent = solvent
        else:
            raise ValueError('Solvent must be a Solvent object')

        if type(solute) == Solute or type(solute) == list[Solute]:
            self._solute = [solute]
        else:
            raise ValueError('Solute must be a Solute object')

        self._concentrations = concentrations
        self._check_concentration_entries(concentrations)

    def _check_concentration_entries(self, concentrations: dict[Cation|Anion|Solute, float]):
        """
        Function to reformat the concentrations dictionary so the keys are lowercase names
        """
        for c in concentrations.keys():
            if c not in self._cation and c not in self._anion and c not in self._solute:
                raise ValueError(f'{c} is not in the electrolyte')
            
    @property
    def concentrations(self):
        return self._concentrations

    @property
    def solute(self):
        if len(self._solute) == 1:
            return self._solute[0]
        else:
            return self._solute

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
        if len(self._cation) == 1:
            return self._cation[0]
        else:
            return self._cation
    
    @property
    def anion(self):
        if len(self._anion) == 1:
            return self._anion[0]
        else:
            return self._anion

