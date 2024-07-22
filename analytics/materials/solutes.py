from analytics.materials.material_lists import common_solutes
from analytics.materials.material_lists import solute_formal_reduction_potentials
from analytics.materials.material_lists import solute_pkas
from analytics.materials.material_lists import solute_standard_reduction_potentials


class Solute():

    def __init__(self, name: str = None) -> None:
        
        self._name = self._get_solute_from_list(name)[0] if name is not None else None
        self._formula = self._get_solute_from_list(name)[1] if name is not None else None

        if self._name is not None and self._formula in solute_formal_reduction_potentials.keys():
            self._formal_reduction_potentials = solute_formal_reduction_potentials[self._formula]
        else:
            self._formal_reduction_potentials = None

        if self._name is not None and self._formula in solute_standard_reduction_potentials.keys():
            self._standard_reduction_potentials = solute_standard_reduction_potentials[self._formula]
        else:
            self._standard_reduction_potentials = None

        if self._name is not None and self._formula in solute_pkas.keys():
            self._pka = solute_pkas[self._formula]
        else:
            self._pka = None

    @classmethod
    def from_custom_inputs(cls, name: str = None, formula: str = None, formal_reduction_potential: dict[str, float] = None, pka: float = None):
        """
        Class method to create an instance of Solute with custom inputs
        """
        solute = cls()
        solute._name = name
        solute._formula = formula
        solute._formal_reduction_potential = formal_reduction_potential
        solute._pka = pka
        return solute
    
    @property
    def standard_reduction_potentials(self) -> dict[str, float]:
        return self._standard_reduction_potentials
    
    @property
    def pka(self) -> float:
        return self._pka

    @property
    def formal_reduction_potentials(self) -> dict[str, float]:
        return self._formal_reduction_potentials
    
    @property
    def formula(self) -> str:
        return self._formula.upper()
    
    @property
    def name(self) -> str:
        return self._name.capitalize()

    @staticmethod
    def _get_solute_from_list(name: str) -> tuple:
        """
        Function to search for a solute in the common_solutes dictionary, given that the name
        can be either the chemistry or the common name
        """
        name = name.lower()
        if name in common_solutes.keys():
            solute_name = name
            solute_formula = common_solutes[name]
        elif name in common_solutes.values():
            solute_name = list(common_solutes.keys())[list(common_solutes.values()).index(name)]
            solute_formula = name
        else:
            raise ValueError(f'Solute not found in common_solutes list, {common_solutes.keys()}')

        return solute_name, solute_formula
    
    def __str__(self) -> str:
        return f'{self.name} solute, {self.formula}'

