from Materials_Data_Analytics.materials.material_lists import common_solvents
from Materials_Data_Analytics.materials.material_lists import solvent_phs
from Materials_Data_Analytics.materials.solutes import Solute


class Solvent():

    def __init__(self, name: str = None) -> None:

        self._name = self._get_solvent_from_list(name)[0] if name is not None else None
        self._formula = self._get_solvent_from_list(name)[1] if name is not None else None
        self._pH = solvent_phs[self._name] if self._name is not None else None 

    @classmethod
    def from_custom_inputs(cls, name: str = None, formula: str = None, pH: float = None):
        """
        Class method to create an instance of Solvent with custom inputs
        """
        solvent = cls()
        solvent._name = name
        solvent._formula = formula
        solvent._pH = pH
        return solvent

    @property
    def pH(self):
        return self._pH

    @property
    def name(self):
        return self._name.capitalize()
    
    @property
    def formula(self):
        return self._formula.upper()

    @staticmethod
    def _get_solvent_from_list(name: str):
        """
        Function to search for a solvent in the common_solvents dictionary, given that the name
        can be either the chemistry or the common name
        """
        name = name.lower()
        if name in common_solvents.keys():
            solvent_name = name
            solvent_formula = common_solvents[name]
        elif name in common_solvents.values():
            solvent_name = list(common_solvents.keys())[list(common_solvents.values()).index(name)]
            solvent_formula = name
        else:
            raise ValueError(f'Solvent not found in common_solvents list, {common_solvents.keys()}')

        return solvent_name, solvent_formula
    
    def __str__(self) -> str:
        return f'{self.name} solvent, {self.formula}'
    
    