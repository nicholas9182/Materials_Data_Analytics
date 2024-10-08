from Materials_Data_Analytics.materials.material_lists import common_cations, common_anions
from Materials_Data_Analytics.materials.material_lists import cation_charges, anion_charges


class Ion():
    """
    Master class for ions
    """
    def __init__(self, name: str = None) -> None:
        self._name = self._get_ion_from_list(name)[0] if name is not None else None
        self._formula = self._get_ion_from_list(name)[1] if name is not None else None

        if self._formula in cation_charges.keys():
            self._charge = cation_charges[self._formula]
        elif self._formula in anion_charges.keys():
            self._charge = anion_charges[self._formula]
        elif self._formula is None:
            self._charge = None
        else:
            raise ValueError(f'Charge not found for {self._formula}')
        
    @classmethod
    def from_custom_inputs(cls, name: str = None, formula: str = None, charge: int = None):
        """
        Class method to create an instance of Ion with custom inputs
        """
        ion = cls()
        ion._name = name
        ion._formula = formula
        ion._charge = charge
        return ion

    @property
    def formula(self):
        return self._formula.upper()

    @property
    def charge(self):
        return self._charge

    @property
    def name(self):
        return self._name.capitalize()

    @staticmethod
    def _get_ion_from_list(name: str):
        """
        Function to search for an ion in the common_cations or common_anions dictionary, given that the name
        can be either the chemistry or the common name
        """
        name = name.lower()
        if name in common_cations.keys():
            ion_name = name
            ion_formula = common_cations[name]
        elif name in common_cations.values():
            ion_name = list(common_cations.keys())[list(common_cations.values()).index(name)]
            ion_formula = name
        elif name in common_anions.keys():
            ion_name = name
            ion_formula = common_anions[name]
        elif name in common_anions.values():
            ion_name = list(common_anions.keys())[list(common_anions.values()).index(name)]
            ion_formula = name
        else:
            raise ValueError(f'Ion not found in common_cations or common_anions list, {common_cations.keys()}')
        
        return ion_name, ion_formula
    
    def __str__(self) -> str:
        return f'{self.name} ion, {self.formula}'
    

class Cation(Ion):
    """
    Class for cations  
    """
    def __init__(self, name: str = None) -> None:
        super().__init__(name=name)
        

class Anion(Ion):
    """
    Class for anions
    """
    def __init__(self, name: str = None) -> None:      
        super().__init__(name=name)
        
