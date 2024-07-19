from analytics.materials.material_lists import common_cations, common_anions


class Ion():
    """
    Master class for ions
    """
    def __init__(self, charge: int, name: str = None) -> None:
        self._charge = charge
        self._name = self._get_ion_from_list(name)[0]
        self._formula = self._get_ion_from_list(name)[1]

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
    

class Cation(Ion):
    """
    Class for cations  
    """
    def __init__(self, charge: int,  name: str = None) -> None:
        super().__init__(charge=charge, name=name)
        if charge <= 0:
            raise ValueError('Cation charge must be positive')


class Anion(Ion):
    """
    Class for anions
    """
    def __init__(self, charge: int, name: str = None) -> None:
        super().__init__(charge=charge, name=name)
        if charge >= 0:
            raise ValueError('Anion charge must be negative')

  