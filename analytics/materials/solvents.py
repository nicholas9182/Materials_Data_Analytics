from analytics.materials.material_lists import common_solvents


class Solvent():

    def __init__(self, name: str = None, pH: float = None) -> None:
        self._name = self._get_solvent_from_list(name)[0]
        self._formula = self._get_solvent_from_list(name)[1]
        self._pH = pH

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
    