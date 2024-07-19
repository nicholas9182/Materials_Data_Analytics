from analytics.materials.material_lists import common_solutes

class Solute():

    def __init__(self, name: str = None) -> None:
        self._name = self._get_solute_from_list(name)[0]
        self._formula = self._get_solute_from_list(name)[1]

    @property
    def name(self) -> str:
        return self._name.capitalize()

    @property
    def formula(self) -> str:
        return self._formula.upper()

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
        