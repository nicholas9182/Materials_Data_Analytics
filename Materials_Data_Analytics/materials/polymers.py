

class Polymer():
    """
    Master class for all polymers.  An instance of this class represents a repeating unit of this materials
    """
    def __init__(self, name: str = None) -> None:
        self._name = name.lower()

    @property
    def name(self):
        return self._name.upper()


class NType(Polymer):
    """
    Class for n-types
    """
    def __init__(self, name: str = None, formal_reduction_potential: float = None) -> None:
        super().__init__(name=name)
        self._formal_reduction_potential = formal_reduction_potential

    @property
    def formal_reduction_potential(self):
        return self._formal_reduction_potential
    
    @formal_reduction_potential.setter
    def formal_reduction_potential(self, value):
        self._formal_reduction_potential = value
    

class PType(Polymer):
    """
    class for p-types
    """
    def __init__(self, name: str = None, formal_oxidation_potential: float = None) -> None:
        super().__init__(name=name)
        self._formal_oxidation_potential = formal_oxidation_potential

    @property
    def formal_oxidation_potential(self):
        return self._formal_oxidation_potential
    
    @formal_oxidation_potential.setter
    def formal_oxidation_potential(self, value):
        self._formal_oxidation_potential = value

