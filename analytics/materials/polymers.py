

class Polymer():
    """
    Master class for all polymers.  An instance of this class represents a single repeat unit of this materials
    """
    def __init__(self, name: str = None) -> None:
        self._name = name

    @property
    def name(self):
        return self._name


class NType(Polymer):
    """
    Class for n-types
    """
    def __init__(self, formal_reduction_potential: float = None) -> None:
        super().__init__()
        self._formal_reduction_potential = formal_reduction_potential

    @property
    def formal_reduction_potential(self):
        return self._formal_reduction_potential
    

class PType(Polymer):
    """
    class for p-types
    """
    def __init__(self, formal_oxidation_potential: float = None) -> None:
        super().__init__()
        self._formal_oxidation_potential = formal_oxidation_potential

