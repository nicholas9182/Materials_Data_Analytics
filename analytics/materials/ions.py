

class Ion():
    """
    Master class for ions
    """
    def __init__(self, charge: int, name: str = None) -> None:
        self._charge = charge
        self._name = name

    @property
    def charge(self):
        return self._charge

    @property
    def name(self):
        return self._name

class Cation(Ion):
    """
    Class for cations  
    """
    def __init__(self, charge: int,  name: str = None) -> None:
        super().__init__(charge=charge, name=name)
        if self._charge < 0:
            raise ValueError('Cation charge must be positive')


class Anion(Ion):
    """
    Class for anions
    """
    def __init__(self, charge: int, name: str = None) -> None:
        super().__init__(charge=charge, name=name)
        if self._charge > 0:
            raise ValueError('Anion charge must be negative')

  