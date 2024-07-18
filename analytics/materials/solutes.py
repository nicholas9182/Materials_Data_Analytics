

class Solute():

    def __init__(self, name: str = None) -> None:
        self._name = name

    @property
    def name(self):
        return self._name
        