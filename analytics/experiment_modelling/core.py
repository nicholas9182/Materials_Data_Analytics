from analytics.materials.electrolytes import Electrolyte

class Measurement():

    def __init__(self, metadata: dict = None) -> None:
        self._metadata = metadata if metadata is not None else {}

    @property
    def metadata(self) -> dict:
        return self._metadata
    
    @metadata.setter
    def metadata(self, value: dict):
        """
        Add items to the metadata dictionary.  If the key is already in the metadata, then it will overwrite the
        existing value
        """
        if type(value) != dict:
            raise ValueError('metadata must be a dictionary')
        
        for k in value.keys():
            self._metadata[k] = value[k]

    @staticmethod
    def get_root(x1, x2, y1, y2):
        """
        Function to get the root of the line given by two points
        """
        slope = (y2 - y1) / (x2 - x1)
        intercept = y1 - slope * x1
        root = -intercept / slope
        return root


class ElectrochemicalMeasurement(Measurement):

    def __init__(self, electrolyte: Electrolyte = None, metadata: dict = None) -> None:
        super().__init__(metadata=metadata)
        self._electrolyte = electrolyte

    @property
    def electrolyte(self) -> Electrolyte:
        return self._electrolyte


class ScatteringMeasurement(Measurement):

    def __init__(self, metadata: dict = None) -> None:
        super().__init__(metadata=metadata)

