from analytics.materials.electrolytes import Electrolyte
import pandas as pd
import scipy as sp


class Measurement():
    """
    A general measurement class
    Main contributors:
    Nicholas Siemons
    Contributors:
    """
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
    def get_root_linear_interpolation(x: pd.Series, y: pd.Series):
        """
        Function to get the root of the line given by two points
        """

        if len(x) != 2 or len(y) != 2:
            raise ValueError('The x and y series must have exactly two points to find the root with linear interpolation')
        
        x1 = x.iloc[0]
        x2 = x.iloc[1]
        y1 = y.iloc[0]
        y2 = y.iloc[1]
        
        slope = (y2 - y1) / (x2 - x1)
        intercept = y1 - slope * x1
        root = -intercept / slope

        return root
    
    @staticmethod
    def get_root_cube_spline(x: pd.Series, y: pd.Series):
        """
        Function to get the root of a data using cubic spline interpolation
        """
        spline = sp.interpolate.CubicSpline(x, y)
        roots = spline.roots()

        if len(roots) == 0 or len(roots) > 1:
            raise ValueError('There are no or multiple roots at a crossing point! Check the noise around current crossing points')
        
        return roots[0]
        

class ElectrochemicalMeasurement(Measurement):
    """
    A general electrochemical measurement class
    Main contributors:
    Nicholas Siemons
    Contributors:
    """
    def __init__(self, electrolyte: Electrolyte = None, metadata: dict = None) -> None:
        super().__init__(metadata=metadata)
        self._electrolyte = electrolyte

    @property
    def electrolyte(self) -> Electrolyte:
        return self._electrolyte


class ScatteringMeasurement(Measurement):

    def __init__(self, metadata: dict = None) -> None:
        super().__init__(metadata=metadata)

