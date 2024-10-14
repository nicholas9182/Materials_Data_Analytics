from Materials_Data_Analytics.materials.electrolytes import Electrolyte
import pandas as pd
import scipy as sp
import numpy as np


class Measurement():
    """
    A general measurement class
    Main contributors:
    Nicholas Siemons
    Contributors:
    """
    def __init__(self, metadata: dict = None) -> None:
        self._data = pd.DataFrame()
        self._metadata = metadata if metadata is not None else {}

    @property
    def metadata(self) -> dict:
        return self._metadata
    
    @property
    def data(self) -> pd.DataFrame:
        return self._data
    
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
    
    @staticmethod
    def find_local_peak(data, y_col: str, x_col: str, initial_guess: float, window = 0.01, polynomial_order = 4) -> pd.DataFrame:
        """
        Function to find a local peak in a data set. 
        1. Fit a polynomial to the data
        2. Find the derivative of the polynomial
        3. Find the roots of the derivative
        4. Find the y values of the roots
        5. Find the maximum y value
        6. Return the x and y values of the peak
        :param data: The data frame with the x and y values
        :param y_col: The column name of the y values
        :param x_col: The column name of the x values
        :param initial_guess: The initial guess of the x value of the peak
        :param window: The x window around the initial guess to search for the peak
        :param polynomial_order: The order of the polynomial to fit to the data
        :return: a data frame with the interpolated polynomial and the peak values
        """

        # check fitting with even polynomial order
        if polynomial_order % 2 != 0:
            raise ValueError('Polynomial order must be an even number')
        
        # check that the length of data is atleast twice the polynomial order
        if len(data) < polynomial_order * 2:
            raise ValueError('Dataframe must have at least twice the polynomial order number of points. Increase the window size or decrease polynomial order')

        # find the window around the initial guess
        data = data[(data[x_col] < initial_guess + window) & (data[x_col] > initial_guess - window)]
        x_min = data[x_col].min()
        x_max = data[x_col].max()
        x_range = (x_min, x_max)

        # get the coefficients and covariance matrix of the polynomial fit
        coefficients = np.polyfit(data[x_col], data[y_col], polynomial_order, cov=True)[0]

        # create a polymonial object and its derivative from those coefficients
        polynomial = np.poly1d(coefficients)
        differential = np.polyder(polynomial)

        # find the roots of the derivative
        roots = np.roots(differential)

        # get the roots which are in the range of the x values
        critical_points_in_range = [np.real(x) for x in roots if x_range[0] <= x <= x_range[1] and np.isreal(x)]

        # get the y values of the critical points and the boundaries
        x_values_to_check = critical_points_in_range + [x_range[0], x_range[1]]
        y_values = [polynomial(x) for x in x_values_to_check]

        # get the maximum y value
        max_y = max(y_values)
        max_x = x_values_to_check[y_values.index(max_y)]

        # check that the max isnt an edge case
        if max_x == x_range[0] or max_x == x_range[1]:
            raise ValueError('The peak is at the edge of the window. Increase the window size or decrease polynomial order')

        data['fit_' + y_col] = polynomial(data[x_col])
        data[y_col + '_peak'] = data['fit_' + y_col].max()

        data[y_col + '_peak'] = max_y
        data[x_col + '_peak'] = max_x

        return data
        

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

