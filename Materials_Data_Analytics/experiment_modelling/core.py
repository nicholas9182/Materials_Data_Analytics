from Materials_Data_Analytics.materials.electrolytes import Electrolyte
import pandas as pd
import scipy as sp
import numpy as np
import plotly.graph_objects as go
import plotly.express as px


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

        try:
            # get the coefficients and covariance matrix of the polynomial fit
            coefficients = np.polyfit(data[x_col], data[y_col], polynomial_order, cov=True)[0]
        except np.linalg.LinAlgError:
            raise ValueError('The polynomial fit failed. Increase the window size or decrease polynomial order')

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
    
    def plot_pixel_map_px(self,
                       data: pd.DataFrame,
                       x: str,
                       y: str,
                       z: str,
                       colorscale: str, 
                       x_label: str = None,
                       y_label: str = None,
                       z_label: str = None,
                       xlim: tuple = None,
                       ylim: tuple = None,
                       log_scale: bool = False,
                       z_lower_cuttoff: float = None,
                       aspect: str = 'equal',
                       **kwargs) -> go.Figure:
        
        """
        Function for plotting a pixel map of data
        :param data: The data frame with the x, y, and z values
        :param x_col: The column name of the x values
        :param y_col: The column name of the y values
        :param z_col: The column name of the z values
        :param colorscale: The colorscale of the pixel map
        :param x_label: The x axis label
        :param y_label: The y axis label
        :param z_label: The z axis label
        :param xlim: The x limits of the plot
        :param ylim: The y limits of the plot
        :param width: The width of the plot
        :param height: The height of the plot
        :param title: The title of the plot
        :param log_scale: Whether to plot the z values on a log scale
        :param template: The plotly template to use
        :return: a plotly figure of the pixel map
        """
        if x not in data.columns or y not in data.columns or z not in data.columns:
            raise ValueError('The x, y, and z columns must be in the data frame')
        
        if z_lower_cuttoff is None and log_scale and min(data[z]) <= 0:
            raise ValueError('The z values must be positive to plot on a log scale. Either add a z_lower_cuttoff or remove the log_scale')
        
        data = data.query(f'{z} > 0')

        if xlim is not None:
            data = data.query(f'{x} >= {xlim[0]} and {x} <= {xlim[1]}')
        if ylim is not None:
            data = data.query(f'{y} >= {ylim[0]} and {y} <= {ylim[1]}')
        if z_lower_cuttoff is not None:
            data = data.query(f'{z} >= {z_lower_cuttoff}')

        square_data = data.pivot(index=y, columns=x, values=z)
        square_data = square_data.fillna(z_lower_cuttoff) if z_lower_cuttoff is not None else square_data.fillna(0)
        square_data = np.log10(square_data) if log_scale else square_data
        z_label = 'log(' + z_label + ')' if log_scale else z_label
        
        z_data = square_data.values
        x_vec = square_data.columns
        y_vec = square_data.index

        figure = px.imshow(z_data, x=x_vec, y=y_vec, color_continuous_scale=colorscale, aspect=aspect, **kwargs)
    
        if x_label is not None:
            figure.update_xaxes(title_text=x_label)
        if y_label is not None:
            figure.update_yaxes(title_text=y_label)
        if z_label is not None:
            figure.update_layout(coloraxis_colorbar={'title': z_label})

        return figure
    
    def plot_pixel_map_hv(self,
                       data: pd.DataFrame,
                       x: str,
                       y: str,
                       z: str,
                       log_scale: bool = False,
                       aspect: str = 'equal',
                       **kwargs):
        
        """
        Function for plotting a pixel map of data
        ùparam data: The data frame with the x, y, and z values
        :param x_col: The column name of the x values
        :param y_col: The column name of the y values
        :param z_col: The column name of the z values
        :param log_scale: Whether to plot the z values on a log scale
        :param aspect: The aspect ratio of the plot
        :return: a holoviews image of the pixel map
        """

        import holoviews as hv
        hv.extension('bokeh')        
        
        if x not in data.columns or y not in data.columns or z not in data.columns:
            raise ValueError('The x, y, and z columns must be in the data frame')
     

        square_data = data.pivot(index=y, columns=x, values=z)
        
        z_data = square_data.values
        x_vec = square_data.columns
        y_vec = square_data.index
        figure = hv.Image((x_vec, y_vec, z_data), kdims=['x', 'y'], vdims=['z']).opts(aspect = aspect,
                                                                                      logz = log_scale,
                                                                                      **kwargs)

        return figure
        
    
    def plot_contour_map(self, 
                         data: pd.DataFrame,
                         x: str,
                         y: str,
                         z: str,
                         colorscale: str,
                         x_label: str = None,
                         y_label: str = None,
                         z_label: str = None,
                         xlim: tuple = None,
                         ylim: tuple = None,
                         width: int = None,
                         height: int = None,
                         title = None,
                         ncontours: int = 200,
                         log_scale: bool = False, 
                         z_lower_cuttoff: float = None,
                         template: str = None) -> go.Figure:
        """
        Function for plotting a contour map of data
        :param data: The data frame with the x, y, and z values
        :param x_col: The column name of the x values
        :param y_col: The column name of the y values
        :param z_col: The column name of the z values
        :param colorscale: The colorscale of the contour map
        :param x_label: The x axis label
        :param y_label: The y axis label
        :param z_label: The z axis label
        :param xlim: The x limits of the plot
        :param ylim: The y limits of the plot
        :param width: The width of the plot
        :param height: The height of the plot
        :param title: The title of the plot
        :param ncontours: The number of contours to plot
        :param template: The plotly template to use
        :param log_scale: Whether to plot the z values on a log scale
        :return: a plotly figure of the contour map
        """
        if x not in data.columns or y not in data.columns or z not in data.columns:
            raise ValueError('The x, y, and z columns must be in the data frame')
        
        data = data.query(f'{z} > 0')

        if xlim is not None:
            data = data.query(f'{x} >= {xlim[0]} and {x} <= {xlim[1]}')
        if ylim is not None:
            data = data.query(f'{y} >= {ylim[0]} and {y} <= {ylim[1]}')
        if z_lower_cuttoff is not None:
            data = data.query(f'{z} >= {z_lower_cuttoff}')
        if log_scale:
            data[z] = np.log10(data[z])
            z_label = 'log(' + z_label + ')'

        figure = go.Figure()   
        figure.add_trace(go.Contour(x=data[x], y=data[y], z=data[z], colorscale=colorscale, contours_showlines=False, ncontours=ncontours))

        if title is not None:
            figure.update_layout(title_text=title)
        if width is not None:
            figure.update_layout(width=width)
        if height is not None:
            figure.update_layout(height=height)
        if x_label is not None:
            figure.update_xaxes(title_text=x_label)
        if y_label is not None:
            figure.update_yaxes(title_text=y_label)
        if template is not None:
            figure.update_layout(template=template)
        if z_label is not None:
            figure.update_traces(colorbar={'title': z_label})

        return figure
        

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

