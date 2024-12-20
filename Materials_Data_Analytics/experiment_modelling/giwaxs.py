from Materials_Data_Analytics.experiment_modelling.core import ScatteringMeasurement
import os
import pandas as pd
import numpy as np
import pyFAI
import pygix
import pickle
from PIL import Image
from datetime import datetime
import re
import plotly.graph_objects as go
import plotly.express as px


class Calibrator():
    ''' 
    A class to store the calibration parameters of a diffraction experiment 

    Main contributors:
    Arianna Magni

    Contributors:
    Nicholas Siemons
    '''
    def __init__(self, 
                 distance: float, 
                 poni1: float, 
                 poni2: float, 
                 rot1: float = 0,
                 rot2: float = 0,
                 rot3: float = 0,
                 energy: float = None,
                 wavelength: float = None,
                 detector = None):
        """Create a calibration object
        :param distance: sample-detector distance in meters
        :param poni1: coordinate of the point of normal incidence on the detector in the detector plane
        :param poni2: coordinate of the point of normal incidence on the detector in the detector plane
        :param rot1: rotation angle around the beam in radians
        :param rot2: rotation angle around the detector in radians
        :param rot3: rotation angle around the normal to the detector in radians
        :param energy: energy of the X-ray beam in eV
        :param wavelength: wavelength of the X-ray beam in meters
        :param detector: detector object or string
        """
        
        if isinstance(detector, str):
            self._detector = pyFAI.detector_factory(detector)
        else:
            self._detector = detector
        
        self._distance = distance
        self._poni1 = poni1
        self._poni2 = poni2
        self._rot1 = rot1
        self._rot2 = rot2
        self._rot3 = rot3

        if energy is not None:
            self._energy = energy
            self._wavelength = 1.239842e-6 / energy
        elif wavelength is not None:
            self._energy = 1.239842e-6 / wavelength
            self._wavelength = wavelength
        else:
            raise ValueError('One of energy or wavelength must be provided')
        
        self._azimuthal_integrator = self._make_azimuthal_integrator()

    @property
    def energy(self):
        return self._energy
    
    @property
    def wavelength(self):
        wavelength_nm = self._wavelength * 1e9
        return np.round(wavelength_nm, 5)
    
    @property
    def detector(self):
        return self._detector
    
    @property
    def distance(self):
        return np.round(self._distance, 5)
    
    @property
    def poni1(self):
        return np.round(self._poni1, 5)
    
    @property
    def poni2(self):
        return np.round(self._poni2, 5)
    
    @property
    def rot1(self):
        return np.round(self._rot1, 7)
    
    @property
    def rot2(self):
        return np.round(self._rot2, 7)
    
    @property
    def rot3(self):
        return np.round(self._rot3, 7)

    @classmethod
    def from_poni_file(cls, poni_file) -> 'Calibrator':
        """Create a calibration object from a .poni file
        :param poni_file: path to the .poni file
        :return: an instance of the Calibrator class
        """
        poni = pyFAI.load(poni_file)

        return cls(distance = poni.dist,
                   poni1 = poni.poni1,
                   poni2 = poni.poni2,
                   rot1 = poni.rot1,
                   rot2 = poni.rot2,
                   rot3 = poni.rot3,
                   detector = poni.detector,
                   wavelength = poni.wavelength)

    def save_to_pickle(self, pickle_file: str) -> 'Calibrator':
        """Save the calibration object to a pickle file
        :param pickle_file: path to the pickle file
        :return: the calibrator object
        """
        with open(pickle_file, 'wb') as file:
            pickle.dump(self, file)
        return self
    
    def _make_azimuthal_integrator (self) -> pyFAI.AzimuthalIntegrator:
        """
        Function to return an Azimuthal Integrator class from the pyFAI class
        """
        
        return pyFAI.azimuthalIntegrator.AzimuthalIntegrator(dist=self._distance, poni1=self._poni1, poni2=self._poni2,
                                                             rot1=self._rot1, rot2=self._rot2, rot3=self._rot3, detector=self._detector, 
                                                             wavelength=self._wavelength)
    

class GIWAXSPixelImage(ScatteringMeasurement):
    ''' 
    A class to store a GIWAXS measurement 

    Main contributors:
    Arianna Magni

    Contributors:
    Nicholas Siemons
    '''
    def __init__(self,
                 image : np.ndarray,
                 incidence_angle : float,
                 exposure_time : float,
                 timestamp : datetime,
                 number_of_averaged_images : int = 1,
                 metadata: dict = None):

        super().__init__(metadata=metadata)
        self._image = image
        self._incidence_angle = incidence_angle
        self._exposure_time = exposure_time
        self._timestamp = timestamp
        self._number_of_averaged_images = number_of_averaged_images
        self._mask = None

    @property
    def incidence_angle(self):
        return self._incidence_angle
    
    @property
    def exposure_time(self):
        return self._exposure_time
    
    @property
    def timestamp(self):
        return self._timestamp
    
    @property
    def number_of_averaged_images(self):
        return self._number_of_averaged_images
    
    @property
    def meta_data(self):
        return self._meta_data
            
    @property
    def image(self):
        return self._image

    @staticmethod   
    def _get_SLAC_BL11_3_parameters(txt_filepath: str) -> pd.DataFrame:
        '''
        Read the txt files from SLAC BL11-3 beamline and return a pandas DataFrame
        :param txt_filepaths: list of filepaths to the txt files
        :return: a pandas DataFrame with temperature, exposure time, i0, and monitor intensity
        '''       
        with open(txt_filepath, "r") as file:
            text = file.read()
            timestamp_str = re.search(r"time:\s*(.*)", text).group(1)
            timestamp = datetime.strptime(timestamp_str, "%a %b %d %H:%M:%S %Y")

            try:
                temperature = float(re.search(r"CTEMP=([\d.]+)", text).group(1))
            except:
                temperature = None
                                    
            incidence_angle = float(re.search(r" th=([\d.]+)", text).group(1))
            exposure_time = float(re.search(r"sec=([\d.]+)", text).group(1))
            intensity_norm = float(re.search(r"i0=([\d.]+)", text).group(1))
            monitor = float(re.search(r"mon=([\d.]+)", text).group(1))

        return {'timestamp': timestamp,
                'incidence_angle_deg': incidence_angle,
                'exposure_time_s': exposure_time,
                'intensity_norm': intensity_norm,
                'monitor': monitor,
                'temperature_c': temperature}
    
    @classmethod
    def from_SLAC_BL11_3(cls, 
                         tif_filepaths: list[str] | str  = None, 
                         txt_filepaths: list[str] | str = None,
                         verbose: bool = False,
                         metadata: dict = {}) -> 'GIWAXSPixelImage':
        
        """Load a GIWAXS measurement from SLAC BL11-3 beamline
        :param tif_filepaths: list of filepaths to the tif files
        :param txt_filepaths: list of filepaths to the txt files
        :param verbose: whether to print the output
        :param metadata: metadata to be stored with the measurement
        :return: an instance of the GIWAXSMeasurement class
        """     
        if txt_filepaths is None:
            txt_filepaths = [os.path.splitext(f)[0]+'.txt' for f in tif_filepaths if f.endswith('.tif')]
            if verbose: print(f'Metadata will be extracted from {txt_filepaths}')

        if isinstance(tif_filepaths, str):
            tif_filepaths = [tif_filepaths]
        if isinstance(txt_filepaths, str):
            txt_filepaths = [txt_filepaths]

        data = (pd
                .DataFrame({'txt_filepath': txt_filepaths, 'img_filepath': tif_filepaths})
                .groupby(['txt_filepath', 'img_filepath'], group_keys=True)
                .apply(lambda df: df.assign(param_dict = lambda x: [cls._get_SLAC_BL11_3_parameters(x) for x in x['txt_filepath']]))
                .assign(
                    timestamp = lambda x: [d['timestamp'] for d in x['param_dict']],
                    incidence_angle_deg = lambda x: [d['incidence_angle_deg'] for d in x['param_dict']],
                    exposure_time_s = lambda x: [d['exposure_time_s'] for d in x['param_dict']],
                    intensity_norm = lambda x: [d['intensity_norm'] for d in x['param_dict']],
                    monitor = lambda x: [d['monitor'] for d in x['param_dict']],
                    temperature_c = lambda x: [d['temperature_c'] for d in x['param_dict']]
                )
                .drop(columns=['param_dict'])
                )

        image, incidence_angle, exposure_time, N = cls._average_multiple_tif_files(tif_filepaths,
                                                                                   data['intensity_norm'].to_list(),
                                                                                   data['exposure_time_s'].to_list(),
                                                                                   data['incidence_angle_deg'].to_list(),
                                                                                   verbose=verbose)
        
        timestamp = data['timestamp'].min()    
        metadata['instrument_parameters'] = data
        metadata['source'] = 'SLAC_BL11_3'

        return cls(image,
                   incidence_angle,
                   exposure_time,
                   timestamp,
                   metadata =metadata,
                   number_of_averaged_images = N)

    @staticmethod
    def _load_tif_file(filepath: str) -> np.ndarray:
        """Load a TIFF file and return it as a NumPy array
        :param filepath: path to the TIFF file
        :return: the image data as a np.ndarray
        """
        with Image.open(filepath) as img:
            return np.array(img)
               
    @staticmethod
    def _average_multiple_tif_files(image_file_list : list[str],
                                    intensity_norm_list : list[float],
                                    exposure_time_list : list[float],
                                    incidence_angle_list : list[float],
                                    verbose: bool = False)  -> np.ndarray:
        
        """ Average multiple tif files and return the averaged image
        :param image_file_list: list of filepaths to the tif files
        :param intensity_list: list of intensities
        :param exposure_time_list: list of exposure times
        :param incidence_angle_list: list of incidence angles
        :param print_output: whether to print the output
        :return: a tuple containing the averaged image as a NumPy array, incidence angle as a float, exposure time as a float, and the number of averaged images as an int
        """
        ## Load tiff file and returns it as a np.array
        if verbose:
            print(" --- List of selected images ")
            for image_file in image_file_list:
                print(image_file)

        if len(set(exposure_time_list)) != 1:
            raise Warning("Not all files have the same exposure time. The exposure time will be averaged.")
        else:
            exposure_time = np.mean(exposure_time_list)

        if len(set(incidence_angle_list)) != 1:
            raise ValueError('Not all files have the same incidence angle. Files cannot be averaged.')
        else:
            incidence_angle = incidence_angle_list[0]            

        images_list = []

        for image_file, intensity_norm in zip(image_file_list, intensity_norm_list):
            image_data = GIWAXSPixelImage._load_tif_file(image_file)
            image_data_norm = (image_data/intensity_norm) *np.mean(intensity_norm_list)
            images_list.append(image_data_norm)

        # Convert the list of images to a NumPy array
        images_array = np.array(images_list)

        # Calculate the average over all the images
        image_data_average = np.squeeze(np.mean(images_array, axis=0))
        N = len(image_file_list)

        return image_data_average, incidence_angle, exposure_time, N 
    

    def apply_mask(self, mask_path: str) -> 'GIWAXSPixelImage':
        """ Apply a mask to the image.
        :param mask_path: path to the mask file
        :return: the masked image
        """   
        img = self._image
        mask = GIWAXSPixelImage._load_tif_file(mask_path)
        self._mask = mask
        img_masked = np.where(mask == 1, np.nan, img)
        self._image = img_masked
        self._image_original = img
        self.metadata['mask_path'] = mask_path
        return self
       
    def get_giwaxs_pattern(self,
                           calibrator: Calibrator,
                           qxy_range = (-3, 3),
                           qz_range = (0, 3),
                           q_range = (0, 3),
                           chi_range = (-95, 95),
                           pixel_q: int = 500,
                           pixel_chi: int = 360,
                           correct_solid_angle: bool = True,
                           polarization_factor: bool = None,
                           unit: str = 'A') -> 'GIWAXSPattern':
        """Transform the data from pixels to q space.
        :param calibrator: the calibrator object
        :param qxy_range: range of qxy values
        :param qz_range: range of qz values
        :param q_range: range of q values
        :param chi_range: range of chi values
        :param pixel_q: number of pixels in q
        :param pixel_chi: number of pixels in chi
        :param correct_solid_angle: whether to correct for solid angle
        :param polarization_factor: polarization factor
        :param unit: unit of the q values
        :param verbose: whether to print the output
        :return: an instance of the GIWAXSPattern class
        """
        azimuthal_integrator = calibrator._azimuthal_integrator
        transformer = pygix.transform.Transform().load(azimuthal_integrator)
        transformer.incident_angle = np.deg2rad(self.incidence_angle)

        [intensity_reciprocal, qxy, qz] = transformer.transform_reciprocal(self._image,
                                                                           npt = (pixel_q, pixel_q),
                                                                           ip_range = qxy_range,
                                                                           op_range = (-qz_range[0], -qz_range[1]),
                                                                           method = 'splitbbox',
                                                                           unit = unit,
                                                                           correctSolidAngle = correct_solid_angle,
                                                                           polarization_factor = polarization_factor)

        qz = -qz

        pixel_chi_corr = int(pixel_chi*360/(chi_range[1] - chi_range[0]))

        [intensity_polar, q, chi] = transformer.transform_polar(self._image,
                                                                npt = (pixel_q, pixel_chi_corr),
                                                                q_range = q_range,
                                                                chi_range = (-180, 180),
                                                                correctSolidAngle = correct_solid_angle,
                                                                polarization_factor = polarization_factor,
                                                                unit = unit,
                                                                method = 'splitbbox')
        
        chi = np.where(chi > 0, -chi + 180, -chi - 180)

        return GIWAXSPattern.from_numpy_arrays(qxy = qxy, 
                                               qz = qz, 
                                               intensity_reciprocal = intensity_reciprocal, 
                                               chi = chi, 
                                               q = q, 
                                               intensity_polar = intensity_polar, 
                                               metadata = self.metadata)

    def save_to_pickle(self, pickle_file: str) -> 'GIWAXSPixelImage':
        """Save the GIWAXS measurement to a pickle file
        :param pickle_file: path to the pickle file
        :return: the GIWAXS measurement object
        """
        with open(pickle_file, 'wb') as file:
            pickle.dump(self, file)
        return self
    
        
class GIWAXSPattern(ScatteringMeasurement):
    ''' 
    A class to store a GIWAXS measurement 

    Main contributors:
    Arianna Magni

    Contributors:
    Nicholas Siemons
    '''

    def __init__(self,
                 data_reciprocal: pd.DataFrame = None,
                 data_polar: pd.DataFrame = None,
                 metadata: dict = None):
        
        super().__init__(metadata=metadata)
        self._data_reciprocal = data_reciprocal
        self._data_polar = data_polar
        
    @classmethod
    def from_numpy_arrays(cls,
                          qxy: np.ndarray = None,
                          qz: np.ndarray = None,
                          intensity_reciprocal: np.ndarray = None,
                          chi: np.ndarray = None,
                          q: np.ndarray = None,
                          intensity_polar: np.ndarray = None,
                          metadata: dict = None):
        """
        Create a GIWAXSPattern object from numpy arrays
        """
        data_reciprocal = (pd
                           .DataFrame(intensity_reciprocal, columns=qxy, index=qz)
                           .reset_index()
                           .melt(id_vars='index')
                           .rename(columns={'index': 'qz', 'variable': 'qxy', 'value': 'intensity'})
                           )

        data_polar = (pd
                     .DataFrame(intensity_polar, columns=q, index=chi)
                     .reset_index()
                     .melt(id_vars='index')
                     .rename(columns={'index': 'chi', 'variable': 'q', 'value': 'intensity'})
                     )

        return cls(data_reciprocal = data_reciprocal, data_polar = data_polar, metadata = metadata)

        
    @property
    def data_reciprocal(self):
        return self._data_reciprocal
    
    @property
    def qxy(self):
        qxy = self._data_reciprocal.sort_values(by='qxy')['qxy'].unique()
        return qxy
    
    @property
    def qz(self):
        qz = self._data_reciprocal.sort_values(by='qz')['qz'].unique()
        return qz
        
    @property
    def data_polar(self):
        return self._data_polar
    
    @property
    def chi(self):
        chi = self._data_polar.sort_values(by='chi')['chi'].unique()
        return chi
    
    @property
    def q(self):
        q = self._data_polar.sort_values(by='q')['q'].unique()
        return q
    
    @property
    def meta_data(self):
        return self._meta_data
    
    def export_reciprocal_data(self, export_filepath: str, format: str = 'wide') -> 'GIWAXSPattern':
        """Export the reciprocal space data to a CSV file.
        :param export_filepath: Filepath to export the data to.
        :param format: Format of the data. Either 'long' or 'wide'.
        :return: the current instance
        """
        if format not in ['long', 'wide']:
            raise ValueError('format must be either "long" or "wide"')

        directory = os.path.dirname(export_filepath)
        if directory and not os.path.exists(directory):
            os.makedirs(directory)
            print(f"Directory {directory} created.")
        elif os.path.exists(export_filepath):
            print(f"File {export_filepath} already exists. It will be overwritten.")

        if format == 'long':
            pd.to_csv(self._data_reciprocal, export_filepath)

        return self
    
    def export_polar_data(self, export_filepath: str, format: str = 'wide') -> 'GIWAXSPattern':
        """Export the polar space data to a CSV file.
        :param export_filepath: Filepath to export the data to.
        :param format of the data, long or wide
        :return: the current instance
        """
        if format not in ['long', 'wide']:
            raise ValueError('format must be either "long" or "wide"')
        
        directory = os.path.dirname(export_filepath)
        if directory and not os.path.exists(directory):
            os.makedirs(directory)
            print(f"Directory {directory} created.")
        elif os.path.exists(export_filepath):
            print(f"File {export_filepath} already exists. It will be overwritten.")

        if format == 'long':
            pd.to_csv(self._data_polar, export_filepath)

        return self
    
    def save_to_pickle(self, pickle_file: str) -> 'GIWAXSPattern':
        """Save the GIWAXS measurement to a pickle file
        :param pickle_file: path to the pickle file
        :return: the GIWAXS measurement object
        """
        with open(pickle_file, 'wb') as file:
            pickle.dump(self, file)
        return self
    
    def plot_reciprocal_map_contour(self, 
                                    colorscale: str = 'blackbody', 
                                    ncontours: int = 100, 
                                    log_scale: bool = True, 
                                    template: str = 'simple_white',
                                    intensity_lower_cuttoff: float = 0.001,
                                    **kwargs) -> go.Figure:
        """Plot the reciprocal space map.
        :param colorscale: The colorscale to use. See plotly colorscales for options
        :param ncontours: The number of contours to use.
        :param log_scale: Whether to use a log scale.
        :param template: The template to use.
        :param intensity_lower_cuttoff: The lower cuttoff for the intensity. Useful if using log values.
        :return: The plot.
        """
        data = self.data_reciprocal.copy()
        fig = self.plot_contour_map(data = data, x='qxy', y='qz', z='intensity', colorscale=colorscale, ncontours=ncontours, z_lower_cuttoff=intensity_lower_cuttoff,
                                    template=template, x_label='qxy [\u212B\u207B\u00B9]', y_label='qz [\u212B\u207B\u00B9]', log_scale=log_scale,
                                    z_label='Intensity', **kwargs)
        return fig
    
    def plot_reciprocal_map(self, 
                            colorscale: str = 'blackbody',
                            log_scale: bool = True,  
                            template: str = 'simple_white',
                            origin: str = 'lower',
                            intensity_lower_cuttoff: float = 0.001,
                            **kwargs) -> go.Figure:
        """Plot the reciprocal space map.
        :param colorscale: The colorscale to use. See plotly colorscales for options
        :param log_scale: Whether to use a log scale.
        :param template: The template to use.
        :param origin: The origin zero point of the plot.
        :param intensity_lower_cuttoff: The lower cuttoff for the intensity. Useful if using log values.
        :return: The plot.
        """
        data = self.data_reciprocal.copy().sort_values(by=['qxy', 'qz'], ascending=[True, False])
        fig = self.plot_pixel_map(data = data, x='qxy', y='qz', z='intensity', colorscale=colorscale, log_scale=log_scale, z_lower_cuttoff=intensity_lower_cuttoff,
                            x_label='qxy [\u212B\u207B\u00B9]', y_label='qz [\u212B\u207B\u00B9]', template=template, origin=origin,
                            z_label='Intensity', **kwargs)
        return fig
    
    def plot_polar_map_contour(self, 
                               colorscale: str = 'blackbody', 
                               ncontours: int = 100, 
                               log_scale: bool = True,
                               template: str = 'simple_white',
                               intensity_lower_cuttoff: float = 0.001,
                               **kwargs) -> go.Figure:
        """Plot the polar space map.
        :param colorscale: The colorscale to use. See plotly colorscales for options
        :param ncontours: The number of contours to use.
        :param log_scale: Whether to use a log scale.
        :param template: The template to use.
        :param intensity_lower_cuttoff: The lower cuttoff for the intensity. Useful if using log values.
        :return: The plot.
        """
        data = self.data_polar.copy()
        fig = self.plot_contour_map(data = data, y='chi', x='q', z='intensity', colorscale=colorscale, ncontours=ncontours, log_scale=log_scale, 
                                    z_lower_cuttoff=intensity_lower_cuttoff, template=template, x_label='q [\u212B\u207B\u00B9]', y_label='chi [\u00B0]', 
                                    z_label='Intensity', **kwargs)
        return fig
    
    def plot_polar_map(self, 
                       colorscale: str = 'blackbody', 
                       log_scale: bool = True,
                       template: str = 'simple_white',
                       origin: str = 'lower',
                       intensity_lower_cuttoff: float = 0.001,
                       **kwargs):
        """Plot the polar space map.
        :param colorscale: The colorscale to use. See plotly colorscales for options
        :return: The plot.
        """
        data = self.data_polar.copy().sort_values(by=['q', 'chi'], ascending=[True, False])
        fig = self.plot_pixel_map(data = data, y='chi', x='q', z='intensity', colorscale=colorscale, aspect='auto', z_lower_cuttoff=intensity_lower_cuttoff,
                                  origin=origin, log_scale=log_scale,x_label='Q [\u212B\u207B\u00B9]', y_label='chi [\u00B0]', 
                                  z_label='Intensity', template=template, **kwargs)
        return fig
    
    def get_linecut(self,
                    chi : tuple | list | pd.Series | float = None,
                    q_range : tuple | list | pd.Series = None) -> pd.DataFrame:
        """Extract a profile from the polar space data.
        :param chi: Range of chi values or a single chi value.
        :param q_range: q_range.
        """
        data = self.data_polar.copy()

        # check if chi is iterable
        try:
            iter(chi)
            chi_iterable = True
            if len(chi) != 2: raise ValueError('If chi is a range it must be two values')
        except TypeError:
            chi_iterable = False
            if chi < data['chi'].min() or chi > data['chi'].max(): raise ValueError('chi value out of range of the data')

        # Filter the data for chi
        if chi_iterable: 
            data = data.query(f'chi >= {min(chi)} and chi <= {max(chi)}')
        else:
            closest_index = data['chi'].sub(chi).abs().idxmin()
            closest_chi = data.loc[closest_index, 'chi']
            data = data.query(f'chi == {closest_chi}')

        # Filter the data for q
        if q_range is not None: 
            data = data.query(f'q >= {min(q_range)} and q <= {max(q_range)}')
        
        data = data.groupby('q').mean().reset_index().filter(['chi', 'q', 'intensity'])

        return data
    
    def plot_linecut(self,
                     chi: tuple | list | pd.Series | float = None,
                     q_range: tuple | list | pd.Series = None, 
                     **kwargs) -> px.line:

        """Plot a profile extracted from the polar space data.
        :param chi: Range of chi values or a single chi value.
        :param q_range: q_range.
        :return: The plot.
        """
        profile = self.get_linecut(chi, q_range)
        figure = px.line(profile, x='q', y='intensity', labels={'q': 'q [\u212B\u207B\u00B9]', 'intensity': 'Intensity'}, **kwargs)
        return figure
