from Materials_Data_Analytics.experiment_modelling.core import ScatteringMeasurement
import os
#import fabio
import pandas as pd
import numpy as np
import pyFAI
import pygix
import pickle
#import holoviews as hv
from PIL import Image
from datetime import datetime
import re
import holoviews as hv


class Calibrator():
    ''' A class to store the calibration parameters of a diffraction experiment '''
    def __init__(self, 
                 distance: float, 
                 poni1: float, 
                 poni2: float, 
                 rot1: float = 0,
                 rot2: float = 0,
                 rot3: float = 0,
                 energy: float = None,
                 wavelegth: float = None,
                 source: str = None, 
                 detector = None):
        
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
        self._source = source

        if energy is not None:
            self._energy = energy
            self._wavelegth = 12.39842 / energy
        elif wavelegth is not None:
            self._energy = 12.39842 / wavelegth
            self._wavelegth = wavelegth
        else:
            raise ValueError('One of energy or wavelegth must be provided')
        
        self._ai = self._make_ai()

    @property
    def energy(self):
        return self._energy
    
    @property
    def wavelegth(self):
        return self._wavelegth
    
    @property
    def detector(self):
        return self._detector
    
    @property
    def distance(self):
        return self._distance
    
    @property
    def poni1(self):
        return self._poni1
    
    @property
    def poni2(self):
        return self._poni2
    
    @property
    def rot1(self):
        return self._rot1
    
    @property
    def rot2(self):
        return self._rot2
    
    @property
    def rot3(self):
        return self._rot3
    
    @property
    def source(self):
        return self._source
    
    @property
    def ai(self):
        return self._ai

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
                   wavelegth = poni.wavelength,
                   source = poni_file)
  

    def save_to_pickle(self, pickle_file: str) -> 'Calibrator':
        """Save the calibration object to a pickle file
        :param pickle_file: path to the pickle file
        :return: the calibrator object
        """
        with open(pickle_file, 'wb') as file:
            pickle.dump(self, file)
        return self
    
    def _make_ai (self) -> pyFAI.AzimuthalIntegrator:
        return pyFAI.AzimuthalIntegrator(dist=self.distance, poni1=self.poni1, poni2=self.poni2,
                                         rot1=self.rot1, rot2=self.rot2, rot3=self.rot3, detector=self.detector, 
                                         wavelength=self.energy)
    

class GIWAXSMeasurementPixelImage(ScatteringMeasurement):
    ''' A class to store a GIWAXS measurement '''

    def __init__(self,
                 row_image : np.ndarray,
                 incidence_angle : float,
                 exposure_time : float,
                 timestamp : datetime,
                 number_of_averaged_images : int = 1,
                 metadata: dict = None):

        super().__init__(metadata=metadata)
        self._image_row = row_image
        self._incidence_angle = incidence_angle
        self._exposure_time = exposure_time
        self._timestamp = timestamp
        self._number_of_averaged_images = number_of_averaged_images

 
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
    def image_row(self):
        return self._image_row
    
    @property
    def data_reciprocal(self):
        return self._data_reciprocal
    
    @property
    def data_polar(self):
        return self._data_polar
    
    @property
    def image_masked(self):
        if hasattr(self, '_row_image_masked'):
            return self._image_masked
        else:
            raise ValueError('No mask applied')
        

    @staticmethod   
    def _read_txt_file_SLAC_BL11_3(txt_filepaths: list[str]) -> pd.DataFrame:
        '''Read the txt files from SLAC BL11-3 beamline and return a pandas DataFrame
        :param txt_filepaths: list of filepaths to the txt files
        :return: a pandas DataFrame with temperature, exposure time, i0, and monitor intensity
        '''       

        timestamp_list = []
        temperature_list = []
        exposure_time_list = []
        incidence_angle_list = []
        i0_list = []
        monitor_list = []

        for filepath in txt_filepaths:  
            with open(filepath, "r") as file:
                text = file.read()
                
                timestamp_str = re.search(r"time:\s*(.*)", text).group(1)
                timestamp_list.append(datetime.strptime(timestamp_str, "%a %b %d %H:%M:%S %Y")) 
                try:
                    temperature = float(re.search(r"CTEMP=([\d.]+)", text).group(1))
                except:
                    try:
                        temperature = float(re.search(r"CTEMP=([\d.]+)", text).group(1))
                    except:
                        temperature = None
                                       
                temperature_list.append(temperature)
                incidence_angle_list.append(float(re.search(r" th=([\d.]+)", text).group(1)))
                exposure_time_list.append(float(re.search(r"sec=([\d.]+)", text).group(1)))
                i0_list.append(float(re.search(r"i0=([\d.]+)", text).group(1)))
                monitor_list.append(float(re.search(r"mon=([\d.]+)", text).group(1)))

        return pd.DataFrame({
            'txt filepath': txt_filepaths,
            'timestamp': timestamp_list,
            'incidence angle [deg]': incidence_angle_list,
            'Exposure time [s]': exposure_time_list,
            'i0': i0_list,
            'monitor': monitor_list,
            'Temperature [C]': temperature_list})
    
   
    @classmethod
    def from_SLAC_BL11_3(cls,
        tif_filepaths: list[str] | str  = None,
        txt_filepaths: list[str] | str = None,
        verbose: bool = False,
        metadata: dict = {}) -> 'GIWAXSMeasurementPixelImage':
        """Load a GIWAXS measurement from SLAC BL11-3 beamline
        :param tif_filepaths: list of filepaths to the tif files
        :param txt_filepaths: list of filepaths to the txt files
        :return: an instance of the GIWAXSMeasurement class
        """
                
        if isinstance(tif_filepaths, str):
            tif_filepaths = [tif_filepaths]
        if isinstance(txt_filepaths, str):
            txt_filepaths = [txt_filepaths]

        if txt_filepaths is None:
            txt_filepaths = [os.path.splitext(f)[0]+'.txt' for f in tif_filepaths if f.endswith('.tif')]
            if verbose:
                print(f'Metadata will be extracted from {txt_filepaths}')
        
        property_df = cls._read_txt_file_SLAC_BL11_3(txt_filepaths)
        
        row_image, incidence_angle, exposure_time, N = cls._average_multiple_tif_files(tif_filepaths,
                                                                  property_df['i0'].to_list(),
                                                                  property_df['Exposure time [s]'].to_list(),
                                                                  property_df['incidence angle [deg]'].to_list(),
                                                                  verbose=verbose)
        source_df = pd.DataFrame({
            'img filepath': tif_filepaths,
            'txt filepath': txt_filepaths})
        
        timestamp = property_df['timestamp'].min()    

        metadata['instrument_parameters'] = pd.concat([property_df, source_df], axis=1)
        metadata['source'] = 'SLAC BL11-3'

        return cls(row_image,
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
    def _average_multiple_tif_files (image_file_list : list[str],
                                    intensity_list : list[float],
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
            print("   --- List of selected images ")
            for image_file in image_file_list:
                print(image_file)

        if len(set(exposure_time_list)) != 1:
            print('Not all files have the same exposure. Files will be averaged anyway.')
        else:
            exposure_time = exposure_time_list[0]

        if len(set(incidence_angle_list)) != 1:
            raise ValueError('Not all files have the same incidence angle. Files cannot be averaged.')
        else:
            incidence_angle = incidence_angle_list[0]            

        images_list = []

        for image_file, intensity in zip(image_file_list, intensity_list):
 
            image_data = GIWAXSMeasurementPixelImage._load_tif_file(image_file)
            image_data_norm = image_data/intensity
            # Append the image data to the list
            images_list.append(image_data_norm)

        # Convert the list of images to a NumPy array
        images_array = np.array(images_list)

        # Calculate the average over all the images
        image_data_average = np.squeeze(np.mean(images_array, axis=0))
        N = len(image_file_list)

        return image_data_average, incidence_angle, exposure_time, N 
    

    def apply_mask(self, maskpath) -> 'GIWAXSMeasurementPixelImage':
        """ Apply a mask to the image.
        :param maskpath: path to the mask file
        :return: the masked image
        """
                
        img = self._image_row
        
        mask = GIWAXSMeasurementPixelImage._load_tif_file(maskpath)
        self._mask = mask

        img_masked = np.where(mask == 1, np.nan, img)
        self._image_masked = img_masked
        self.metadata['maskpath'] = maskpath
        
        return self
       

    def transform (
        self,
        calibrator: Calibrator,
        qxy_range = (-3, 3),
        qz_range = (0, 3),
        q_range = (0, 3),
        chi_range = (-95, 95),
        pixel_Q: int = 500,
        pixel_chi: int = 360,
        correct_solid_angle=True,
        polarization_factor=None,
        unit='A') -> 'GIWAXSPattern':
    
        """Transform the data from pixels to q space.
        Args:
        - poni_file_path (str): Filepath to the PONI file for calibration.
        - incidence_angle (float): Incident angle.
        - pixel_Q (int): Pixel for Q.
        - pixel_Chi (int): Pixel for Chi.

        Returns:
        - self: The current instance.
        """
        ai = calibrator.ai
        pg = pygix.transform.Transform(
            dist=ai.dist,
            poni1=ai.poni1,
            poni2=ai.poni2,
            rot1=ai.rot1,
            rot2=ai.rot2,
            rot3=ai.rot3,
            pixel1=ai.pixel1,
            pixel2=ai.pixel2,
            wavelength=ai.wavelength)
        
        pg.incident_angle = self.incidence_angle
        
        if hasattr(self, '_row_image_masked'):
            image = self._image_masked
        else:
            print('No mask applied')
            image = self._image_row

        
        IntQ, qxy, qz = pg.transform_reciprocal(image,
                                                npt = (pixel_Q, pixel_Q),
                                                ip_range = qxy_range,
                                                op_range = qz_range,
                                                method = 'splitbbox',
                                                unit = unit,
                                                correctSolidAngle = correct_solid_angle,
                                                polarization_factor = polarization_factor)
    
             

        IntChiQ, Q, chi = pg.transform_polar (image,
                                              npt = (pixel_Q, pixel_chi),
                                              q_range = q_range,
                                              chi_range = chi_range,
                                              correctSolidAngle = correct_solid_angle,
                                              polarization_factor = polarization_factor,
                                              unit = unit,
                                              method = 'splitbbox')
        

        data_reciprocal = pd.DataFrame({
            'qxy':  np.repeat(qxy, len(qz)),  # Repeat qxy values
            'qz': np.tile(qz, len(qxy)),       # Tile qz values
            'I': IntQ.flatten()
        })

        data_polar = pd.DataFrame({
            'Q':  np.repeat(Q, len(chi)),  # Repeat Q values
            'chi': np.tile(chi, len(Q)),       # Tile chi values
            'I': IntChiQ.flatten()
        })

        self._data_reciprocal = data_reciprocal
        self._data_polar = data_polar
        self._ai = ai    
   
        return self
    
    def to_GIWAXSPattern(self) -> 'GIWAXSPattern':
        """Convert the GIWAXS measurement to a GIWAXS pattern.
        :return: the GIWAXS pattern
        """
        return GIWAXSPattern(data_reciprocal = self._data_reciprocal,
                             data_polar = self._data_polar,
                             metadata = self._metadata)

    def save_to_pickle(self, pickle_file: str) -> 'GIWAXSMeasurementPixelImage':
        """Save the GIWAXS measurement to a pickle file
        :param pickle_file: path to the pickle file
        :return: the GIWAXS measurement object
        """
        with open(pickle_file, 'wb') as file:
            pickle.dump(self, file)
        return self
    
        

class GIWAXSPattern(ScatteringMeasurement):
    ''' A class to store a GIWAXS measurement '''

    def __init__(self,
                 data_reciprocal: pd.DataFrame = None,
                 data_polar: pd.DataFrame = None,
                 qxy: np.ndarray = None,
                 qz: np.ndarray = None,
                 Intensity_reciprocal: np.ndarray = None,
                 chi: np.ndarray = None,
                 Q: np.ndarray = None,
                 Intensity_polar: np.ndarray = None,
                 metadata: dict = None):
        
        super().__init__(metadata=metadata)
        
        if ((data_reciprocal is not None) & 
            (data_polar is not None)):
            
            self._data_reciprocal = data_reciprocal
            self._data_polar = data_polar
        
        elif ((qxy is not None) &
              (qz is not None) &
              (Intensity_reciprocal is not None) &
              (chi is not None) &
              (Q is not None) &
              (Intensity_polar is not None)):
            
            self._data_reciprocal = pd.DataFrame({
                'qxy':  np.repeat(qxy, len(qz)),  
                'qz': np.tile(qz, len(qxy)),       
                'I': Intensity_reciprocal.flatten()
            })

            self._data_polar = pd.DataFrame({
                'Q':  np.repeat(Q, len(chi)),
                'chi': np.tile(chi, len(Q)),  
                'I': Intensity_polar.flatten()
            })
        
        else:
            raise ValueError('Either data_reciprocal and data_polar or qxy, qz, Intensity_reciprocal, chi, Q, Intensity_polar must be provided')
                
        self._metadata = metadata

        
    @property
    def data_reciprocal(self):
        return self._data_reciprocal
    
    @property
    def qxy(self):
        return np.sort(self._data_reciprocal['qxy'].unique())
    
    @property
    def qz(self):
        return np.sort(self._data_reciprocal['qz'].unique())
    
    @property
    def Intensity_reciprocal(self):
        return self._data_reciprocal.pivot(index='qxy', columns='qz', values='I').to_numpy()
        
    @property
    def data_polar(self):
        return self._data_polar
    
    @property
    def chi(self):
        return np.sort(self._data_polar['chi'].unique())
    
    @property
    def Q(self):
        return np.sort(self._data_polar['Q'].unique())
    
    @property
    def Intensity_polar(self):
        return self._data_polar.pivot(index='Q', columns='chi', values='I').to_numpy()
    
    @property
    def meta_data(self):
        return self._meta_data
    
    def export_reciprocal_data(self, export_filepath: str) -> 'GIWAXSPattern':
        """Export the reciprocal space data to a CSV file.
        Args:
        - export_filepath (str): Filepath to export the data to.
        Returns:
        - self: The current instance.
        """
  
        qxy = np.sort(self.data_reciprocal['qxy'].unique())
        qz = np.sort(self.data_reciprocal['qz'].unique())
        Intensity = self.data_reciprocal.pivot(index='qxy', columns='qz', values='I').to_numpy()

        matrix = pd.DataFrame(Intensity, index=qz, columns=qxy)
        matrix.index.name = "qz"
        matrix.columns.name = "qxy"

        directory = os.path.dirname(export_filepath)
        if directory and not os.path.exists(directory):
            os.makedirs(directory)
            print(f"Directory {directory} created.")
        elif os.path.exists(export_filepath):
            print(f"File {export_filepath} already exists. It will be overwritten.")

        matrix.to_csv(export_filepath)

        return self
    

    def export_polar_data(self, export_filepath: str) -> 'GIWAXSPattern':
        """Export the polar space data to a CSV file.
        Args:
        - export_filepath (str): Filepath to export the data to.
        Returns:
        - self: The current instance.
        """

        Q = np.sort(self.data_polar['Q'].unique())
        chi = np.sort(self.data_polar['chi'].unique())
        Intensity = self.data_polar.pivot(index='Q', columns='chi', values='I').to_numpy()

        matrix = pd.DataFrame(Intensity, index=chi, columns=Q)
        matrix.index.name = "chi"
        matrix.columns.name = "Q"

        directory = os.path.dirname(export_filepath)
        if directory and not os.path.exists(directory):
            os.makedirs(directory)
            print(f"Directory {directory} created.")
        elif os.path.exists(export_filepath):
            print(f"File {export_filepath} already exists. It will be overwritten.")

        matrix.to_csv(export_filepath)
        
        return self
    
    def save_to_pickle(self, pickle_file: str) -> 'GIWAXSPattern':
        """Save the GIWAXS measurement to a pickle file
        :param pickle_file: path to the pickle file
        :return: the GIWAXS measurement object
        """
        with open(pickle_file, 'wb') as file:
            pickle.dump(self, file)
        return self
    
    def plot_map_reciprocal_space(self,
                                  title = r'Q map', 
                                  xlim : tuple = (-0.5, 2.3),
                                  ylim : tuple = (-0.5, 2.3),
                                  clim : tuple = None,
                                  binning : int = 1,  
                                  frame_width : int = 800,
                                  frame_height : int = 600,
                                  cmap : str = 'viridis',
                                  colorbar : bool = True,
                                  logz : bool = True) -> hv.Image:  
        """Plot the reciprocal space map.
        Args:
        - title (str): Title of the plot.
        - xlim (tuple): X-axis limits.
        - ylim (tuple): Y-axis limits.
        - clim (tuple): Color limits.
        - binning (int): Binning factor.
        - frame_width (int): Width of the frame.
        - frame_height (int): Height of the frame.
        - cmap (str): Colormap.
        - colorbar (bool): Whether to show the colorbar.
        - logz (bool): Whether to use a logarithmic scale for the z-axis.
        Returns:
        - hv.Image: The plot.
        """
                         
        qxy = self.qxy
        qz = self.qz
        Intensity = self.Intensity_reciprocal.T

        if binning > 1:
            qxy = qxy[::binning]
            qz = qz[::binning]
            Intensity_binned = np.zeros((len(qxy), len(qz)))
            for i in range(len(qxy)):
                for j in range(len(qz)):
                    Intensity_binned[i,j] = np.mean(Intensity[i*binning:(i+1)*binning, j*binning:(j+1)*binning])
            Intensity = Intensity_binned

        return hv.Image((qxy, qz, Intensity), kdims=['qxy', 'qz'], vdims=['Intensity']).opts(
            title=title,
            cmap=cmap,
            colorbar=colorbar,
            xlim = xlim,
            ylim = ylim,
            clim = clim,
            logz=logz,
            width=frame_width,
            height=frame_height,
            xlabel = 'Qxy [Å⁻¹]',
            ylabel = 'Qz [Å⁻¹]')
    
    def plot_map_polar_space(self,
                             title = r'Polar map', 
                             xlim : tuple = (0, 2.3),
                             ylim : tuple = (-90, 90),
                             clim : tuple = None,
                             binning : int = 1,  
                             frame_width : int = 800,
                             frame_height : int = 600,
                             cmap : str = 'viridis',
                             colorbar : bool = True,
                             logz : bool = True) -> hv.Image:  
        """Plot the polar space map.
        Args:
        - title (str): Title of the plot.
        - xlim (tuple): X-axis limits.
        - ylim (tuple): Y-axis limits.
        - clim (tuple): Color limits.
        - binning (int): Binning factor.
        - frame_width (int): Width of the frame.
        - frame_height (int): Height of the frame.
        - cmap (str): Colormap.
        - colorbar (bool): Whether to show the colorbar.
        - logz (bool): Whether to use a logarithmic scale for the z-axis.
        Returns:
        - hv.Image: The plot.
        """
        
        Q = self.Q
        chi = self.chi
        Intensity = self.Intensity_polar.T
        
        if binning > 1:
            Q = Q[::binning]
            chi = chi[::binning]
            Intensity_binned = np.zeros((len(Q), len(chi)))
            for i in range(len(Q)):
                for j in range(len(chi)):
                    Intensity_binned[i,j] = np.mean(Intensity[i*binning:(i+1)*binning, j*binning:(j+1)*binning])
            Intensity = Intensity_binned

        return hv.Image((Q, chi, Intensity), kdims=['Q', 'chi'], vdims=['Intensity']).opts(
            title=title,
            cmap=cmap,
            colorbar=colorbar,
            xlim = xlim,
            ylim = ylim,
            clim = clim,
            logz=logz,
            width=frame_width,
            height=frame_height,
            xlabel = 'Q [Å⁻¹]',
            ylabel = 'chi [°]') 
    
    def extract_profile (self,
                         chi_min : float = None,
                         chi_max : float = None,
                         q_range : tuple = None) -> pd.DataFrame:
        """Extract a profile from the polar space data.
        Args:
        - chi_min (float): Minimum chi value.
        - chi_max (float): Maximum chi value.
        - q_range (tuple): Range of Q values.
        Returns:
        - pd.DataFrame: The profile.
        """
    
        data = self.data_polar
        if chi_min is not None:
            data = data[data['chi'] >= chi_min]
        if chi_max is not None:
            data = data[data['chi'] <= chi_max]
        
        if q_range is not None:
            data = data[(data['Q'] >= q_range[0]) & (data['Q'] <= q_range[1])]
            
        profile = data.groupby('Q').mean().reset_index()

        return pd.DataFrame({
            'Q': profile['Q'],
            'I': profile['I']})
    
    def plot_profile (self,
                      chi_min : float = None,
                      chi_max : float = None,
                      q_range : tuple = None,
                      label : str = '',
                      xlim : tuple = None,
                      ylim : tuple = None,
                      frame_width : int = 600,
                      frame_height : int = 400,
                      color : str = 'blue',
                      logy : bool = False,
                      ) -> hv.Curve:

        """Plot a profile extracted from the polar space data.
        Args:
        - chi_min (float): Minimum chi value.
        - chi_max (float): Maximum chi value.
        - q_range (tuple): Range of Q values.
        - xlim (tuple): X-axis limits.
        - ylim (tuple): Y-axis limits.
        - frame_width (int): Width of the frame.
        - frame_height (int): Height of the frame.
        - color (str): Color of the curve.
        Returns:
        - hv.Curve: The plot.
        """

        profile = self.extract_profile(chi_min, chi_max, q_range)

        return hv.Curve(profile, kdims = ['Q'], vdims = ['I'], label = label).opts(
            line_width=2,
            color=color,
            width=frame_width,
            height=frame_height,
            xlim=xlim,
            ylim=ylim,
            xlabel='Q [Å⁻¹]',
            ylabel='Intensity [arb. units]',
            logy=logy)
    
    def export_profile (self,
                        export_filepath: str,
                        chi_min : float = None,
                        chi_max : float = None,
                        q_range : tuple = None) -> 'GIWAXSPattern':
        """Export a profile extracted from the polar space data to a CSV file.
        Args:
        - export_filepath (str): Filepath to export the data to.
        - chi_min (float): Minimum chi value.
        - chi_max (float): Maximum chi value.
        - q_range (tuple): Range of Q values.
        Returns:
        - self: The current instance.
        """
        
        profile = self.extract_profile(chi_min, chi_max, q_range)
        
        directory = os.path.dirname(export_filepath)
        if directory and not os.path.exists(directory):
            os.makedirs(directory)
            print(f"Directory {directory} created.")
        elif os.path.exists(export_filepath):
            print(f"File {export_filepath} already exists. It will be overwritten.")

        profile.to_csv(export_filepath)
        return self
                    