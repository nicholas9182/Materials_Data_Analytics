from Materials_Data_Analytics.experiment_modelling.core import ScatteringMeasurement
import os
import pandas as pd
import numpy as np
import pickle
from datetime import datetime
import re
import plotly.graph_objects as go
import plotly.express as px
import lmfit
import importlib.util


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
        if importlib.util.find_spec('pyFAI') is None:
            raise ImportError('pyFAI is required to run this function. Please install pyFAI using pip install pyFAI')
        else:
            import pyFAI

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
        if importlib.util.find_spec('pyFAI') is None:
            raise ImportError('pyFAI is required to run this function. Please install pyFAI using pip install pyFAI')
        else:
            import pyFAI

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
    
    def _make_azimuthal_integrator (self):
        """
        Function to return an Azimuthal Integrator class from the pyFAI class
        """
        if importlib.util.find_spec('pyFAI') is None:
            raise ImportError('pyFAI is required to run this function. Please install pyFAI using pip install pyFAI')
        else:
            import pyFAI

        return pyFAI.AzimuthalIntegrator(dist=self._distance, poni1=self._poni1, poni2=self._poni2,
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
    def metadata(self):
        return self._metadata
            
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
        if importlib.util.find_spec('PIL') is None:
            raise ImportError('PIL is required to run this function. Please install PIL using pip install PIL')
        else:
            from PIL import Image

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
            image_data_norm = (image_data/intensity_norm) * exposure_time
            images_list.append(image_data_norm)

        # Convert the list of images to a NumPy array
        images_array = np.array(images_list)

        # Calculate the average over all the images
        image_data_average = np.squeeze(np.mean(images_array, axis=0))
        N = len(image_file_list)

        return image_data_average, incidence_angle, exposure_time, N 
    
    @classmethod
    def from_NSLS_II_CMS(cls,
                         filepaths: list[str] | str  = None,
                         verbose: bool = False,
                         stiching_offset: int = 30,
                         timestamp: datetime = None,
                         metadata: dict = {})-> 'GIWAXSPixelImage':
        
        if isinstance(filepaths, list) and len(filepaths) == 1:
            filepaths = filepaths[0]

        if isinstance(filepaths, str):
            # single image
            metadata = cls._get_NSLS_II_CMS_parameters(filepaths, verbose=verbose)
            image = cls._load_tif_file(filepaths)
            incidence_angle = metadata['incidence_angle']
            exposure_time = metadata['exposure_time_s']
            timestamp = timestamp
            N = 1
        
        else:
            # multiple images
            metadata_list = [cls._get_NSLS_II_CMS_parameters(file, verbose=verbose) for file in filepaths]
            metadata_df = pd.DataFrame(metadata_list)

            if metadata_df['sample'].nunique() > 1:
                raise ValueError('Not all files have the same sample. Files cannot be averaged.')
            else:
                sample = metadata_df['sample'].unique()[0]

            if metadata_df['incidence_angle'].nunique() > 1:
                raise ValueError('Not all files have the same incidence angle. Files cannot be averaged.')
            else:
                incidence_angle = metadata_df['incidence_angle'].unique()[0]

            if metadata_df['exposure_time_s'].nunique() > 1:
                raise ValueError('Not all files have the same exposure time. Files cannot be averaged.')
            else:
                exposure_time = metadata_df['exposure_time_s'].unique()[0]

            if metadata_df['x_position'].nunique() > 1:
                raise ValueError('Not all files have the same x position. Files cannot be averaged.')
            
            if 'pos' in metadata_df.columns:
                if (len(metadata_df) == 2) and (metadata_df['pos'].nunique() == 2) and (1 in metadata_df['pos'].values) and (2 in metadata_df['pos'].values):
                    filename1 = metadata_df[metadata_df['pos'] == 1]['filepath'].values[0]
                    filename2 = metadata_df[metadata_df['pos'] == 2]['filepath'].values[0]
                    image = cls._stitch_images(filename1, filename2, offset = stiching_offset)
                    N = 1
                    metadata = {'sample': sample,
                               'filepaths': [filename1, filename2],
                               'relative humidity': metadata_df['relative humidity'].values.mean(),
                               'x_position': metadata_df['x_position'].values.mean()
                    }
                    
                else:
                    raise ValueError(f"""
                                     It seems like you need to stitch the files, but they are not compatible. \n
                                     The length of your metadata is {len(metadata_df)} and it should be 2 \n
                                     You may be trying to load more than two files to stitch together \n
                                     The files you are trying to stitch are {metadata_df['filepath'].values()} \n
                                     Each image needs to have 'pos1' and 'pos2' in the file name \n
                                     """)
                
            else:
                images_list = [cls._load_tif_file(f) for f in filepaths]
                images_array = np.array(images_list)
                image = np.squeeze(np.mean(images_array, axis=0))
                N = len(images_list)

                metadata = {'sample': sample,
                            'filepaths': filepaths,
                            'relative humidity': metadata_df['relative humidity'].values.mean(),
                            'x_position': metadata_df['x_position'].values,
                            }

        metadata['source'] = 'NSLS_II_CMS'

        return cls(image,
                   incidence_angle,
                   exposure_time,
                   timestamp,
                   metadata = metadata,
                   number_of_averaged_images = N)

    @staticmethod
    def _get_NSLS_II_CMS_parameters(filepath: str, verbose: bool = False) -> dict:
        """Get the parameters from the NSLS-II CMS beamline
        :param filepath: path to the image file
        :param verbose: whether to print the output
        :return: a dictionary with the parameters
        """
        filename = os.path.basename(filepath)
        parameters_from_file_name = {}
        parameters_from_file_name['filepath'] = filepath

        #Extract time duration
        time_matches = re.findall(r"\d+\.\d+s", filename)
        time_duration = time_matches[-1] if time_matches else None
        if time_duration is not None:
            time_duration = float(time_duration[:-1])
            parameters_from_file_name['exposure_time_s'] = time_duration

        # Extract angle (th)
        th_match = re.search(r"th(\d+\.\d+)", filename)
        th_angle = float(th_match.group(1)) if th_match else None
        if th_angle is not None:
            parameters_from_file_name['incidence_angle'] = th_angle
        else:
            raise ValueError('Incidence angle not found in the file name')

        # Extract x position
        x_match = re.search(r"x(-?\d+\.\d+)", filename)
        x_position = float(x_match.group(1)) if x_match else None
        if x_position is not None:
            parameters_from_file_name['x_position'] = x_position

        # Extract pos
        pos_match = re.search(r"_pos(\d+)", filename)
        pos = int(pos_match.group(1)) if pos_match else None
        if pos is not None:
            parameters_from_file_name['pos'] = pos

        # Extract RH
        rh_match = re.search(r"RH(\d+\.\d+)", filename)
        rh = float(rh_match.group(1)) if rh_match else None
        if rh is not None:
            parameters_from_file_name['relative humidity'] = rh

        # Extract series number
        series_match = re.search(r"_series_(\d+)_", filename)
        series = int(series_match.group(1)) if series_match else None
        if series is not None:
            # find progressive number
            prog_match = re.search(r"_(\d+)_waxs", filename)
            prog = int(prog_match.group(1)) if prog_match else None
            parameters_from_file_name['series'] = series
            parameters_from_file_name['progressive'] = prog
        else:
            prog = None

        # Extract sample name. This is the first part of the file name before am underscore which is followed by a number, "series", "pos", or "waxs"
        sample_match = re.search(r"^(.*?)_(?=\d+|series|pos|waxs)", filename)
        sample = sample_match.group(1) if sample_match else None
        if sample is not None:
            parameters_from_file_name['sample'] = sample

        if verbose:
            print(f"""
            Sample: {sample} \n
            Time Duration: {time_duration} \n
            Angle (th): {th_angle} \n
            X Position: {x_position} \n
            Pos: {pos} \n
            RH: {rh} \n
            Series Number: {series} \n
            Progressive Number: {prog} 
            """)
        
        return parameters_from_file_name

    @staticmethod
    def _stitch_images(file1: str,
                       file2:str,
                       offset:int = 30):
        """Merge two images with an offset # pixels in the y direction
        :param file1: path to the first image file
        :param file2: path to the second image file
        :param offset: offset in pixels in the y direction
        :return: the merged image as a NumPy array
        """
      
        array1 = GIWAXSPixelImage._load_tif_file(file1)
        array2 = GIWAXSPixelImage._load_tif_file(file2)
        merged_image = np.zeros_like(array1)

        # Loop through each row index i in array1
        for i in range(array1.shape[0]-offset):
            # Calculate the corresponding index in array2 (i + 10)
            j = i + offset

            for k in range(array1[i].size):
                if array1[i][k] == -1:
                    merged_image[i][k] = array2[j][k]
                elif array2[j][k] == -1:
                    merged_image[i][k] = array1[i][k]
                else:
                    merged_image[i][k] = (array1[i][k] + array2[j][k])/2

        return merged_image

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
        if importlib.util.find_spec('pygix') is None:
            raise ImportError('pygix is required to run this function. Please install pygix using pip install pygix')
        else:
            import pygix

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
    
    # Nick fix this
    def show(self, 
             engine:str = 'px', 
             **kwargs):
        """Plot the image.
        :param engine: The engine to use for plotting. Either plotly or hvplot.
        :return: The plot.
        """
        if engine == 'px':
            return self._show_px(**kwargs)

        elif engine == 'hv':
            return self._show_hv(**kwargs)

        else:
            raise ValueError('engine must be either px or hv')
        
    def _show_px(self, **kwargs):
        """Plot the image using plotly express.
        :param kwargs: additional arguments to pass to the plot
        :return: The plot.
        """
        fig = px.imshow(self._image, **kwargs)
        return fig
    
    def _show_hv(self, **kwargs):
        """Plot the image using holoviews.
        :param kwargs: additional arguments to pass to the plot
        :return: The plot.
        """
        if importlib.util.find_spec('holoviews') is None:
            raise ImportError('holoviews is required to run this function. Please install holoviews using pip install holoviews')
        else:
            import holoviews as hv
            hv.extension('bokeh')
        
        img = hv.Image(self._image, kdims=['x', 'y']).opts(**kwargs)
        return img
    
        
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
    def metadata(self):
        return self._metadata
    
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
                            engine:str = 'px',
                            **kwargs):
        """Plot the reciprocal space map.
        :param engine: The engine to use for plotting. Either plotly or hvplot.
        :return: The plot.
        """
        if engine == 'px':
            return self._plot_reciprocal_map_px(**kwargs)

        elif engine == 'hv':
            return self._plot_reciprocal_map_hv(**kwargs)

        else:
            raise ValueError('engine must be either px or hv')

    def _plot_reciprocal_map_px(self, 
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
        fig = self.plot_pixel_map_px(data = data, x='qxy', y='qz', z='intensity', colorscale=colorscale, log_scale=log_scale, z_lower_cuttoff=intensity_lower_cuttoff,
                            x_label='qxy [\u212B\u207B\u00B9]', y_label='qz [\u212B\u207B\u00B9]', template=template, origin=origin,
                            z_label='Intensity', **kwargs)
        return fig
    
    def _plot_reciprocal_map_hv(self, 
                            **kwargs) -> go.Figure:
       """Plot the reciprocal space map using holoviews
       :param kwargs: additional arguments to pass to the plot
       :return: The plot.
       """

       data = self.data_reciprocal.copy().sort_values(by=['qxy', 'qz'], ascending=[True, False])
       figure = self.plot_pixel_map_hv(data = data, x='qxy', y='qz', z='intensity',
                                     xlabel='qxy [\u212B\u207B\u00B9]',
                                     ylabel='qz [\u212B\u207B\u00B9]',
                                     clabel='Intensity [arb. units]', **kwargs)
       return figure
    

    def _plot_polar_map_contour_px(self, 
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
                                    z_lower_cuttoff=intensity_lower_cuttoff, template=template, x_label='q [\u212B\u207B\u00B9]', y_label='\u03C7 [\u00B0]', 
                                    z_label='Intensity', **kwargs)
        return fig
    
    def plot_polar_map(self, engine:str = 'px', **kwargs):
        """Plot the polar space map.
        :param engine: The engine to use for plotting. Either plotly or hvplot.
        :return: The plot.
        """
        if engine == 'px':
            return self._plot_polar_map_px(**kwargs)
        elif engine == 'hv':
            return self._plot_polar_map_hv(**kwargs)
        else:
            raise ValueError('engine must be either px or hv')
    
    def _plot_polar_map_px(self, 
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

        fig = self.plot_pixel_map_px(data = data, y='chi', x='q', z='intensity', colorscale=colorscale, aspect='auto', z_lower_cuttoff=intensity_lower_cuttoff,
                                     origin=origin, log_scale=log_scale,x_label='Q [\u212B\u207B\u00B9]', y_label='\u03C7 [\u00B0]', 
                                     z_label='Intensity', template=template, **kwargs)
        return fig
    
    def _plot_polar_map_hv(self, 
                       **kwargs):
        """Plot the polar space map.
        :param kwargs: additional arguments to pass to the plot
        :return: The plot.
        """
        data = self.data_polar.copy().sort_values(by=['q', 'chi'], ascending=[True, False])
        figure = self.plot_pixel_map_hv(data = data, y='chi', x='q', z='intensity',
                                     xlabel='Q [\u212B\u207B\u00B9]', ylabel='\u03C7 [\u00B0]',
                                     clabel='Intensity [arb. units]', **kwargs)
        return figure
    
    def get_linecut(self,
                    chi : tuple | list | pd.Series | float = None,
                    q_range : tuple | list | pd.Series = None) -> 'Linecut':
        """Extract a profile from the polar space data.
        :param chi: Range of chi values or a single chi value.
        :param q_range: q_range.
        :return: Lincut object.
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

        metadata = self.metadata.copy()
        metadata['chi'] = chi
        metadata['q_range'] = q_range
        
        return Linecut(data, metadata = metadata)
      
    def get_polar_linecut(self,
                    q : tuple | list | pd.Series | float = None,
                    chi_range : tuple | list | pd.Series = None) -> 'Polar_linecut':
        """Extract a profile from the polar space data.
        :param q: Range of q values or a single q value.
        :param chi_range: chi_range.
        :return: Polar_linecut object.
        """

        data = self.data_polar.copy()

        # check if q is iterable
        try:
            iter(q)
            q_iterable = True
            if len(q) != 2: raise ValueError('If q is a range it must be two values')
        except TypeError:
            q_iterable = False
            if q < data['q'].min() or q > data['q'].max(): raise ValueError('q value out of range of the data')

        # Filter the data for q
        if q_iterable: 
            data = data.query(f'q >= {min(q)} and q <= {max(q)}')
        else:
            closest_index = data['q'].sub(q).abs().idxmin()
            closest_q = data.loc[closest_index, 'q']
            data = data.query(f'q == {closest_q}')

        # Filter the data for chi
        if chi_range is not None: 
            data = data.query(f'chi >= {min(chi_range)} and chi <= {max(chi_range)}')
        
        data = data.groupby('chi').mean().reset_index().filter(['q', 'chi', 'intensity'])

        metadata = self.metadata.copy()
        metadata['chi_range'] = chi_range
        metadata['q'] = q
        
        return Polar_linecut(data, metadata = metadata)
         
    

class Linecut():
    ''' 
    A class to store a linecut from a GIWAXS measurement

    Main contributors:
    Arianna Magni
    '''

    def __init__(self,
                 data: pd.DataFrame,
                 metadata: dict = None):
        
        self._data = data
        self._metadata = metadata

    @property
    def data(self):
        return self._data
    
    @property
    def metadata(self):
        return self._metadata
    
    @property
    def chi(self):
        try:
            return self.metadata['chi']
        except:
            raise AttributeError('No chi value has been set.')
    
    @property
    def fit_results(self):
        if hasattr(self, '_fit_results'):
            return self._fit_results
        else:
            raise AttributeError('No fit has been ran.')
    
    @property
    def fit_model(self):
        try:
            return self.fit_results.model
        except:
            raise AttributeError('No fit model has been set.')
    
    @property
    def fit_params(self):
        try:
            return self.fit_results.params
        except:
            raise AttributeError('No fit has been ran.')
    
    @property
    def fit_report(self):
        try:
            return self.fit_results.fit_report
        except:
            raise AttributeError('No fit has been ran.')
    
    @property
    def x(self):
        return self._x
    
    @property
    def y(self):
        return self._y
    
    @property
    def y_fit(self):
        try:
            return self.fit_results.best_fit
        except:
            raise AttributeError('No fit has been ran.')
    
    def subtract_background(self,
                            background,
                            background_metadata: dict = {}) -> 'Linecut':
        """Subtract the background from the linecut.
        :param background: The background data as a pandas DataFrame or Linecut object
        :param background_metadata: The metadata for the background data.
        :return: The linecut with the background subtracted.
        """
        if isinstance(background, Linecut):
            background_metadata.update(background.metadata)
            background_df = background.data
        elif isinstance(background, pd.DataFrame):
            background_df = background.copy()
        else:
            raise ValueError('background must be a pandas DataFrame or Linecut object')
        
        data = self.data.copy()

        #check background has q and intensity columns
        if 'q' not in background_df.columns or 'intensity' not in background_df.columns:
            raise ValueError('background must have q and intensity columns')
        
        # make sure the background has the same q values as the linecut, if not interpolate values
        # Nick look into this 
        if not set(self.data['q']).issubset(set(background_df['q'])):
            background_df = background_df.set_index('q').reindex(data['q']).reset_index().interpolate()
               
        background_df = background_df.rename(columns={'intensity': 'background_intensity'})
        data = data.merge(background_df, on='q', how='left')
        data['intensity_raw'] = data['intensity']
        data['intensity'] = data['intensity'] - data['background_intensity']
        self._data = data
        self._metadata['background_metadata'] = background_metadata

        return self

    def plot(self,
             engine: str = 'px',
             show_background: bool = False,
                **kwargs) -> px.line:
        """Plot the profile.
        :param engine: The engine to use for plotting. Either plotly or hvplot.
        :param show_background: Whether to show the background.
        :return: The plot.
        """
        if engine == 'px':
            return self._plot_px(show_background, **kwargs)
        elif engine == 'hv':
            return self._plot_hv(show_background, **kwargs)
        else:
            raise ValueError('engine must be either px or hv')
        
    def _plot_hv(self,
                 show_background: bool = False,
                label: str = '',
                **kwargs):
        """Plot the profile.
        :param show_background: Whether to show the background.
        :param label: The label for the plot.
        :param kwargs: additional arguments to pass to the plot
        :return: The hv plot.
        """
        if importlib.util.find_spec('holoviews') is None:
            raise ImportError('holoviews is required to run this function. Please install holoviews using pip install holoviews')
        else:
            import holoviews as hv
            hv.extension('bokeh')

        if show_background:
            if 'background_intensity' not in self.data.columns:
                raise ValueError('No background has been subtracted.')
            background_curve = hv.Curve(self.data, kdims='q', vdims='background_intensity', label = label + ' background').opts(
                xlabel='q [\u212B\u207B\u00B9]',
                ylabel='Intensity [arb. units]',
                color = 'black',
                **kwargs)
            
            row_curve = hv.Curve(self.data, kdims='q', vdims='intensity_raw', label = label + ' raw').opts(
                xlabel='q [\u212B\u207B\u00B9]',
                ylabel='Intensity [arb. units]',
                color = 'blue',
                **kwargs)
            
            curve = hv.Curve(self.data, kdims='q', vdims='intensity', label = label + ' after background subtraction').opts(
                xlabel='q [\u212B\u207B\u00B9]',
                ylabel='Intensity [arb. units]',
                color = 'red',
                **kwargs)
            
            return hv.Overlay([row_curve, background_curve, curve])

        else:
            curve = hv.Curve(self.data, kdims='q', vdims='intensity', label = label).opts(
                xlabel='q [\u212B\u207B\u00B9]',
                ylabel='Intensity [arb. units]',
            **kwargs)
            return curve
    
    def _plot_px(self, 
                    show_background: bool = False,
                    **kwargs) -> px.line:
        """Plot the profile.
        :param show_background: Whether to show the background.
        :param kwargs: additional arguments to pass to the plot
        :return: The plot.
        """
        profile = self.data
        
        if show_background:
            if 'background_intensity' not in profile.columns:
                raise ValueError('No background has been subtracted.')
            profile_melted = profile.melt(id_vars="q", var_name="Linecut", value_name="Intensity")
            figure = px.line(profile_melted, x="q", y="Intensity", color="Linecut", labels={'q': 'q [\u212B\u207B\u00B9]', 'intensity': 'Intensity [arb. units]'}, **kwargs)
        else:
            figure = px.line(profile, x='q', y='intensity', labels={'q': 'q [\u212B\u207B\u00B9]', 'intensity': 'Intensity [arb. units]'}, **kwargs)
        
        return figure
    
    def remove_spikes(self,
                       q_range: tuple = None,
                       threshold: float = None,
                       window: int = 3) -> 'Linecut':
        """Remove cosmic rays from the linecut.
        :param q_range: The range of q values to consider.
        :param threshold: The threshold for the z-score.
        :param window: The window size for the rolling mean.
        """
        data = self.data.copy()
        if q_range is None:
            q_range = (data['q'].min(), data['q'].max())
        data_q_range = data.query(f'q >= {q_range[0]} and q <= {q_range[1]}')
        
        data_q_range.loc[:, 'intensity_mean'] = data_q_range['intensity'].rolling(window, center=True).mean()
        data_q_range.loc[:, 'z_score'] = (data_q_range['intensity'] - data_q_range['intensity_mean']) 
    
        # Mask for keeping only non-spikes (z-score within threshold)
        mask = data_q_range['z_score'].abs() <= threshold
        data_filtered__q_range = data_q_range[mask]

        #remove newly created colums
        data_filtered__q_range = data_filtered__q_range.drop(columns=['intensity_mean', 'z_score'])

        # Update the data, keep everything outside the q_range the same
        data_ext_q_range = data[~data['q'].isin(data_q_range['q'])]
        data_filtered = pd.concat([data_ext_q_range, data_filtered__q_range]).sort_values(by='q')
        self._data = data_filtered    

        return self
    
    def fit_linecut(self,
                    peak_model: str,
                    background_model: str,
                    q_range: tuple,
                    initial_parameters: dict = {}) -> 'Linecut':          
        """Fit the linecut to a model
        :param peak_model: The peak model to use. Options are 'GaussianModel', 'LorentzianModel', 'VoigtModel', 'PseudoVoigtModel', 'SkewedVoigtModel'
        :param background_model: The background model to use. Options are 'ExponentialModel', 'LinearModel', 'ConstantModel', PowerLawModel
        :param q_range: The range of q values to fit
        :param initial_parameters: The initial parameters for the fit
        :return: The fit results.
        """        
        default_fit_parameters = {
            'peak_center_value': 1.0,
            'peak_center_vary': True,
            'peak_center_min': 0.1,
            'peak_center_max': 2.5,

            'peak_sigma_value': 0.1,
            'peak_sigma_vary': True,
            'peak_sigma_min': 0.001,
            'peak_sigma_max': 0.3,

            'peak_amplitude_value': 1,
            'peak_amplitude_vary': True,
            'peak_amplitude_min': 0.00001,
            'peak_amplitude_max': 5000,

            'peak_gamma_value': 0.1,
            'peak_gamma_vary': True,
            'peak_gamma_min': 0.001,
            'peak_gamma_max': 0.3,

            'peak_fraction_value': 0.5,
            'peak_fraction_vary': True,
            'peak_fraction_min': 0.000,
            'peak_fraction_max': 1.0,

            'peak_skew_value': 0,
            'peak_skew_vary': True,
            'peak_skew_min': -1000,
            'peak_skew_max': 1000,
            
            'bkg_slope_value': 0,
            'bkg_slope_vary': True,
            'bkg_slope_min': -1000,
            'bkg_slope_max': 1000,

            'bkg_intercept_value': 0,
            'bkg_intercept_vary': True,
            'bkg_intercept_min': -1000,
            'bkg_intercept_max': 1000,

            'bkg_value': 0,
            'bkg_vary': True,
            'bkg_min': 0,
            'bkg_max': 1000,

            'bkg_amplitude_value': 0,
            'bkg_amplitude_vary': True,
            'bkg_amplitude_min': -1000,
            'bkg_amplitude_max': 1000,

            'bkg_decay_value': 1,
            'bkg_decay_vary': True,
            'bkg_decay_min': -1000,
            'bkg_decay_max': 1000,

            'bkg_exponent_value': 1,
            'bkg_exponent_vary': True,
            'bkg_exponent_min': -1000,
            'bkg_exponent_max': 1000
        }

        for key in initial_parameters.keys():
            if key not in default_fit_parameters.keys():
                raise ValueError(f'{key} is not a valid parameter. Available parameters are {default_fit_parameters.keys()}')
   
        default_fit_parameters.update(initial_parameters)
                
        from lmfit.models import (ExponentialModel,
                                  PseudoVoigtModel,
                                  SkewedVoigtModel,
                                  VoigtModel,
                                  LinearModel,
                                  ConstantModel,
                                  GaussianModel,
                                  LorentzianModel,
                                  PowerLawModel)
        
        # select values in q_range
        data = self.data.query(f'q >= {q_range[0]} and q <= {q_range[1]}')
        x = data['q']
        y = data['intensity']
        
        peak_model_dict = {
            'GaussianModel': GaussianModel(prefix='peak_'),
            'LorentzianModel': LorentzianModel(prefix='peak_'),
            'VoigtModel': VoigtModel(prefix='peak_'),
            'PseudoVoigtModel': PseudoVoigtModel(prefix='peak_'),
            'SkewedVoigtModel': SkewedVoigtModel(prefix='peak_')
        }

        if peak_model not in peak_model_dict:
            raise ValueError('peak_model must be one of GaussianModel, LorentzianModel, VoigtModel, PseudoVoigtModel, SkewedVoigtModel')

        selected_peak_model = peak_model_dict[peak_model]
        
        background_model_dict = {
            'ExponentialModel': ExponentialModel(prefix='bkg_') + ConstantModel(prefix='bkg_'),
            'LinearModel': LinearModel(prefix='bkg_'),
            'ConstantModel': ConstantModel(prefix='bkg_'),
            'PowerLawModel': PowerLawModel(prefix='bkg_') + ConstantModel(prefix='bkg_')
        }

        if background_model not in background_model_dict:
            raise ValueError('background_model must be one of ExponentialModel, LinearModel, ConstantModel, PowerLawModel')

        selected_background_model = background_model_dict[background_model]
        
        model = selected_peak_model + selected_background_model
        pars = model.make_params()
        
        pars['peak_center'].set(value=default_fit_parameters['peak_center_value'],
                                min=default_fit_parameters['peak_center_min'],
                                max=default_fit_parameters['peak_center_max'],
                                vary=default_fit_parameters['peak_center_vary'])
                                       
        pars['peak_sigma'].set(value=default_fit_parameters['peak_sigma_value'],
                               min=default_fit_parameters['peak_sigma_min'],
                               max=default_fit_parameters['peak_sigma_max'],
                               vary=default_fit_parameters['peak_sigma_vary'])
    
        pars['peak_amplitude'].set(value=default_fit_parameters['peak_amplitude_value'],
                                   min=default_fit_parameters['peak_amplitude_min'],
                                   max=default_fit_parameters['peak_amplitude_max'],
                                   vary=default_fit_parameters['peak_amplitude_vary'])

        if peak_model == 'PseudoVoigtModel':
            pars['peak_fraction'].set(value=default_fit_parameters['peak_fraction_value'],
                                   min=default_fit_parameters['peak_fraction_min'],
                                   max=default_fit_parameters['peak_fraction_max'],
                                   vary=default_fit_parameters['peak_fraction_vary'])
                                   
        elif (peak_model == 'VoigtModel') or (peak_model == 'SkewedVoigtModel'):
            pars['peak_gamma'].set(value=default_fit_parameters['peak_gamma_value'],
                                   min=default_fit_parameters['peak_gamma_min'],
                                   max=default_fit_parameters['peak_gamma_max'],
                                   vary=default_fit_parameters['peak_gamma_vary'])

            if peak_model == 'SkewedVoigtModel':
                pars['peak_skew'].set(value=default_fit_parameters['peak_skew_value'],
                                   min=default_fit_parameters['peak_skew_min'],
                                   max=default_fit_parameters['peak_skew_max'],
                                   vary=default_fit_parameters['peak_skew_vary'])

        pars.add('peak_d_spacing', expr='2*pi/peak_center')
        if peak_model != 'SkewedVoigtModel':
            pars.add('peak_coherence_length', expr='2*pi*0.9/peak_fwhm')
                 
        if background_model == 'ExponentialModel':
            pars['bkg_amplitude'].set(value=default_fit_parameters['bkg_amplitude_value'],
                                      min=default_fit_parameters['bkg_amplitude_min'],
                                      max=default_fit_parameters['bkg_amplitude_max'],
                                      vary=default_fit_parameters['bkg_amplitude_vary'])    

            pars['bkg_decay'].set(value=default_fit_parameters['bkg_decay_value'],
                                    min=default_fit_parameters['bkg_decay_min'],
                                    max=default_fit_parameters['bkg_decay_max'],
                                    vary=default_fit_parameters['bkg_decay_vary'])

            pars['bkg_'].set(value=default_fit_parameters['bkg_value'],
                             min=default_fit_parameters['bkg_min'],
                             max=default_fit_parameters['bkg_max'],
                             vary=default_fit_parameters['bkg_vary'])
        
        elif background_model == 'PowerLawModel':
            pars['bkg_amplitude'].set(value=default_fit_parameters['bkg_amplitude_value'],
                                      min=default_fit_parameters['bkg_amplitude_min'],
                                      max=default_fit_parameters['bkg_amplitude_max'],
                                      vary=default_fit_parameters['bkg_amplitude_vary'])

            pars['bkg_exponent'].set(value=default_fit_parameters['bkg_exponent_value'],
                                    min=default_fit_parameters['bkg_exponent_min'],
                                    max=default_fit_parameters['bkg_exponent_max'],
                                    vary=default_fit_parameters['bkg_exponent_vary'])

            pars['bkg_'].set(value=default_fit_parameters['bkg_value'],
                             min=default_fit_parameters['bkg_min'],
                             max=default_fit_parameters['bkg_max'],
                             vary=default_fit_parameters['bkg_vary'])
                        
        elif background_model == 'LinearModel':
            pars['bkg_slope'].set(value=default_fit_parameters['bkg_slope_value'],
                                  min=default_fit_parameters['bkg_slope_min'],
                                  max=default_fit_parameters['bkg_slope_max'],
                                  vary=default_fit_parameters['bkg_slope_vary'])

            pars['bkg_intercept'].set(value=default_fit_parameters['bkg_intercept_value'],
                                      min=default_fit_parameters['bkg_intercept_min'],
                                      max=default_fit_parameters['bkg_intercept_max'],
                                      vary=default_fit_parameters['bkg_intercept_vary'])
        
        elif background_model == 'ConstantModel':
            pars['bkg_'].set(value=default_fit_parameters['bkg_value'],
                             min=default_fit_parameters['bkg_min'],
                             max=default_fit_parameters['bkg_max'],
                             vary=default_fit_parameters['bkg_vary'])
            
        result = model.fit(y, pars, x=x)

        if peak_model == 'SkewedVoigtModel':   
            fitted_y = result.best_fit
            half_max = max(fitted_y) / 2
            indices = np.where(fitted_y >= half_max)[0]
            fwhm = x.to_numpy()[indices[-1]] - x.to_numpy()[indices[0]]
            # Add calculated FWHM to the parameters
            result.params.add('peak_fwhm', value=fwhm, vary=False)
            result.params.add('peak_coherence_length', expr='2*pi*0.9/peak_fwhm')

        self._x = x
        self._y = y
        self._fit_results = result

        
        return self

    def plot_fitted(self,
                    engine: str = 'px',
                    **kwargs):
        """Plot the fitted linecut
        :param engine: The engine to use for plotting. Either plotly or hvplot.
        :return: The plot.
        """
        if engine == 'px':
            return self._plot_fitted_px(**kwargs)
        elif engine == 'hv':
            return self._plot_fitted_hv(**kwargs)
        else:
            raise ValueError('engine must be either px or hv')
            
    def _plot_fitted_hv(self, **kwargs):
        """Plot the fitted linecut using holoviews
        :param kwargs: additional arguments to pass to the plot
        :return: The plot.
        """
        if importlib.util.find_spec('holoviews') is None:
            raise ImportError('holoviews is required to run this function. Please install holoviews using pip install holoviews')
        else:
            import holoviews as hv
            hv.extension('bokeh')

        curve_data = self._plot_hv(label = 'Data').opts(
            line_width=3,
            color='black',
            **kwargs)
        
        curve_fit = hv.Curve((self.x, self.y_fit), kdims='q', vdims='intensity', label = 'Fit').opts(
            xlabel='q [\u212B\u207B\u00B9]',
            ylabel='Intensity [arb. units]',
            color='red',
            line_width=2,
            **kwargs)
        
        curve_peak = hv.Curve((self.x, self.fit_results.eval_components()['peak_']), kdims='q', vdims='intensity', label = 'Peak').opts(
            xlabel='q [\u212B\u207B\u00B9]',
            ylabel='Intensity [arb. units]',
            line_dash='dashed',
            color='green',
            line_width=2,
            **kwargs)
        
        curve_bkg = hv.Curve((self.x, self.fit_results.eval_components()['bkg_']), kdims='q', vdims='intensity', label = 'Background').opts(
            xlabel='q [\u212B\u207B\u00B9]',
            ylabel='Intensity [arb. units]',
            color='blue',
            line_dash='dashed',
            line_width=2,
            **kwargs)
        
        return hv.Overlay([curve_data, curve_fit, curve_peak, curve_bkg])
    
    def _plot_fitted_px(self, **kwargs) -> px.line:
        """Plot the fitted linecut using plotly express
        :param kwargs: additional arguments to pass to the plot
        :return: The plot.
        """
        data = self.data
        data['fitted'] = self.fit_results.best_fit
        data['peak'] = self.fit_results.eval_components()['peak_']
        data['bkg'] = self.fit_results.eval_components()['bkg_']
        figure = px.line(data, x='q', y=['intensity', 'fitted', 'peak', 'bkg'], labels={'value': 'Intensity [a.u.]', 'variable': 'Fit components'})
        return figure
    

class Polar_linecut():
    ''' 
    A class to store a polar linecut from a GIWAXS measurement

    Main contributors:
    Arianna Magni

    '''

    def __init__(self,
                 data: pd.DataFrame,
                 metadata: dict = None):
        
        self._data = data
        self._metadata = metadata

    @property
    def data(self):
        return self._data
    
    @property
    def metadata(self):
        return self._metadata
    
    @property
    def q(self):
        try:
            return self.metadata['q']
        except:
            raise AttributeError('No q value has been set.')
    
    @property
    def chi(self):
        return self.data['chi']

    @property
    def intensity(self):
        return self.data['intensity']   

    def plot(self,
             engine: str = 'px',
                **kwargs) -> px.line:
        """Plot the profile.
        :param engine: The engine to use for plotting. Either plotly or hvplot.
        :return: The plot.
        """
        if engine == 'px':
            return self._plot_px(**kwargs)
        elif engine == 'hv':
            return self._plot_hv(**kwargs)
        else:
            raise ValueError('engine must be either px or hv')
        
    def _plot_hv(self,
                label: str = '',
                **kwargs):
        """Plot the profile.
        :param label: The label for the plot.
        :param kwargs: additional arguments to pass to the plot
        :return: The hv plot.
        """
        if importlib.util.find_spec('holoviews') is None:
            raise ImportError('holoviews is required to run this function. Please install holoviews using pip install holoviews')
        else:
            import holoviews as hv
            hv.extension('bokeh')
        
        curve = hv.Curve(self.data, kdims='chi', vdims='intensity', label = label).opts(
            xlabel='\u03C7 [\u00B0]',
            ylabel='Intensity [arb. units]',
            **kwargs)
        
        return curve
    
    def _plot_px(self, 
             **kwargs) -> px.line:
        """Plot the profile.
        :param kwargs: additional arguments to pass to the plot
        :return: The plot.
        """
        profile = self.data
        figure = px.line(profile, x='q', y='intensity', labels={'chi': '\u03C7 [\u00B0]', 'intensity': 'Intensity'}, **kwargs)
        return figure