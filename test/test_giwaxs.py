import unittest
import tracemalloc
import pandas as pd
import numpy as np
import plotly.graph_objects as go
from Materials_Data_Analytics.experiment_modelling.giwaxs import Calibrator
from Materials_Data_Analytics.experiment_modelling.giwaxs import GIWAXSPixelImage, GIWAXSPattern, Linecut, Polar_linecut
import plotly.express as px
import plotly as pl
import holoviews as hv
hv.extension('bokeh')


class TestCalibration(unittest.TestCase):
    ''' Test the Calibrator class '''
    def setUp(self):
        self.my_calibrator = Calibrator.from_poni_file('./test_trajectories/giwaxs/calibration_SLAC_BL11_3.poni')

    def test_attributes(self):
        ''' Test the attributes of the Calibration class '''
        self.assertEqual(self.my_calibrator.wavelength, 0.09919)
        self.assertEqual(self.my_calibrator.distance, 0.28556)
        self.assertEqual(self.my_calibrator.poni1, 0.21226)
        self.assertEqual(self.my_calibrator.poni2, 0.11560)
        self.assertEqual(self.my_calibrator.rot1, 0.0025858)
        self.assertEqual(self.my_calibrator.rot2, 0.0093694)
        self.assertEqual(self.my_calibrator.rot3, 0.0000000)


class TestGIWAXSPixelImage(unittest.TestCase):
    ''' Test the GIWAXS class '''
    def setUp(self): 
        self.data_SLAC_BL113 = GIWAXSPixelImage.from_SLAC_BL11_3(tif_filepaths= ['./test_trajectories/giwaxs/GIWAXS_image_SLAC_1.tif'])
        self.ai_SLAC_BL113 = Calibrator.from_poni_file('./test_trajectories/giwaxs/calibration_SLAC_BL11_3.poni')
        self.ai_NSLS_II_CMS = Calibrator.from_poni_file('./test_trajectories/giwaxs/calibration_NSLS_II_CMS.poni')

        self.data_NSLS_II_CMS_stiched = GIWAXSPixelImage.from_NSLS_II_CMS(filepaths= [
            r'./test_trajectories/giwaxs/GIWAXS_image_NSLS_II_CMS_pos1_31_1563.0s_RH1.396_x-1.500_th0.100_10.00s_1711351_waxs.tiff',
            r'./test_trajectories/giwaxs/GIWAXS_image_NSLS_II_CMS_pos2_36_1641.3s_RH1.028_x-1.500_th0.100_10.00s_1711356_waxs.tiff'])

    def test_from_SLAC_BL113(self):
        ''' Test the attributes of the GIWAXS class '''
        self.assertTrue(self.data_SLAC_BL113.image.shape == (3072, 3072))
        self.assertTrue(np.round(self.data_SLAC_BL113.image[5][2], 3) == 24.229)
        self.assertTrue(np.round(self.data_SLAC_BL113.image[15][27], 3) == 23.279)
        self.assertTrue(np.round(self.data_SLAC_BL113.image[257][43], 3) == 28.505)
        self.assertTrue(self.data_SLAC_BL113.incidence_angle == 0.12)
        self.assertTrue(self.data_SLAC_BL113.exposure_time == 120.0)

        ''' Test the mask method of the GIWAXS class '''
        self.data_SLAC_BL113.apply_mask(mask_path='./test_trajectories/giwaxs/mask_SLAC.tif')
        self.assertTrue(np.isnan(self.data_SLAC_BL113.image[3000][43]))
        
        ''' Test the transform method of the GIWAXS class '''
        my_giwaxs_pattern = self.data_SLAC_BL113.get_giwaxs_pattern(calibrator = self.ai_SLAC_BL113,
                                                                    qxy_range = (-3, 3),
                                                                    qz_range = (0, 3),
                                                                    q_range = (0, 3),
                                                                    chi_range = (-95, 95),
                                                                    pixel_q = 100,
                                                                    pixel_chi = 80,
                                                                    correct_solid_angle = True,
                                                                    polarization_factor = None,
                                                                    unit = 'A')
        
        self.assertTrue(type(my_giwaxs_pattern == GIWAXSPattern))
        self.assertTrue('qxy' in my_giwaxs_pattern.data_reciprocal.columns)
        self.assertTrue('qz' in my_giwaxs_pattern.data_reciprocal.columns)
        self.assertTrue('intensity' in my_giwaxs_pattern.data_reciprocal.columns)
        self.assertTrue('chi' in my_giwaxs_pattern.data_polar.columns)
        self.assertTrue('q' in my_giwaxs_pattern.data_polar.columns)
        self.assertTrue('intensity' in my_giwaxs_pattern.data_polar.columns)
    
    def test_from_NSLS_II_CMS(self):
        ''' Test the attributes of the GIWAXS class '''
        self.assertTrue(self.data_NSLS_II_CMS_stiched.image.shape == (1043, 981))
        self.assertTrue(np.round(self.data_NSLS_II_CMS_stiched.image[5][2], 3) == 0)
        self.assertTrue(np.round(self.data_NSLS_II_CMS_stiched.image[15][27], 3) == 4)
        self.assertTrue(np.round(self.data_NSLS_II_CMS_stiched.image[257][43], 3) == 5)
        self.assertTrue(np.round(self.data_NSLS_II_CMS_stiched.image[257][500], 3) == 845)
        self.assertTrue(self.data_NSLS_II_CMS_stiched.incidence_angle == 0.1)
        self.assertTrue(self.data_NSLS_II_CMS_stiched.exposure_time == 10.0)

        ''' Test the mask method of the GIWAXS class '''
        self.data_NSLS_II_CMS_stiched.apply_mask(mask_path='./test_trajectories/giwaxs/mask_NSLS_II_CMS_stiched.tif')
        self.assertTrue(np.isnan(self.data_NSLS_II_CMS_stiched.image[257][43]))
        self.assertTrue(np.round(self.data_NSLS_II_CMS_stiched.image[257][500], 3) == 845)

        
        ''' Test the transform method of the GIWAXS class '''
        my_giwaxs_pattern = self.data_NSLS_II_CMS_stiched.get_giwaxs_pattern(calibrator = self.ai_NSLS_II_CMS,
                                                                    qxy_range = (-2, 2),
                                                                    qz_range = (0, 2),
                                                                    q_range = (0, 2),
                                                                    chi_range = (-95, 95),
                                                                    pixel_q = 100,
                                                                    pixel_chi = 80,
                                                                    correct_solid_angle = True,
                                                                    polarization_factor = None,
                                                                    unit = 'A')
        
        self.assertTrue(type(my_giwaxs_pattern == GIWAXSPattern))
        self.assertTrue('qxy' in my_giwaxs_pattern.data_reciprocal.columns)
        self.assertTrue('qz' in my_giwaxs_pattern.data_reciprocal.columns)
        self.assertTrue('intensity' in my_giwaxs_pattern.data_reciprocal.columns)
        self.assertTrue('chi' in my_giwaxs_pattern.data_polar.columns)
        self.assertTrue('q' in my_giwaxs_pattern.data_polar.columns)
        self.assertTrue('intensity' in my_giwaxs_pattern.data_polar.columns)


class TestGIWAXSPattern(unittest.TestCase):
    ''' Test the GIWAXSPattern class '''   
    def setUp(self):                        
        ai = Calibrator.from_poni_file('./test_trajectories/giwaxs/calibration_SLAC_BL11_3.poni')
        self.data_SLAC_BL113 = (GIWAXSPixelImage
                                .from_SLAC_BL11_3(tif_filepaths= ['./test_trajectories/giwaxs/GIWAXS_image_SLAC_1.tif'])
                                .get_giwaxs_pattern(calibrator = ai,
                                                    qxy_range = (-3, 3),
                                                    qz_range = (0, 3),
                                                    q_range = (0, 3),
                                                    chi_range = (-95, 95),
                                                    pixel_q = 100,
                                                    pixel_chi = 80,
                                                    correct_solid_angle = True,
                                                    polarization_factor = None,
                                                    unit = 'A'
                                                    )
                                )
    
    def test_attributes(self):
        ''' Test the attributes of the GIWAXSPattern class '''
        self.assertTrue(self.data_SLAC_BL113.qxy.shape == (80,))
        self.assertTrue(self.data_SLAC_BL113.qz.shape == (100,))
        self.assertTrue(self.data_SLAC_BL113.q.shape == (100,))
        self.assertTrue(self.data_SLAC_BL113.chi.shape == (151,))
        self.assertTrue(len(self.data_SLAC_BL113.data_polar) == 8179)
        self.assertTrue(len(self.data_SLAC_BL113.data_reciprocal) == 7206)

    def test_plotting(self):
        
        figure = self.data_SLAC_BL113.plot_reciprocal_map_contour(intensity_lower_cuttoff = 30)
        self.assertTrue(type(figure) == go.Figure)   

        figure = self.data_SLAC_BL113.plot_reciprocal_map(width=800, height=500, intensity_lower_cuttoff = 30)
        self.assertTrue(type(figure) == go.Figure)   

        figure = self.data_SLAC_BL113.plot_reciprocal_map(engine='hv')
        self.assertTrue(type(figure) == hv.Image)

        figure = self.data_SLAC_BL113._plot_polar_map_contour_px(width=800, height=500, intensity_lower_cuttoff = 30)
        self.assertTrue(type(figure) == go.Figure)     

        figure = self.data_SLAC_BL113.plot_polar_map(intensity_lower_cuttoff = 30)
        self.assertTrue(type(figure) == go.Figure)

        figure = self.data_SLAC_BL113.plot_polar_map(engine='hv')
        self.assertTrue(type(figure) == hv.Image)      

    def test_get_linecuts(self):
        ''' Test the Linecut classes ''' 
        linecut = self.data_SLAC_BL113.get_linecut(q_range = (0, 2.5), chi = (0,20))
        self.assertTrue(type(linecut) == Linecut)
        self.assertTrue(len(linecut.data) == 83)
        self.assertTrue('q' in linecut.data.columns)
        self.assertTrue('intensity' in linecut.data.columns)

        polar_linecut = self.data_SLAC_BL113.get_polar_linecut(q = (1), chi_range = (0,90))
        self.assertTrue(type(polar_linecut) == Polar_linecut)
        self.assertTrue(len(polar_linecut.data) == 37)
        self.assertTrue('chi' in polar_linecut.data.columns)
        self.assertTrue('intensity' in polar_linecut.data.columns)

class TestLinecut(unittest.TestCase):
    ''' Test the Linecut class '''
    def setUp(self):
        x = np.linspace(0, 2, 103)
        amplitude = 10
        center = 1
        sigma = 0.1
        y = (amplitude/(sigma*np.sqrt(2*np.pi))) * np.exp(-0.5 * (x - center) ** 2 / sigma ** 2)+1+0.5*x
        # add noise to the y vector
        y += np.random.normal(0, 0.1, y.shape)                        
        self.my_linecut = Linecut(data = pd.DataFrame({'q':x,'intensity': y}))

    def test_attributes(self):
        ''' Test the attributes of the Linecut class '''
        self.assertTrue(len(self.my_linecut.data) == 103)
        self.assertTrue('q' in self.my_linecut.data.columns)
        self.assertTrue('intensity' in self.my_linecut.data.columns)
    
    def test_plotting(self):
        ''' Test the plotting methods of the Linecut class '''
        figure = self.my_linecut.plot()
        self.assertTrue(type(figure) == go.Figure)

        figure = self.my_linecut.plot(engine='hv')
        self.assertTrue(type(figure) == hv.Curve)
    
    def test_fitting(self):
        from lmfit.model import ModelResult, Parameters 
        ''' Test the fitting methods of the Linecut class '''
        self.my_linecut = self.my_linecut.fit_linecut(peak_model = 'GaussianModel', background_model = 'LinearModel', q_range = (0.5, 1.5),
                                                 initial_parameters = {'peak_center_value': 0.9,
                                                                       'peak_amplitude_value': 9,
                                                                       'peak_sigma_value': 0.2,
                                                                       'bkg_slope_value': 0, 
                                                                       'bkg_intercept_value': 0})
        self.assertTrue(type(self.my_linecut.fit_results) == ModelResult)
        self.assertTrue(type(self.my_linecut.fit_params) == Parameters)
        self.assertTrue(type(self.my_linecut.fit_report) == str)
        self.assertTrue(len(self.my_linecut.y_fit) == 51)
        self.assertTrue(len(self.my_linecut.y) == 51)
        self.assertTrue(bool(self.my_linecut.fit_params['peak_center'].value > 0.95) and bool(self.my_linecut.fit_params['peak_center'].value < 1.05))
        self.assertTrue(bool(self.my_linecut.fit_params['peak_amplitude'].value > 9) and bool(self.my_linecut.fit_params['peak_amplitude'].value < 11))
        self.assertTrue(bool(self.my_linecut.fit_params['peak_sigma'].value > 0.085) and bool(self.my_linecut.fit_params['peak_sigma'].value < 0.115))
        self.assertTrue(bool(self.my_linecut.fit_params['bkg_slope'].value > 0.4) and bool(self.my_linecut.fit_params['bkg_slope'].value < 0.6))
        self.assertTrue(bool(self.my_linecut.fit_params['bkg_intercept'].value > 0.9) and bool(self.my_linecut.fit_params['bkg_intercept'].value < 1.1))

        self.assertTrue(type(self.my_linecut.plot_fitted(engine='hv') == hv.Curve))
        self.assertTrue(type(self.my_linecut.plot_fitted() == go.Figure))

    def test_other_methods(self):
        ''' Test the other methods of the Linecut class '''
        x = np.linspace(0, 2.5, 207)
        y = np.random.normal(0, 0.1, x.shape) +1
        bg_df = pd.DataFrame({'q':x,'intensity': y})
        my_linecut_bg_subtract = self.my_linecut.subtract_background(bg_df)
        self.assertTrue(len(my_linecut_bg_subtract.data) == 103)

class TestPolar_linecut(unittest.TestCase):
    ''' Test the Polar_linecut class '''
    def setUp(self):
        chi = np.linspace(0, 90, 103)
        amplitude = 10
        center = 0
        sigma = 10
        y = (amplitude/(sigma*np.sqrt(2*np.pi))) * np.exp(-0.5 * (chi - center) ** 2 / sigma ** 2)+1+0.5*chi
        # add noise to the y vector
        y += np.random.normal(0, 0.1, y.shape)                        
        self.my_polar_linecut = Polar_linecut(data = pd.DataFrame({'chi':chi,'intensity': y}))

    def test_attributes(self):
        ''' Test the attributes of the Polar_linecut class '''
        self.assertTrue(len(self.my_polar_linecut.data) == 103)
        self.assertTrue('chi' in self.my_polar_linecut.data.columns)
        self.assertTrue('intensity' in self.my_polar_linecut.data.columns)
        self.assertTrue(len(self.my_polar_linecut.chi) == 103)
        self.assertTrue(len(self.my_polar_linecut.intensity) == 103)
    
    def test_plotting(self):
        ''' Test the plotting methods of the Polar_linecut class '''
        figure = self.my_polar_linecut.plot()
        self.assertTrue(type(figure) == go.Figure)

        figure = self.my_polar_linecut.plot(engine='hv')
        self.assertTrue(type(figure) == hv.Curve)
    