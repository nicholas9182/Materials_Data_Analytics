import unittest
import tracemalloc
import holoviews as hv
import pandas as pd
hv.extension('bokeh')
import numpy as np
from Materials_Data_Analytics.experiment_modelling.Xray_diffraction import Calibrator
from Materials_Data_Analytics.experiment_modelling.Xray_diffraction import GIWAXSPixelImage, GIWAXSPattern


class TestCalibration(unittest.TestCase):
    ''' Test the Calibrator class '''
    my_calibrator = Calibrator.from_poni_file('./test_trajectories/diffraction/calibration.poni')

    def test_attributes(self):
        ''' Test the attributes of the Calibration class '''
        self.assertEqual(self.my_calibrator._wavelegth, 9.91873587465602e-11)
        self.assertEqual(self.my_calibrator.distance, 0.2855582907344562)
        self.assertEqual(self.my_calibrator.poni1, 0.21226031843207446)
        self.assertEqual(self.my_calibrator.poni2, 0.11560087180581186)
        self.assertEqual(self.my_calibrator.rot1, 0.002585838591972173)
        self.assertEqual(self.my_calibrator.rot2, 0.009369430517925924)
        self.assertEqual(self.my_calibrator.rot3,3.116397003897587e-10)


class TestGIWAXSPixelImage(unittest.TestCase):
    ''' Test the GIWAXS class '''
    data_SLAC_BL113 = GIWAXSPixelImage.from_SLAC_BL11_3(tif_filepaths= ['./test_trajectories/diffraction/GIWAXS_image_SLAC_1.tif',
                                                                  './test_trajectories/diffraction/GIWAXS_image_SLAC_2.tif'])

    def test_from_SLAC_BL113(self):
        ''' Test the attributes of the GIWAXS class '''
        self.assertTrue(self.data_SLAC_BL113.image_row.shape == (3072, 3072))
        self.assertTrue(self.data_SLAC_BL113.image_row[5][2] == 51.004597211078504)
        self.assertTrue(self.data_SLAC_BL113.image_row[15][27] == 45.97566237392315)
        self.assertTrue(self.data_SLAC_BL113.image_row[257][43] == 60.510200909349756)
        self.assertTrue(self.data_SLAC_BL113.incidence_angle == 0.12)
        self.assertTrue(self.data_SLAC_BL113.exposure_time == 120.0)

    def test_mask(self):
        ''' Test the mask method of the GIWAXS class '''
        self.data_SLAC_BL113.apply_mask(maskpath='./test_trajectories/diffraction/mask.tif')
        self.assertTrue(self.data_SLAC_BL113.image_masked[5][2] == 51.004597211078504)
        self.assertTrue(self.data_SLAC_BL113.image_masked[15][27] == 45.97566237392315)
        self.assertTrue(self.data_SLAC_BL113.image_masked[257][43] == 60.510200909349756)
        self.assertTrue(np.isnan(self.data_SLAC_BL113.image_masked[3000][43]))
        
    def test_transform(self):
        ''' Test the transform method of the GIWAXS class '''
        ai = Calibrator.from_poni_file('./test_trajectories/diffraction/calibration.poni')

        self.data_SLAC_BL113.transform(calibrator = ai,
                                    qxy_range = (-3, 3),
                                    qz_range = (0, 3),
                                    q_range = (0, 3),
                                    chi_range = (-95, 95),
                                    pixel_Q = 500,
                                    pixel_chi = 360,
                                    correct_solid_angle = True,
                                    polarization_factor = None,
                                    unit = 'A')
        
        self.assertTrue(type(self.data_SLAC_BL113 == GIWAXSPixelImage))
        self.assertTrue('qxy' in self.data_SLAC_BL113.data_reciprocal.columns)
        self.assertTrue('qz' in self.data_SLAC_BL113.data_reciprocal.columns)
        self.assertTrue('I' in self.data_SLAC_BL113.data_reciprocal.columns)
        self.assertTrue('chi' in self.data_SLAC_BL113.data_polar.columns)
        self.assertTrue('Q' in self.data_SLAC_BL113.data_polar.columns)
        self.assertTrue('I' in self.data_SLAC_BL113.data_polar.columns)
    
    def test_to_GIWAXSPattern(self):
        ''' Test the to_GIWAXSPattern method of the GIWAXS class '''
        ai = Calibrator.from_poni_file('./test_trajectories/diffraction/calibration.poni')
        output = self.data_SLAC_BL113.transform(calibrator = ai, qxy_range = (-3, 3), qz_range = (0, 3),
                                      q_range = (0, 3), chi_range = (-95, 95), pixel_Q = 500, pixel_chi = 360,
                                      correct_solid_angle = True, polarization_factor = None, unit = 'A').to_GIWAXSPattern()
        
        self.assertTrue(type(output) == GIWAXSPattern)

class TestGIWAXSPattern(unittest.TestCase):
    ''' Test the GIWAXSPattern class '''                             
    ai = Calibrator.from_poni_file('./test_trajectories/diffraction/calibration.poni')
    
    data_SLAC_BL113 = (GIWAXSPixelImage
                       .from_SLAC_BL11_3(tif_filepaths= ['./test_trajectories/diffraction/GIWAXS_image_SLAC_1.tif', './test_trajectories/diffraction/GIWAXS_image_SLAC_2.tif'])
                       .transform(calibrator = ai,
                                  qxy_range = (-3, 3),
                                  qz_range = (0, 3),
                                  q_range = (0, 3),
                                  chi_range = (-95, 95),
                                  pixel_Q = 500,
                                  pixel_chi = 360,
                                  correct_solid_angle = True,
                                  polarization_factor = None,
                                  unit = 'A')
                       .to_GIWAXSPattern()
                       )
    
    def test_attributes(self):
        ''' Test the attributes of the GIWAXSPattern class '''
        self.assertTrue(self.data_SLAC_BL113.qxy.shape == (500,))
        self.assertTrue(self.data_SLAC_BL113.qz.shape == (500,))
        self.assertTrue(self.data_SLAC_BL113.Intensity_reciprocal.shape == (500, 500))
        self.assertTrue(self.data_SLAC_BL113.chi.shape == (360,))
        self.assertTrue(self.data_SLAC_BL113.Q.shape == (500,))
        self.assertTrue(self.data_SLAC_BL113.Intensity_polar.shape == (500, 360))

    def test_plot_map_reciprocal_space(self):
        ''' Test the plot_map_reciprocal_space method of the GIWAXSPattern class '''
        output = self.data_SLAC_BL113.plot_map_reciprocal_space()
        self.assertTrue(type(output) == hv.Image)   

    def test_plot_map_polar_space(self):
        ''' Test the plot_map_reciprocal_space method of the GIWAXSPattern class '''
        output = self.data_SLAC_BL113.plot_map_polar_space()
        self.assertTrue(type(output) == hv.Image)                                   

    def test_extract_profile(self):
        ''' Test the extract_profile method of the GIWAXSPattern class '''
        profile_df = self.data_SLAC_BL113.extract_profile(chi_min = 20, chi_max = 40, q_range = (0.2, 2))
        self.assertTrue(type(profile_df) == pd.DataFrame)
        self.assertTrue('Q' in profile_df.columns)
        self.assertTrue('I' in profile_df.columns)
        
    def test_plot_profile(self):
        ''' Test the plot_profile method of the GIWAXSPattern class '''
        output = self.data_SLAC_BL113.plot_profile(chi_min = 20, chi_max = 40, q_range = (0.2, 2))
        self.assertTrue(type(output) == hv.Curve)
