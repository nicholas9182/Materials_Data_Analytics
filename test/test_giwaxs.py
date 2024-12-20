import unittest
import tracemalloc
import pandas as pd
import numpy as np
import plotly.graph_objects as go
from Materials_Data_Analytics.experiment_modelling.giwaxs import Calibrator
from Materials_Data_Analytics.experiment_modelling.giwaxs import GIWAXSPixelImage, GIWAXSPattern
import plotly.express as px
import plotly as pl


class TestCalibration(unittest.TestCase):
    ''' Test the Calibrator class '''
    my_calibrator = Calibrator.from_poni_file('./test_trajectories/giwaxs/calibration.poni')

    def test_attributes(self):
        ''' Test the attributes of the Calibration class '''
        self.assertEqual(self.my_calibrator.wavelength, 0.09919)
        self.assertEqual(self.my_calibrator.distance, 0.28556)
        self.assertEqual(self.my_calibrator.poni1, 0.21226)
        self.assertEqual(self.my_calibrator.poni2, 0.11560)
        self.assertEqual(self.my_calibrator.rot1, 0.0025858)
        self.assertEqual(self.my_calibrator.rot2, 0.0093694)
        self.assertEqual(self.my_calibrator.rot3, 0.0000000)


# class TestGIWAXSPixelImage(unittest.TestCase):
#     ''' Test the GIWAXS class '''
#     data_SLAC_BL113 = GIWAXSPixelImage.from_SLAC_BL11_3(tif_filepaths= ['./test_trajectories/giwaxs/GIWAXS_image_SLAC_1.tif'])
#     ai = Calibrator.from_poni_file('./test_trajectories/giwaxs/calibration.poni')

#     def test_from_SLAC_BL113(self):
#         ''' Test the attributes of the GIWAXS class '''
#         self.assertTrue(self.data_SLAC_BL113.image.shape == (3072, 3072))
#         self.assertTrue(np.round(self.data_SLAC_BL113.image[5][2], 3) == 51.0)
#         self.assertTrue(np.round(self.data_SLAC_BL113.image[15][27], 3) == 49.0)
#         self.assertTrue(np.round(self.data_SLAC_BL113.image[257][43], 3) == 60.0)
#         self.assertTrue(self.data_SLAC_BL113.incidence_angle == 0.12)
#         self.assertTrue(self.data_SLAC_BL113.exposure_time == 120.0)

#     def test_mask(self):
#         ''' Test the mask method of the GIWAXS class '''
#         self.data_SLAC_BL113.apply_mask(mask_path='./test_trajectories/giwaxs/mask.tif')
#         self.assertTrue(np.round(self.data_SLAC_BL113.image[5][2], 3) == 51.0)
#         self.assertTrue(np.round(self.data_SLAC_BL113.image[15][27], 3) == 49.0)
#         self.assertTrue(np.round(self.data_SLAC_BL113.image[257][43], 3) == 60.0)
#         self.assertTrue(np.isnan(self.data_SLAC_BL113.image[3000][43]))
        
#     def test_transform(self):
#         ''' Test the transform method of the GIWAXS class '''
#         my_giwaxs_pattern = self.data_SLAC_BL113.get_giwaxs_pattern(calibrator = self.ai,
#                                                                     qxy_range = (-3, 3),
#                                                                     qz_range = (0, 3),
#                                                                     q_range = (0, 3),
#                                                                     chi_range = (-95, 95),
#                                                                     pixel_q = 100,
#                                                                     pixel_chi = 80,
#                                                                     correct_solid_angle = True,
#                                                                     polarization_factor = None,
#                                                                     unit = 'A')
        
#         self.assertTrue(type(my_giwaxs_pattern == GIWAXSPattern))
#         self.assertTrue('qxy' in my_giwaxs_pattern.data_reciprocal.columns)
#         self.assertTrue('qz' in my_giwaxs_pattern.data_reciprocal.columns)
#         self.assertTrue('intensity' in my_giwaxs_pattern.data_reciprocal.columns)
#         self.assertTrue('chi' in my_giwaxs_pattern.data_polar.columns)
#         self.assertTrue('q' in my_giwaxs_pattern.data_polar.columns)
#         self.assertTrue('intensity' in my_giwaxs_pattern.data_polar.columns)
    

# class TestGIWAXSPattern(unittest.TestCase):
#     ''' Test the GIWAXSPattern class '''                             
#     ai = Calibrator.from_poni_file('./test_trajectories/giwaxs/calibration.poni')

#     data_SLAC_BL113 = (GIWAXSPixelImage
#                         .from_SLAC_BL11_3(tif_filepaths= ['./test_trajectories/giwaxs/GIWAXS_image_SLAC_1.tif'])
#                         .get_giwaxs_pattern(calibrator = ai,
#                                             qxy_range = (-3, 3),
#                                             qz_range = (0, 3),
#                                             q_range = (0, 3),
#                                             chi_range = (-95, 95),
#                                             pixel_q = 100,
#                                             pixel_chi = 80,
#                                             correct_solid_angle = True,
#                                             polarization_factor = None,
#                                             unit = 'A'
#                                             )
#                         )
    
#     def test_attributes(self):
#         ''' Test the attributes of the GIWAXSPattern class '''
#         self.assertTrue(self.data_SLAC_BL113.qxy.shape == (100,))
#         self.assertTrue(self.data_SLAC_BL113.qz.shape == (100,))
#         self.assertTrue(len(self.data_SLAC_BL113.data_polar) == 15100)
#         self.assertTrue(len(self.data_SLAC_BL113.data_reciprocal) == 10000)

#     def test_plot_map_reciprocal_space_contour(self):
#         ''' Test the plot_map_reciprocal_space method of the GIWAXSPattern class '''
#         figure = self.data_SLAC_BL113.plot_reciprocal_map_contour(intensity_lower_cuttoff = 30)
#         # figure.show()
#         self.assertTrue(type(figure) == go.Figure)   

#     def test_plot_map_reciprocal_space_image(self):
#         ''' Test the plot_map_reciprocal_space method of the GIWAXSPattern class '''
#         figure = self.data_SLAC_BL113.plot_reciprocal_map(width=800, height=500, intensity_lower_cuttoff = 30)
#         # figure.show()
#         self.assertTrue(type(figure) == go.Figure)   

#     def test_plot_map_polar_space_contour(self):
#         ''' Test the plot_map_reciprocal_space method of the GIWAXSPattern class '''
#         figure = self.data_SLAC_BL113.plot_polar_map_contour(width=800, height=500, intensity_lower_cuttoff = 30)
#         # figure.show()
#         self.assertTrue(type(figure) == go.Figure)     

#     def test_plot_map_polar_space_image(self):
#         ''' Test the plot_map_reciprocal_space method of the GIWAXSPattern class '''
#         figure = self.data_SLAC_BL113.plot_polar_map(intensity_lower_cuttoff = 30)
#         # figure.show()
#         self.assertTrue(type(figure) == go.Figure)                              

#     def test_extract_profile_with_range(self):
#         ''' Test the extract_profile method of the GIWAXSPattern class using a range '''
#         profile_df = self.data_SLAC_BL113.get_linecut(chi = (20, 40), q_range = (0.2, 2))
#         # px.line(profile_df, x = 'q', y = 'intensity').show()
#         self.assertTrue(type(profile_df) == pd.DataFrame)
#         self.assertTrue('q' in profile_df.columns)
#         self.assertTrue('intensity' in profile_df.columns)
#         self.assertTrue(profile_df.round(3)['q'].iloc[10] == 0.525)
#         self.assertTrue(profile_df.round(3)['q'].iloc[20] == 0.825)
#         self.assertTrue(profile_df.round(0)['intensity'].iloc[10] == 91.0)
#         self.assertTrue(profile_df.round(0)['intensity'].iloc[20] == 88.0)
#         self.assertTrue(profile_df.round(0)['chi'].iloc[30] == 30.0)

#     def test_extract_profile_with_qnone(self):
#         ''' Test the extract_profile method of the GIWAXSPattern class using None as q range '''
#         profile_df = self.data_SLAC_BL113.get_linecut(chi = (20, 40))
#         # px.line(profile_df, x = 'q', y = 'intensity').show()
#         self.assertTrue(type(profile_df) == pd.DataFrame)
#         self.assertTrue('q' in profile_df.columns)
#         self.assertTrue('intensity' in profile_df.columns)
#         self.assertTrue(profile_df.round(3)['q'].iloc[10] == 0.315)
#         self.assertTrue(profile_df.round(3)['q'].iloc[20] == 0.615)
#         self.assertTrue(profile_df.round(0)['intensity'].iloc[10] == 257.0)
#         self.assertTrue(profile_df.round(0)['intensity'].iloc[20] == 93.0)
#         self.assertTrue(profile_df.round(0)['chi'].iloc[30] == 30.0)

#     def test_extract_profile_with_chi_float(self):
#         profile_df = self.data_SLAC_BL113.get_linecut(chi = 30)
#         # px.line(profile_df, x = 'q', y = 'intensity').show()
#         self.assertTrue(type(profile_df) == pd.DataFrame)
#         self.assertTrue('q' in profile_df.columns)
#         self.assertTrue('intensity' in profile_df.columns)
#         self.assertTrue(profile_df.round(3)['q'].iloc[10] == 0.315)
#         self.assertTrue(profile_df.round(3)['q'].iloc[20] == 0.615)
#         self.assertTrue(profile_df.round(0)['intensity'].iloc[10] == 196.0)
#         self.assertTrue(profile_df.round(0)['intensity'].iloc[20] == 91.0)
#         self.assertTrue(profile_df.round(0)['chi'].iloc[30] == 30.0)

#     def test_plot_linecut(self):
#         ''' Test the plot_linecut method of the GIWAXSPattern class '''
#         output = self.data_SLAC_BL113.plot_linecut(chi = (20, 40), q_range = (0.2, 2), template = 'presentation')
#         # output.show()
#         self.assertTrue(type(output) == pl.graph_objs._figure.Figure)
