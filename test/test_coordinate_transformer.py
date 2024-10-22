import unittest
import pandas as pd
from Materials_Data_Analytics.core.coordinate_transformer import CoordinateTransformer


class TestCoordinateTransformer(unittest.TestCase):

    data = pd.DataFrame({
        'x': [-0.931, -0.468, -2.329, -2.972],
        'y': [2.537, 3.513, 2.402, 3.270],
        'z': [0, 0, 0, 0]
    })

    def test_rotate_x(self):
        """ Test the rotation of the x axis """
        rot_data = CoordinateTransformer(self.data).rotate(theta_x=82).data

        correct_data = pd.DataFrame({
            'x': [-0.931, -0.468, -2.329, -2.972],
            'y': [0.353, 0.489, 0.334, 0.455],
            'z': [2.512, 3.479, 2.379, 3.238]
        })

        self.assertTrue(type(rot_data) == pd.DataFrame)
        pd.testing.assert_frame_equal(rot_data, correct_data, atol=0.001)
