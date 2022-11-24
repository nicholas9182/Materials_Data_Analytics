import unittest
import plumed as pl
from analytics.metadynamics.free_energy import FreeEnergyLandscape


class TestFreeEnergyLandscape(unittest.TestCase):

    def test_make_landscape(self):

        data = pl.read_as_pandas("../test_trajectories/ndi_na_binding/HILLS")
        landscape = FreeEnergyLandscape(data)
        self.assertTrue('height' in landscape.hills.columns.to_list())
        self.assertTrue('time' in landscape.hills.columns.to_list())
        self.assertEqual(type(landscape), FreeEnergyLandscape)
