import unittest
import tracemalloc
from analytics.core.universe import Universe
from analytics.metadynamics.free_energy import FreeEnergySpace
tracemalloc.start()


class TestUniverse(unittest.TestCase):

    def test_universe_creation(self):
        """
        testing if the universe creation works as expected
        :return:
        """
        my_universe = Universe()
        self.assertTrue(type(my_universe) == Universe)

    def test_universe_attributes(self):
        """
        testing if the attributes work as expected
        :return:
        """
        landscape = FreeEnergySpace("./test_trajectories/ndi_na_binding/HILLS")
        my_universe = Universe(fes=landscape)
        self.assertTrue(my_universe._fes == [landscape])

