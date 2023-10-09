import unittest
import tracemalloc
import pandas as pd
from analytics.quantum_chemistry.gaussian import GaussianParser
tracemalloc.start()


class TestGaussianParser(unittest.TestCase):
    """
    Test class for the GaussianParser class
    """
    pedot_log = GaussianParser("../test_trajectories/pedot_raman/step1.log")

    def test_charge(self):
        self.assertTrue(self.pedot_log.charge == 0)

    def test_multiplicity(self):
        self.assertTrue(self.pedot_log.multiplicity == 1)

    def test_keywords(self):
        self.assertTrue(type(self.pedot_log.keywords) == list)

    def test_raman(self):
        self.assertTrue(self.pedot_log.raman is True)

    def test_opt(self):
        self.assertTrue(self.pedot_log.opt is True)

    def test_complete(self):
        self.assertTrue(self.pedot_log.complete is True)

    def test_atomcount(self):
        self.assertTrue(self.pedot_log.atomcount == 28)

    def test_raman(self):
        raman_frequencies = self.pedot_log.get_raman_frequencies()
        self.assertTrue(type(raman_frequencies) == pd.DataFrame)
        self.assertTrue(len(raman_frequencies) == 78)

    def test_raman_spectra(self):
        raman_spectra = self.pedot_log.get_raman_spectra()
        self.assertTrue(len(raman_spectra) == 3000)
