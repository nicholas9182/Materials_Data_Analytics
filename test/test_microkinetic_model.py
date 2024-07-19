import unittest
import tracemalloc
from analytics.continuum_modelling.microkinetic_modelling import MicroKineticModel, OxygenReductionModel

class TestMicrokinetic(unittest.TestCase):

    top_model = MicroKineticModel(pH=7, rotation_rate=1600, bulk_concentration=0.008, temperature=298)
    orr_model = OxygenReductionModel(pH=14.2, rotation_rate=1600, bulk_concentration=0.008, temperature=298)

    def test_top_model_attributes(self):
        self.assertTrue(self.top_model.pH == 7)
        self.assertTrue(self.top_model.rotation_rate == 1600)
        self.assertTrue(self.top_model.bulk_concentration == 0.008)
        self.assertTrue(self.top_model.temperature == 298)
    
    def test_ORR_model_attributes(self):
        self.assertTrue(self.orr_model.pH == 14.2)
        self.assertTrue(self.orr_model.rotation_rate == 1600)
        self.assertTrue(self.orr_model.bulk_concentration == 0.008)
        self.assertTrue(self.orr_model.temperature == 298)
        self.assertTrue(self.orr_model.diffusion_layer_thickness == 0.00154)
        self.assertTrue(self.orr_model.mass_transfer_coefficient == 0.01233)
        self.assertTrue(self.orr_model.potential_formal_O2_HXO2 == -0.071)

        