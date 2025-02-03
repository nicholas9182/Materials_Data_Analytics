import unittest
import tracemalloc
import pandas as pd
import numpy as np
import plotly.express as px
from Materials_Data_Analytics.continuum_modelling.microkinetic_modelling import OxygenReductionModel, ECpD
from Materials_Data_Analytics.materials.electrolytes import Electrolyte
from Materials_Data_Analytics.materials.solutes import Solute, MolecularOxygen
from Materials_Data_Analytics.materials.ions import Cation, Anion
from Materials_Data_Analytics.materials.solvents import Solvent
from Materials_Data_Analytics.materials.polymers import NType


class TestOrrModelHighPh(unittest.TestCase):

    def setUp(self):
        self.my_polymer = NType('BBL')
        self.na_cation = Cation(name='Na+')
        self.cl_anion = Anion(name='Cl-')
        self.water_solvent = Solvent('water')
        self.oxygen_solute = MolecularOxygen()

        self.my_electrolye = Electrolyte(solvent=self.water_solvent, 
                                        cation=self.na_cation, 
                                        anion=self.cl_anion, 
                                        concentrations={self.na_cation: 0.1, self.cl_anion: 0.1, self.oxygen_solute: 0.00008}, 
                                        solute=self.oxygen_solute, 
                                        pH=14.2, 
                                        temperature=298,
                                        diffusivities={self.oxygen_solute: 0.000019},
                                        viscosity=0.01)
        
        self.my_ORR_model = OxygenReductionModel(electrolyte=self.my_electrolye, polymer=self.my_polymer, rotation_rate=1600)
        self.my_polymer = NType('BBL', formal_reduction_potential=-0.3159)
        self.my_ECpD_model = ECpD(electrolyte=self.my_electrolye, polymer=self.my_polymer, rotation_rate=1600)

    def test_orr_model(self):
        """ Test the ORR model """
        self.assertTrue(self.my_ORR_model.pH == 14.2)
        self.assertTrue(self.my_ORR_model.rotation_rate == 1600)
        self.assertTrue(self.my_ORR_model.temperature == 298)
        self.assertTrue(self.my_ORR_model.cation == self.na_cation)
        self.assertTrue(self.my_ORR_model.anion == self.cl_anion)
        self.assertTrue(self.my_ORR_model.solvent == self.water_solvent)
        self.assertTrue(self.my_ORR_model.solute == self.oxygen_solute)
        self.assertTrue(self.my_ORR_model.x == 1)
        self.assertTrue(self.my_ORR_model.diffusivities[self.oxygen_solute] == 0.000019)
        self.assertTrue(self.my_ORR_model.viscosity == 0.01)
        self.assertTrue(round(self.my_ORR_model.diffusion_layer_thickness, 4) == 0.0015)
        self.assertTrue(round(self.my_ORR_model._o2._diffusion_layer_thickness, 4) == 0.0015)
        self.assertTrue(round(self.my_ORR_model.mass_transfer_coefficient, 4) == 0.0123)

    def test_ECpD_model(self):
        """ Test the ECpD model """
        self.assertTrue(self.my_ECpD_model.pH == 14.2)
        self.assertTrue(self.my_ECpD_model.rotation_rate == 1600)
        self.assertTrue(self.my_ECpD_model.temperature == 298)
        self.assertTrue(self.my_ECpD_model.cation == self.na_cation)
        self.assertTrue(self.my_ECpD_model.anion == self.cl_anion)
        self.assertTrue(self.my_ECpD_model.solvent == self.water_solvent)
        self.assertTrue(self.my_ECpD_model.solute == self.oxygen_solute)
        self.assertTrue(self.my_ECpD_model.x == 1)

    def test_k3_calculation(self):
        """ Test the k3 calculation """
        k3 = self.my_ECpD_model.calculate_k3()
        k3 = round(k3, 0)
        self.assertTrue(k3 == 16380)

    def test_solve_parameters(self):
        """ Test the solve parameters method """
        my_ECpD_model = ECpD(electrolyte=self.my_electrolye, polymer=self.my_polymer, rotation_rate=1600)
        parameters = my_ECpD_model.solve_parameters(E=0.7, k01=np.exp(-5.906), beta=0.4999, kf2=np.exp(1.0807), kf3=np.exp(4.688))
        self.assertTrue(type(parameters) == dict)

    def test_get_e_sweep(self):
        """ Test the get e sweep method """
        my_ECpD_model = ECpD(electrolyte=self.my_electrolye, polymer=self.my_polymer, rotation_rate=1600)
        e_sweep = my_ECpD_model.get_e_sweep(E_max=1, E_min=-1, E_n=200, k01=(-5.908), beta=0.5, kf2=(1.082), kf3=(4.70), guess=[1,0,0,0])
        e_sweep = e_sweep.assign(potential = lambda x: x['potential'] + 0.059 * 14.2) # convert to RHE
        # px.line(e_sweep, x='potential', y=['disk_current_density', 'ring_current_density']).show()
        self.assertTrue(type(e_sweep) == pd.DataFrame)


class TestOrrModelLowPh(unittest.TestCase):

    def setUp(self):
        self.my_polymer = NType('BBL')
        self.na_cation = Cation(name='Na+')
        self.cl_anion = Anion(name='Cl-')
        self.water_solvent = Solvent('water')
        self.oxygen_solute = MolecularOxygen()

        self.my_electrolye = Electrolyte(solvent=self.water_solvent, 
                                        cation=self.na_cation, 
                                        anion=self.cl_anion, 
                                        concentrations={self.na_cation: 0.1, self.cl_anion: 0.1, self.oxygen_solute: 0.00025}, 
                                        solute=self.oxygen_solute, 
                                        pH=0, 
                                        temperature=298,
                                        diffusivities={self.oxygen_solute: 0.000019},
                                        viscosity=0.01)
        
        self.my_ORR_model = OxygenReductionModel(electrolyte=self.my_electrolye, polymer=self.my_polymer, rotation_rate=1600)

    def test_orr_model(self):
        """ Test the ORR model """
        self.assertTrue(self.my_ORR_model.pH == 0)
        self.assertTrue(self.my_ORR_model.rotation_rate == 1600)
        self.assertTrue(self.my_ORR_model.temperature == 298)
        self.assertTrue(self.my_ORR_model.cation == self.na_cation)
        self.assertTrue(self.my_ORR_model.anion == self.cl_anion)
        self.assertTrue(self.my_ORR_model.solvent == self.water_solvent)
        self.assertTrue(self.my_ORR_model.solute == self.oxygen_solute)
        self.assertTrue(self.my_ORR_model.x == 0)
        self.assertTrue(self.my_ORR_model.diffusivities[self.oxygen_solute] == 0.000019)
        self.assertTrue(self.my_ORR_model.viscosity == 0.01)
        self.assertTrue(round(self.my_ORR_model.diffusion_layer_thickness, 4) == 0.0015)
        self.assertTrue(round(self.my_ORR_model._o2._diffusion_layer_thickness, 4) == 0.0015)
        self.assertTrue(round(self.my_ORR_model.mass_transfer_coefficient, 4) == 0.0123)


class TestOrrModelNeutralPh(unittest.TestCase):

    def setUp(self):
        self.my_polymer = NType('BBL')
        self.na_cation = Cation(name='Na+')
        self.cl_anion = Anion(name='Cl-')
        self.water_solvent = Solvent('water')
        self.oxygen_solute = MolecularOxygen()

        self.my_electrolye = Electrolyte(solvent=self.water_solvent, 
                                        cation=self.na_cation, 
                                        anion=self.cl_anion, 
                                        concentrations={self.na_cation: 0.1, self.cl_anion: 0.1, self.oxygen_solute: 0.00025}, 
                                        solute=self.oxygen_solute, 
                                        pH=7, 
                                        temperature=298,
                                        diffusivities={self.oxygen_solute: 0.000019},
                                        viscosity=0.01)
        
        self.my_ORR_model = OxygenReductionModel(electrolyte=self.my_electrolye, polymer=self.my_polymer, rotation_rate=1600)

    def test_orr_model(self):
        """ Test the ORR model """
        self.assertTrue(self.my_ORR_model.pH == 7)
        self.assertTrue(self.my_ORR_model.rotation_rate == 1600)
        self.assertTrue(self.my_ORR_model.temperature == 298)
        self.assertTrue(self.my_ORR_model.cation == self.na_cation)
        self.assertTrue(self.my_ORR_model.anion == self.cl_anion)
        self.assertTrue(self.my_ORR_model.solvent == self.water_solvent)
        self.assertTrue(self.my_ORR_model.solute == self.oxygen_solute)
        self.assertTrue(self.my_ORR_model.x == 2)
        self.assertTrue(self.my_ORR_model.diffusivities[self.oxygen_solute] == 0.000019)
        self.assertTrue(self.my_ORR_model.viscosity == 0.01)
        self.assertTrue(round(self.my_ORR_model.diffusion_layer_thickness, 4) == 0.0015)
        self.assertTrue(round(self.my_ORR_model._o2._diffusion_layer_thickness, 4) == 0.0015)
        self.assertTrue(round(self.my_ORR_model.mass_transfer_coefficient, 4) == 0.0123)
