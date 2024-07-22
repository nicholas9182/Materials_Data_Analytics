import unittest
import tracemalloc
from analytics.continuum_modelling.microkinetic_modelling import MicroKineticModel
from analytics.materials.electrolytes import Electrolyte
from analytics.materials.solutes import Solute, MolecularOxygen
from analytics.materials.ions import Cation, Anion
from analytics.materials.solvents import Solvent
from analytics.materials.polymers import NType
from analytics.continuum_modelling.microkinetic_modelling import OxygenReductionModel


class TestMicrokinetic(unittest.TestCase):

    na_cation = Cation(name='Na+')
    cl_anion = Anion(name='Cl-')
    water_solvent = Solvent('water')
    oxygen_solute = MolecularOxygen()

    my_electrolye = Electrolyte(solvent=water_solvent, 
                                cation=na_cation, 
                                anion=cl_anion, 
                                concentrations={na_cation: 0.1, cl_anion: 0.1, oxygen_solute: 0.00025}, 
                                solute=oxygen_solute, 
                                pH=7, 
                                temperature=298)
    
    my_polymer = NType('BBL')

    def test_electrolyte(self):
        self.assertTrue(self.my_electrolye.pH == 7)
        self.assertTrue(self.my_electrolye.temperature == 298)
        self.assertTrue(self.my_electrolye.solvent.name == 'Water')
        self.assertTrue(self.my_electrolye.cation.name == 'Sodium')
        self.assertTrue(self.my_electrolye.anion.name == 'Chloride')
        self.assertTrue(self.my_electrolye.solute.name == 'Oxygen')
        self.assertTrue(self.my_electrolye.cation.formula == 'NA+')
        self.assertTrue(self.my_electrolye.anion.formula == 'CL-')
        self.assertTrue(self.my_electrolye.solvent.formula == 'H2O')
        self.assertTrue(self.my_electrolye.concentrations[self.na_cation] == 0.1)
        self.assertTrue(self.my_electrolye.concentrations[self.cl_anion] == 0.1)
        self.assertTrue(self.my_electrolye.concentrations[self.oxygen_solute] == 0.00025)

    def test_polymer(self):
        self.assertTrue(self.my_polymer.name == 'BBL')

    my_microkinetic_model = MicroKineticModel(electrolyte=my_electrolye, polymer=my_polymer, rotation_rate=1600, temperature=298)

    def test_microkinetic_model_attributes(self):
        self.assertTrue(self.my_microkinetic_model.pH == 7)
        self.assertTrue(self.my_microkinetic_model.rotation_rate == 1600)
        self.assertTrue(self.my_microkinetic_model.temperature == 298)
        self.assertTrue(self.my_microkinetic_model.cation == self.na_cation)
        self.assertTrue(self.my_microkinetic_model.anion == self.cl_anion)
        self.assertTrue(self.my_microkinetic_model.solvent == self.water_solvent)


class TestOrrModel(unittest.TestCase):

    na_cation = Cation(name='Na+')
    cl_anion = Anion(name='Cl-')
    water_solvent = Solvent('water')
    oxygen_solute = MolecularOxygen()

    my_electrolye = Electrolyte(solvent=water_solvent, 
                                cation=na_cation, 
                                anion=cl_anion, 
                                concentrations={na_cation: 0.1, cl_anion: 0.1, oxygen_solute: 0.00025}, 
                                solute=oxygen_solute, 
                                pH=7, 
                                temperature=298)
    
    my_polymer = NType('BBL')

    def test_orr_model(self):
        my_ORR_model = OxygenReductionModel(electrolyte=self.my_electrolye, polymer=self.my_polymer, rotation_rate=1600, temperature=298, pH=7)
        self.assertTrue(my_ORR_model.pH == 7)
        self.assertTrue(my_ORR_model.rotation_rate == 1600)
        self.assertTrue(my_ORR_model.temperature == 298)
        self.assertTrue(my_ORR_model.cation == self.na_cation)
        self.assertTrue(my_ORR_model.anion == self.cl_anion)
        self.assertTrue(my_ORR_model.solvent == self.water_solvent)
        self.assertTrue(my_ORR_model.solute == self.oxygen_solute)
        self.assertTrue(my_ORR_model.x == 2)

    def test_orr_low_pH(self):
        my_ORR_model = OxygenReductionModel(electrolyte=self.my_electrolye, polymer=self.my_polymer, rotation_rate=1600, temperature=298, pH=0)
        self.assertTrue(my_ORR_model.x == 0)

    def test_orr_high_pH(self):
        my_ORR_model = OxygenReductionModel(electrolyte=self.my_electrolye, polymer=self.my_polymer, rotation_rate=1600, temperature=298, pH=13)
        self.assertTrue(my_ORR_model.x == 1)

