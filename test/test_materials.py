import unittest
import tracemalloc
from analytics.materials.electrolytes import Electrolyte
from analytics.materials.ions import Ion, Cation, Anion
from analytics.materials.polymers import Polymer, NType, PType
from analytics.materials.solvents import Solvent
from analytics.materials.solutes import Solute


class TestSolvent(unittest.TestCase):
    """
    Test the solvent class
    """
    def test_solvent_name(self):
        generic_solvent = Solvent(name='Water')
        self.assertTrue(generic_solvent.name == 'Water')
        self.assertTrue(generic_solvent.formula == 'H2O')

    def test_solvent_formula(self):
        generic_solvent = Solvent(name='Water')
        self.assertTrue(generic_solvent.formula == 'H2O')
        self.assertTrue(generic_solvent.name == 'Water')

    def test_solvent_name_2(self):
        generic_solvent = Solvent(name='H2O')
        self.assertTrue(generic_solvent.name == 'Water')
        self.assertTrue(generic_solvent.formula == 'H2O')


class TestCation(unittest.TestCase):

    def test_cation_name(self):
        generic_cation = Cation(charge=1, name='NA+')
        self.assertTrue(generic_cation.name == 'Sodium')
        self.assertTrue(generic_cation.formula == 'NA+')
        self.assertTrue(generic_cation.charge == 1)

    def test_cation_name_2(self):
        generic_cation = Cation(charge=1, name='Sodium')
        self.assertTrue(generic_cation.name == 'Sodium')
        self.assertTrue(generic_cation.formula == 'NA+')
        self.assertTrue(generic_cation.charge == 1)
    
    def test_cation_charge_positive(self):
        with self.assertRaises(ValueError):
            Cation(charge=-1, name='Sodium')
        
    def test_cation_charge_name_3(self):
        generic_cation = Cation(charge=1, name='na+')
        self.assertTrue(generic_cation.name == 'Sodium')
        self.assertTrue(generic_cation.formula == 'NA+')
        self.assertTrue(generic_cation.charge == 1)


class TestAnion(unittest.TestCase):

    def test_anion_name(self):
        generic_anion = Anion(charge=-1, name='Cl-')
        self.assertTrue(generic_anion.name == 'Chloride')
        self.assertTrue(generic_anion.formula == 'CL-')
        self.assertTrue(generic_anion.charge == -1)

    def test_anion_name_2(self):
        generic_anion = Anion(charge=-1, name='Chloride')
        self.assertTrue(generic_anion.name == 'Chloride')
        self.assertTrue(generic_anion.formula == 'CL-')
        self.assertTrue(generic_anion.charge == -1)
    
    def test_anion_charge_positive(self):
        with self.assertRaises(ValueError):
            Anion(charge=1, name='Chloride')
        
    def test_anion_charge_name_3(self):
        generic_anion = Anion(charge=-1, name='cl-')
        self.assertTrue(generic_anion.name == 'Chloride')
        self.assertTrue(generic_anion.formula == 'CL-')
        self.assertTrue(generic_anion.charge == -1)


class TestSolute(unittest.TestCase):
    
    def test_solute_name(self):
        generic_solute = Solute(name='O2')
        self.assertTrue(generic_solute.name == 'Oxygen')
        self.assertTrue(generic_solute.formula == 'O2')

    def test_solute_name_2(self):
        generic_solute = Solute(name='Oxygen')
        self.assertTrue(generic_solute.name == 'Oxygen')
        self.assertTrue(generic_solute.formula == 'O2')

    def test_solute_name_3(self):
        generic_solute = Solute(name='o2')
        self.assertTrue(generic_solute.name == 'Oxygen')
        self.assertTrue(generic_solute.formula == 'O2')

