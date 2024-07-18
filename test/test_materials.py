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
    generic_solvent = Solvent(name='generic_solvent')

    def test_solvent_name(self):
        self.assertTrue(self.generic_solvent.name == 'generic_solvent')


class TestCation(unittest.TestCase):

    generic_cation = Cation(charge=1, name = 'Na+')

    def test_cation_charge(self):
        self.assertTrue(self.generic_cation.charge == 1)

    def test_cation_name(self):
        self.assertTrue(self.generic_cation.name == 'Na+')

    def test_cation_charge_negative(self):
        with self.assertRaises(ValueError):
            Cation(charge=-1)


class TestAnion(unittest.TestCase):
    """
    Test the anion class
    """
    generic_anion = Anion(charge=-1, name = 'Cl-')

    def test_anion_name(self):
        self.assertTrue(self.generic_anion.name == 'Cl-')

    def test_anion_charge(self):
        self.assertTrue(self.generic_anion.charge == -1)

    def test_anion_charge_positive(self):
        with self.assertRaises(ValueError):
            Anion(charge=1)


class TestSolute(unittest.TestCase):
    """
    Test the solute class
    """
    generic_solute = Solute(name='O2')

    def test_solute_name(self):
        self.assertTrue(self.generic_solute.name == 'O2')

