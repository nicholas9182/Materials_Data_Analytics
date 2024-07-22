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
        generic_cation = Cation(name='NA+')
        self.assertTrue(generic_cation.name == 'Sodium')
        self.assertTrue(generic_cation.formula == 'NA+')
        self.assertTrue(generic_cation.charge == 1)

    def test_cation_name_2(self):
        generic_cation = Cation(name='Sodium')
        self.assertTrue(generic_cation.name == 'Sodium')
        self.assertTrue(generic_cation.formula == 'NA+')
        self.assertTrue(generic_cation.charge == 1)
        
    def test_cation_charge_name_3(self):
        generic_cation = Cation(name='na+')
        self.assertTrue(generic_cation.name == 'Sodium')
        self.assertTrue(generic_cation.formula == 'NA+')
        self.assertTrue(generic_cation.charge == 1)

    def test_default_charge(self):
        generic_cation = Cation(name='Sodium')
        self.assertTrue(generic_cation.charge == 1)

    def test_default_charge_2(self):
        generic_cation = Cation(name='NA+')
        self.assertTrue(generic_cation.charge == 1)


class TestAnion(unittest.TestCase):

    def test_anion_name(self):
        generic_anion = Anion(name='Cl-')
        self.assertTrue(generic_anion.name == 'Chloride')
        self.assertTrue(generic_anion.formula == 'CL-')
        self.assertTrue(generic_anion.charge == -1)

    def test_anion_name_2(self):
        generic_anion = Anion(name='Chloride')
        self.assertTrue(generic_anion.name == 'Chloride')
        self.assertTrue(generic_anion.formula == 'CL-')
        self.assertTrue(generic_anion.charge == -1)
        
    def test_anion_charge_name_3(self):
        generic_anion = Anion(name='cl-')
        self.assertTrue(generic_anion.name == 'Chloride')
        self.assertTrue(generic_anion.formula == 'CL-')
        self.assertTrue(generic_anion.charge == -1)

    def test_default_charge(self):
        generic_anion = Anion(name='Chloride')
        self.assertTrue(generic_anion.charge == -1)

    def test_default_charge_2(self):
        generic_anion = Anion(name='CL-')
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


class TestElectrolyte(unittest.TestCase):

    def test_electrolyte_pH(self):
        cation = Cation(name='Sodium')
        anion = Anion(name='Chloride')
        solute = Solute(name='Oxygen')
        solvent = Solvent(name='Water')
        electrolyte = Electrolyte(solvent=solvent, cation=cation, anion=anion, pH=7, solute=solute, 
                                  temperature=298, concentrations={cation: 0.001, anion: 0.001, solute: 0.001})
        self.assertTrue(electrolyte.pH == 7)
        self.assertTrue(electrolyte.temperature == 298)
        self.assertTrue(electrolyte.solvent.name == 'Water')
        self.assertTrue(electrolyte.cation.name == 'Sodium')
        self.assertTrue(electrolyte.anion.name == 'Chloride')
        self.assertTrue(electrolyte.solute.name == 'Oxygen')
        self.assertTrue(electrolyte._concentrations == {cation: 0.001, anion: 0.001, solute: 0.001})
        self.assertTrue(type(electrolyte._cation == list[Cation]))
        self.assertTrue(type(electrolyte._anion == list[Anion]))

    def test_electrolyte_multiple_ions(self):
        cation1 = Cation(name='Sodium')
        cation2 = Cation(name='Potassium')
        anion1 = Anion(name='Chloride')
        anion2 = Anion(name='Bromide')
        solute = Solute(name='Oxygen')
        solvent = Solvent(name='Water')
        electrolyte = Electrolyte(solvent=solvent, cation=[cation1, cation2], anion=[anion1, anion2], pH=7, solute=solute,
                                    temperature=298, concentrations={cation1: 0.001, cation2: 0.001, anion1: 0.001, anion2: 0.001, solute: 0.001})
        self.assertTrue(electrolyte.pH == 7)
        self.assertTrue(electrolyte.temperature == 298)
        self.assertTrue(electrolyte.solvent.name == 'Water')
        self.assertTrue(electrolyte.cation[0].name == 'Sodium')
        self.assertTrue(electrolyte.cation[1].name == 'Potassium')
        self.assertTrue(electrolyte.anion[0].name == 'Chloride')
        self.assertTrue(electrolyte.anion[1].name == 'Bromide')
        self.assertTrue(electrolyte.solute.name == 'Oxygen')
        self.assertTrue(electrolyte._concentrations == {cation1: 0.001, cation2: 0.001, anion1: 0.001, anion2: 0.001, solute: 0.001})
        self.assertTrue(type(electrolyte._cation == list[Cation]))
        self.assertTrue(type(electrolyte._anion == list[Anion]))


class TestPolymer(unittest.TestCase):

    def test_ntype_name(self):
        ntype = NType(name='BBL')
        self.assertTrue(ntype.name == 'BBL')
        self.assertTrue(ntype._name == 'bbl')

    def test_ntype_name_2(self):
        ntype = NType(name='bbl')
        self.assertTrue(ntype.name == 'BBL')
        self.assertTrue(ntype._name == 'bbl')

    def test_ptype_name(self):
        ptype = PType(name='PEDOT')
        self.assertTrue(ptype.name == 'PEDOT')
        self.assertTrue(ptype._name == 'pedot')

    def test_ptype_name_2(self):
        ptype = PType(name='pedot')
        self.assertTrue(ptype.name == 'PEDOT')
        self.assertTrue(ptype._name == 'pedot')

    def test_generic_polymer(self):
        polymer = Polymer(name='P3HT')
        self.assertTrue(polymer.name == 'P3HT')
        self.assertTrue(polymer._name == 'p3ht')

    def test_generic_polymer_2(self):
        polymer = Polymer(name='p3ht')
        self.assertTrue(polymer.name == 'P3HT')
        self.assertTrue(polymer._name == 'p3ht')
        