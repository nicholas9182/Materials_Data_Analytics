import unittest
import tracemalloc
from Materials_Data_Analytics.materials.electrolytes import Electrolyte
from Materials_Data_Analytics.materials.ions import Ion, Cation, Anion
from Materials_Data_Analytics.materials.polymers import Polymer, NType, PType
from Materials_Data_Analytics.materials.solvents import Solvent
from Materials_Data_Analytics.materials.solutes import MolecularOxygen, Solute


class TestSolvent(unittest.TestCase):
    """
    Test the solvent class
    """
    def setUp(self):
        self.generic_solvent = Solvent(name='Water')
        self.custom_solvent = Solvent.from_custom_inputs(name='Some_solvent', formula='SS', pH=7.5)

    def test_solvent_name(self):
        """ Test the name attribute of the solvent class """
        self.assertTrue(self.generic_solvent.name == 'Water')
        self.assertTrue(self.generic_solvent.formula == 'H2O')
        self.assertTrue(self.generic_solvent.pH == 7)

    def test_custom_solvent(self):
        """ Test the custom solvent """
        self.assertTrue(self.custom_solvent.name == 'Some_solvent')
        self.assertTrue(self.custom_solvent.formula == 'SS')
        self.assertTrue(self.custom_solvent.pH == 7.5)


class TestCation(unittest.TestCase):

    def setUp(self):
        self.generic_cation_1 = Cation(name='NA+')
        self.generic_cation_2 = Cation(name='Sodium')
        self.generic_cation_3 = Cation(name='na+')
        self.custom_cation_1 = Cation.from_custom_inputs(name='Some_cation', formula='SC+', charge=1.5)

    def test_cation_properties(self):
        """ Test the name attribute of the cation class when constructed with formula """
        self.assertTrue(self.generic_cation_1.name == 'Sodium')
        self.assertTrue(self.generic_cation_1.formula == 'NA+')
        self.assertTrue(self.generic_cation_1.charge == 1)
        self.assertTrue(self.generic_cation_2.name == 'Sodium')
        self.assertTrue(self.generic_cation_2.formula == 'NA+')
        self.assertTrue(self.generic_cation_2.charge == 1)
        self.assertTrue(self.generic_cation_3.name == 'Sodium')
        self.assertTrue(self.generic_cation_3.formula == 'NA+')
        self.assertTrue(self.generic_cation_3.charge == 1)
        self.assertTrue(self.custom_cation_1.name == 'Some_cation')
        self.assertTrue(self.custom_cation_1.formula == 'SC+')
        self.assertTrue(self.custom_cation_1.charge == 1.5)


class TestAnion(unittest.TestCase):

    def setUp(self):
        self.generic_anion_1 = Anion(name='Cl-')
        self.generic_anion_2 = Anion(name='Chloride')
        self.generic_anion_3 = Anion(name='cl-')

    def test_anion_name(self):
        """ Test the name attribute of the anion class when constructed with formula """
        self.assertTrue(self.generic_anion_1.name == 'Chloride')
        self.assertTrue(self.generic_anion_1.formula == 'CL-')
        self.assertTrue(self.generic_anion_1.charge == -1)
        self.assertTrue(self.generic_anion_2.name == 'Chloride')
        self.assertTrue(self.generic_anion_2.formula == 'CL-')
        self.assertTrue(self.generic_anion_2.charge == -1)
        self.assertTrue(self.generic_anion_3.name == 'Chloride')
        self.assertTrue(self.generic_anion_3.formula == 'CL-')
        self.assertTrue(self.generic_anion_3.charge == -1)


class TestSolute(unittest.TestCase):

    def setUp(self):
        self.generic_solute_1 = MolecularOxygen()
        self.generic_solute_2 = Solute(name='H2O2')
        self.generic_solute_3 = Solute(name='HO2')
    
    def test_attributes(self):
        """ Test the name attribute of the solute class """
        self.assertTrue(self.generic_solute_1.name == 'Oxygen')
        self.assertTrue(self.generic_solute_1.formula == 'O2')
        self.assertTrue(self.generic_solute_1.formal_reduction_potentials == {"O2_superoxide": -0.160})
        self.assertTrue(self.generic_solute_1.standard_reduction_potentials == {"h2o2": 0.695})
        self.assertTrue(self.generic_solute_2.pka == 11.7)
        self.assertTrue(self.generic_solute_3.pka == 4.88)


class TestElectrolyte(unittest.TestCase):

    def setUp(self):
        self.cation = Cation(name='Sodium')
        self.anion = Anion(name='Chloride')
        self.solute = MolecularOxygen()
        self.solvent = Solvent(name='Water')

        self.electrolyte = Electrolyte(solvent=self.solvent, 
                                       cation=self.cation, 
                                       anion=self.anion, 
                                       pH=7, 
                                       solute=self.solute, 
                                       temperature=298, 
                                       concentrations={self.cation: 0.001, self.anion: 0.001, self.solute: 0.001})

        self.cation2 = Cation(name='Potassium')
        self.anion2 = Anion(name='Bromide')

        self.electrolyte2 = Electrolyte(solvent=self.solvent, 
                                   cation=[self.cation, self.cation2], 
                                   anion=[self.anion, self.anion2], 
                                   pH=7, 
                                   solute=self.solute, 
                                   temperature=298, 
                                   concentrations={self.cation: 0.001, self.cation2: 0.001, self.anion: 0.001, self.anion2: 0.001, self.solute: 0.001})

    def test_attributes(self):
        """ test attributes of the electrolyte class """
        self.assertTrue(self.electrolyte.pH == 7)
        self.assertTrue(self.electrolyte.temperature == 298)
        self.assertTrue(self.electrolyte.solvent.name == 'Water')
        self.assertTrue(self.electrolyte.cation.name == 'Sodium')
        self.assertTrue(self.electrolyte.anion.name == 'Chloride')
        self.assertTrue(self.electrolyte.solute.name == 'Oxygen')
        self.assertTrue(self.electrolyte._concentrations == {self.cation: 0.001, self.anion: 0.001, self.solute: 0.001})
        self.assertTrue(self.electrolyte.concentrations == {"Sodium": 0.001, "Chloride": 0.001, "Oxygen": 0.001})
        self.assertTrue(type(self.electrolyte._cation == list[Cation]))
        self.assertTrue(type(self.electrolyte._anion == list[Anion]))
        self.assertTrue(self.electrolyte2.pH == 7)
        self.assertTrue(self.electrolyte2.temperature == 298)
        self.assertTrue(self.electrolyte2.solvent.name == 'Water')
        self.assertTrue(self.electrolyte2.cation[0].name == 'Sodium')
        self.assertTrue(self.electrolyte2.cation[1].name == 'Potassium')
        self.assertTrue(self.electrolyte2.anion[0].name == 'Chloride')
        self.assertTrue(self.electrolyte2.anion[1].name == 'Bromide')
        self.assertTrue(self.electrolyte2.solute.name == 'Oxygen')
        self.assertTrue(self.electrolyte2._concentrations == {self.cation: 0.001, self.cation2: 0.001, self.anion: 0.001, self.anion2: 0.001, self.solute: 0.001})
        self.assertTrue(type(self.electrolyte2._cation == list[Cation]))
        self.assertTrue(type(self.electrolyte2._anion == list[Anion]))


class TestPolymer(unittest.TestCase):

    def setUp(self):
        self.ntype = NType(name='BBL', formal_reduction_potential=0.5)
        self.ptype = PType(name='PEDOT', formal_oxidation_potential=0.5)
        self.ntype2 = NType(name='bbl')
        self.ptype2 = PType(name='pedot')
        self.polymer = Polymer(name='P3HT')
        self.polymer2 = Polymer(name='p3ht')

    def test_attributes(self):
        """ Test the name attribute of the ntype class """
        self.assertTrue(self.ntype.name == 'BBL')
        self.assertTrue(self.ntype._name == 'bbl')
        self.assertTrue(self.ntype2.name == 'BBL')
        self.assertTrue(self.ntype2._name == 'bbl')
        self.assertTrue(self.ptype.name == 'PEDOT')
        self.assertTrue(self.ptype._name == 'pedot')
        self.assertTrue(self.ptype2.name == 'PEDOT')
        self.assertTrue(self.ptype2._name == 'pedot')
        self.assertTrue(self.polymer.name == 'P3HT')
        self.assertTrue(self.polymer._name == 'p3ht')
        self.assertTrue(self.polymer2.name == 'P3HT')
        self.assertTrue(self.polymer2._name == 'p3ht')
        self.ntype.formal_reduction_potential = 0.6
        self.assertTrue(self.ntype.formal_reduction_potential == 0.6)
        self.ptype.formal_oxidation_potential = 0.6
        self.assertTrue(self.ptype.formal_oxidation_potential == 0.6)
        