import unittest
import tracemalloc
import pandas as pd
import plotly.express as px
import numpy as np
from Materials_Data_Analytics.quantum_chemistry.gaussian import GaussianParser
tracemalloc.start()


class TestGaussianParser(unittest.TestCase):
    """
    Test class for the GaussianParser class
    """
    pedot_log = GaussianParser("./test_trajectories/pedot_raman/step1.log")
    bbl_log_double_digit_charge = GaussianParser("./test_trajectories/bbl/step1.log")
    bbl_log = GaussianParser("./test_trajectories/bbl/step3.log")
    raman_log = GaussianParser("./test_trajectories/bbl/raman.log")
    bbl_spe_log = GaussianParser("./test_trajectories/bbl/step_37.log")
    bbl_step6_log = GaussianParser("./test_trajectories/bbl/step6.log")

    def test_multiline_keyword_parsing(self):
        """ Test that the parser can handle multiline keywords """
        result = self.bbl_step6_log.keywords
        self.assertTrue(result == ['scrf=(smd,solvent=water)', 'uwb97xd/6-311(d,p)', 'geom=allcheck', 'guess=read', 'stable'])

    def test_get_bonds_from_log(self):
        """ Test that the parser can extract bond information from the log file """
        result = self.bbl_log.get_bonds_from_log()
        self.assertTrue(type(result))
        self.assertTrue(result.iloc[0, 0] == 1)
        self.assertTrue(result.iloc[10, 2] == 1.3935)
        self.assertTrue(result.iloc[30, 3] == "C")
        self.assertTrue(result.iloc[50, 1] == 41)

    def test_get_bonds_from_coordinates(self):
        """ Test that the parser can extract bond information from the coordinates in the log file """
        result = self.bbl_log.get_bonds_from_coordinates()
        self.assertTrue(type(result))
        self.assertTrue(result.iloc[0, 0] == 1)
        self.assertTrue(result.iloc[10, 2] == 1.3935)
        self.assertTrue(result.iloc[30, 3] == "C")
        self.assertTrue(result.iloc[50, 1] == 41)

    def test_get_bonds_from_coordinates_pre_opt(self):
        """ Test that the parser can extract bond information from the coordinates in the log file before an optimisation """
        result = self.bbl_log.get_bonds_from_coordinates(pre_optimisation=True)
        self.assertTrue(type(result))
        self.assertTrue(result.iloc[0, 0] == 1)
        self.assertTrue(result.iloc[10, 2] == 1.3935)
        self.assertTrue(result.iloc[30, 3] == "C")
        self.assertTrue(result.iloc[50, 1] == 41)

    def test_two_bond_methods(self):
        """ Test that the two bond methods return the same results """
        result1 = self.bbl_log.get_bonds_from_log().sort_values(["atom_id_1", "atom_id_2"])
        result2 = self.bbl_log.get_bonds_from_coordinates().sort_values(["atom_id_1", "atom_id_2"])
        pd.testing.assert_frame_equal(result1.round(3), result2.round(3))

    def test_attributes(self):
        """ test the charge attributes of the parsers """
        self.assertTrue(self.pedot_log.charge == 0)
        self.assertTrue(self.bbl_log.charge == -2)
        self.assertTrue(self.bbl_log_double_digit_charge.charge == -12)
        self.assertTrue(self.pedot_log.multiplicity == 1)
        self.assertTrue(self.bbl_log.multiplicity == 3)
        self.assertTrue(type(self.pedot_log.keywords) == list)
        self.assertTrue(type(self.raman_log.keywords) == list)
        self.assertTrue(self.raman_log._raman is True)
        self.assertTrue(self.pedot_log.raman is True)
        self.assertTrue(self.bbl_log.raman is False)
        self.assertTrue(self.pedot_log.opt is True)
        self.assertTrue(self.pedot_log.complete is True)
        self.assertTrue(self.pedot_log.atomcount == 28)
        self.assertTrue(self.bbl_log.atomcount == 140)
        self.assertTrue(self.pedot_log.energy == -4096904.424959145)
        self.assertTrue(self.bbl_log.energy == -12531769.403127551)
        self.assertTrue(self.pedot_log.functional == "B3LYP")
        self.assertTrue(self.bbl_log.functional == "WB97XD")
        self.assertTrue(self.pedot_log.basis == "6-311g")
        self.assertTrue(self.bbl_log.basis == "6-311(d,p)")
        self.assertTrue(self.pedot_log.heavyatomcount == 18)
        self.assertTrue(self.bbl_log.heavyatomcount == 110)
        self.assertTrue(self.pedot_log.esp is False)
        self.assertTrue(self.bbl_log.esp is True)

    def test_raman_frequencies_pedot(self):
        """ Test that the parser can extract the raman frequencies from the log file """
        raman_frequencies = self.pedot_log.get_raman_frequencies()
        self.assertTrue(type(raman_frequencies) == pd.DataFrame)
        self.assertTrue(len(raman_frequencies) == 19)

    def test_raman_frequencies_bbl(self):
        """ Test that the parser can extract the raman frequencies from the log file for bbl """
        with self.assertRaises(ValueError):
            self.bbl_log.get_raman_frequencies()

    def test_raman_spectra(self):
        """ Test that the parser can extract the raman spectra from the log file for pedot """
        raman_spectra = self.pedot_log.get_raman_spectra()
        self.assertTrue(len(raman_spectra) == 2000)

    def test_energy_bbl_no_SCF(self):
        """ Test that the parser can extract the energy as None from a log file that did not converge """
        energy = GaussianParser("./test_trajectories/bbl/step1_no_scf.log").energy
        self.assertTrue(energy is None)

    def test_unrestricted(self):
        """ Test that the parser can extract the unrestricted attribute from the log file """
        unrestricted = self.pedot_log.unrestricted
        self.assertTrue(unrestricted is False)

    def test_atoms(self):
        """ Test that the parser can extract the atoms from the log file """
        atoms = self.pedot_log.atoms
        self.assertTrue(atoms == ['C', 'C', 'C', 'S', 'C', 'O', 'O', 'C', 'C', 'H', 'H', 'H', 'H', 'H', 'C', 'C', 'C',
                                  'C', 'S', 'O', 'O', 'C', 'C', 'H', 'H', 'H', 'H', 'H'])

    def test_heavy_atoms(self):
        """ Test that the parser can extract the heavy atoms from the log file """
        heavyatoms = self.pedot_log.heavyatoms
        self.assertTrue(heavyatoms == ['C', 'C', 'C', 'S', 'C', 'O', 'O', 'C', 'C', 'C', 'C', 'C', 'C', 'S', 'O', 'O',
                                       'C', 'C'])

    def test_get_mullikens_pedot(self):
        """ Test that the parser can extract the mulliken charges from the log file """
        charges = self.pedot_log.get_mulliken_charges()
        self.assertTrue(charges['atom_id'].tolist() == [i for i in range(1, 29)])
        self.assertTrue(charges['element'].tolist() == self.pedot_log.atoms)

    def test_get_mullikens_bbl(self):
        """ Test that the parser can extract the mulliken charges from the log file for bbl """
        charges = self.bbl_log.get_mulliken_charges()
        self.assertTrue(charges['atom_id'].tolist() == [i for i in range(1, 141)])
        self.assertTrue(charges['element'].tolist() == self.bbl_log.atoms)

    def test_get_mullikens_heavies(self):
        """ Test that the parser can extract the mulliken charges from the log file for heavy atoms """
        charges = self.pedot_log.get_mulliken_charges(heavy_atoms=True)
        self.assertTrue(charges['element'].tolist() == self.pedot_log.heavyatoms)
        self.assertTrue(charges['partial_charge'].tolist() == [0.360925, 0.321741, -0.335369, 0.413119, -0.26947,
                                                               -0.505062, -0.511759, 0.258057, 0.267784, -0.335384,
                                                               0.321773, 0.360884, -0.26943, 0.413122, -0.511813,
                                                               -0.505008, 0.268105, 0.257783])

    def test_esp_partials_pedot(self):
        """ Test that the parser can extract the ESP partial charges from the log file """
        with self.assertRaises(ValueError):
            self.pedot_log.get_esp_charges()

    def test_esp_partials_bbl(self):
        """ Test that the parser can extract the ESP partial charges from the log file for bbl """
        charges = self.bbl_log.get_esp_charges()
        self.assertTrue(type(charges) == pd.DataFrame)
        self.assertTrue(charges['partial_charge'].tolist()[0] == -0.043317)
        self.assertTrue(charges['partial_charge'].tolist()[5] == 0.625322)
        self.assertTrue(charges['partial_charge'].tolist()[10] == -0.112575)
        self.assertTrue(charges['partial_charge'].tolist()[15] == -0.142594)

    def test_esp_partials_bbl_heavies(self):
        """ Test that the parser can extract the ESP partial charges from the log file for bbl heavy atoms """
        charges = self.bbl_log.get_esp_charges(heavy_atoms=True)
        self.assertTrue(type(charges) == pd.DataFrame)
        self.assertTrue(charges['element'].tolist() == self.bbl_log.heavyatoms)
        self.assertTrue(charges['partial_charge'].tolist()[0] == 0.044489)
        self.assertTrue(charges['partial_charge'].tolist()[5] == -0.284547)
        self.assertTrue(charges['partial_charge'].tolist()[10] == 0.049955)
        self.assertTrue(charges['partial_charge'].tolist()[15] == -0.080793)

    def test_get_coordinates_pedot(self):
        """ Test that the parser can extract the coordinates from the log file """
        coordinates = self.pedot_log.get_coordinates()
        self.assertTrue(coordinates['x'].to_list() == [-2.896713, -1.906019, -0.602732, -0.582447, -2.389927,
                                                       -4.265714, -2.257699, -4.571961, -3.656809, -2.920781,
                                                       -5.613404, -4.452347, -3.755911, -3.846831, 0.602789,
                                                       1.906062, 2.896732, 2.389966, 0.582486, 2.257544,
                                                       4.265750, 3.657369, 4.571336, 3.846797, 3.758349,
                                                       4.449836, 5.613228, 2.920821])
        self.assertTrue(coordinates['y'].to_list() == [-0.992769, 0.044153, -0.383270, -2.223480, -2.249765,
                                                       -0.700007, 1.396997, 0.678906, 1.651647, -3.182215,
                                                       0.831367, 0.768250, 1.546480, 2.680598, 0.383293,
                                                       -0.044139, 0.992787, 2.249771, 2.223486, -1.396968,
                                                       0.699883, -1.652235, -0.678299, -2.680703, -1.548920,
                                                       -0.765647, -0.831226, 3.182215])
        self.assertTrue(coordinates['z'].to_list() == [-0.031111, -0.012812, -0.011038, -0.044588, -0.049912,
                                                       -0.050345, 0.022211, 0.360569, -0.359029, -0.067273,
                                                       0.091302, 1.442186, -1.440908, -0.067731, 0.012085,
                                                       0.013091, 0.031796, 0.051678, 0.047337, -0.022861,
                                                       0.049976, 0.355428, -0.364014, 0.062057, 1.437312,
                                                       -1.445597, -0.096827, 0.069261])

    def test_get_coordinates_pedot_pre_eq(self):
        """ Test that the parser can extract the coordinates from the log file before an optimisation for pedot """
        coordinates = self.pedot_log.get_coordinates(pre_optimisation=True).round(3)
        self.assertTrue(coordinates['x'].to_list()[0] == -2.98)
        self.assertTrue(coordinates['x'].to_list()[1] == -2.00)
        self.assertTrue(coordinates['x'].to_list()[2] == -0.632)

    def test_get_coordinates_pedot_heavies(self):
        """ Test that the parser can extract the coordinates from the log file for heavy atoms for pedot """
        coordinates = self.pedot_log.get_coordinates(heavy_atoms=True)
        self.assertTrue(len(coordinates) == self.pedot_log.heavyatomcount)

    def test_get_coordinates_bbl_heavies(self):
        """ Test that the parser can extract the coordinates from the log file for heavy atoms for bbl """
        coordinates = self.bbl_log.get_coordinates(heavy_atoms=True)
        self.assertTrue(len(coordinates) == self.bbl_log.heavyatomcount)

    def test_get_coordinates_bbl(self):
        """ Test that the parser can extract the coordinates from the log file for bbl """
        coordinates = self.bbl_log.get_coordinates()
        self.assertTrue(len(coordinates) == self.bbl_log.atomcount)
        self.assertTrue(coordinates['x'].iloc[0] == 1.817506)

    def test_get_mulliken_charges_with_coordinates_pedot(self):
        """ Test that the parser can extract the mulliken charges with coordinates from the log file for pedot """
        data = self.pedot_log.get_mulliken_charges(with_coordinates=True)
        self.assertTrue(type(data) == pd.DataFrame)

    def test_get_mulliken_charges_with_coordinates_bbl(self):
        """ Test that the parser can extract the mulliken charges with coordinates from the log file for bbl """
        data = self.pedot_log.get_mulliken_charges(with_coordinates=True, heavy_atoms=True)
        self.assertTrue(type(data) == pd.DataFrame)

    def test_get_spin_contamination(self):
        """ Test that the parser can extract the spin contamination from the log file """
        data = self.bbl_log.get_spin_contamination()
        self.assertTrue(type(data) == pd.DataFrame)
        self.assertTrue(data['before_annihilation'].iloc[0] == 2.0630)

    def test_get_spin_contamination_rDFT(self):
        """ Test that the parser can extract the spin contamination from the log file for rDFT """
        data = self.raman_log.get_spin_contamination()
        self.assertTrue(type(data) == pd.DataFrame)

    def test_get_thermochemistry(self):
        """ Test that the parser can extract the thermochemistry numbers from the log file """
        data = self.raman_log.get_thermo_chemistry()
        self.assertTrue(type(data) == pd.DataFrame)
        self.assertTrue(data['g_corr'].iloc[0] == 3342.10134)
        self.assertTrue(data['e_elec_zp'].iloc[0] == -18887381.23105)
        self.assertTrue(data['g_elec_therm'].iloc[0] == -18887783.57055)


class TestGaussianParserRestart(unittest.TestCase):

    def test_is_parser(self):
        with self.assertRaises(ValueError):
            bbl_log = GaussianParser("./test_trajectories/bbl/step5_raman_restart.log")


class TestGaussianParserRestart2(unittest.TestCase):
    """ Unittest class to test gaussian parsers with restart files """

    path1 = './test_trajectories/bbl/step4.log'
    path2 = './test_trajectories/bbl/step4_restart.log'
    bbl_log = GaussianParser([path1, path2])
    bbl_log_tuple = GaussianParser((path1, path2))

    def test_parser_restart(self):
        self.assertTrue(type(self.bbl_log) == GaussianParser)
        self.assertTrue(self.bbl_log.energy == -18889245.002059143)
        self.assertTrue('opt' in self.bbl_log.keywords)
        self.assertTrue(self.bbl_log.raman is False)
        self.assertTrue(self.bbl_log.esp is False)
        self.assertTrue(self.bbl_log.complete is True)
        self.assertTrue(self.bbl_log.opt is True)
        self.assertTrue(self.bbl_log.functional == 'WB97XD')
        self.assertTrue(self.bbl_log.basis == '6-311(d,p)')
        self.assertTrue(self.bbl_log.restart == True)

    def test_parser_restart_tuple(self):
        self.assertTrue(type(self.bbl_log_tuple) == GaussianParser)
        self.assertTrue(self.bbl_log_tuple.energy == -18889245.002059143)
        self.assertTrue('opt' in self.bbl_log_tuple.keywords)
        self.assertTrue(self.bbl_log_tuple.raman is False)
        self.assertTrue(self.bbl_log_tuple.esp is False)
        self.assertTrue(self.bbl_log_tuple.complete is True)
        self.assertTrue(self.bbl_log_tuple.opt is True)
        self.assertTrue(self.bbl_log_tuple.functional == 'WB97XD')
        self.assertTrue(self.bbl_log_tuple.basis == '6-311(d,p)')
        self.assertTrue(self.bbl_log_tuple.restart == True)


class TestGaussianInstabilities(unittest.TestCase):

    stable_log = GaussianParser("./test_trajectories/bbl/step2.log")
    not_tested_log = GaussianParser("./test_trajectories/bbl/step3.log")
    rhf_unstable = GaussianParser("./test_trajectories/bbl/rhf_instability.log")
    internal_log = GaussianParser("./test_trajectories/bbl/internal_instability.log")

    def test_stable(self):
        report = self.stable_log.stable
        self.assertTrue(report == 'stable')

    def test_not_tested(self):
        report = self.not_tested_log.stable
        self.assertTrue(report == 'untested')

    def test_rhf(self):
        report = self.rhf_unstable.stable
        self.assertTrue(report == 'RHF instability')

    def test_internal(self):
        report = self.internal_log.stable
        self.assertTrue(report == 'internal instability')
