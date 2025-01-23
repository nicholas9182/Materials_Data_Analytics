import unittest
import tracemalloc
import pandas as pd
import plotly.express as px
import numpy as np
from Materials_Data_Analytics.quantum_chemistry.gaussian import GaussianParser
from Materials_Data_Analytics.core.coordinate_transformer import PdbParser
tracemalloc.start()


class TestPedotRaman(unittest.TestCase):

    def setUp(self):
        self.log = GaussianParser("./test_trajectories/pedot_raman/step1.log")

    def test_attributes(self):
        """ test the charge attributes of the parsers """
        self.assertTrue(self.log.charge == 0)
        self.assertTrue(self.log.multiplicity == 1)
        self.assertTrue(type(self.log.keywords) == list)
        self.assertTrue(self.log.raman is True)
        self.assertTrue(self.log.freq is True)
        self.assertTrue(self.log.opt is True)
        self.assertTrue(self.log.scf_iterations == 11)
        self.assertTrue(self.log.complete is True)
        self.assertTrue(self.log.atomcount == 28)
        self.assertTrue(self.log.energy == -4096904.424959145)
        self.assertTrue(self.log.functional == "B3LYP")
        self.assertTrue(self.log.basis == "6-311g")
        self.assertTrue(self.log.heavyatomcount == 18)
        self.assertTrue(self.log.esp is False)
        self.assertTrue(self.log.time_stamp == '2023-09-29 16:34:29')
        self.assertTrue(self.log.unrestricted is False)

    def test_raman(self):
        """ Test that the parser can extract the raman frequencies from the log file """
        raman_frequencies = self.log.get_raman_frequencies()
        self.assertTrue(type(raman_frequencies) == pd.DataFrame)
        self.assertTrue(len(raman_frequencies) == 59)

        raman_spectra = self.log.get_raman_spectra()
        self.assertTrue(len(raman_spectra) == 2000)

    def test_atoms(self):
        """ Test that the parser can extract the atoms from the log file """
        atoms = self.log.atoms
        self.assertTrue(atoms == ['C', 'C', 'C', 'S', 'C', 'O', 'O', 'C', 'C', 'H', 'H', 'H', 'H', 'H', 'C', 'C', 'C',
                                  'C', 'S', 'O', 'O', 'C', 'C', 'H', 'H', 'H', 'H', 'H'])
        
        heavyatoms = self.log.heavyatoms
        self.assertTrue(heavyatoms == ['C', 'C', 'C', 'S', 'C', 'O', 'O', 'C', 'C', 'C', 'C', 'C', 'C', 'S', 'O', 'O',
                                       'C', 'C'])
        
    def test_mulliken(self):
        """ Test that the parser can extract the mulliken charges from the log file """
        charges = self.log.get_mulliken_charges()
        self.assertTrue(charges['atom_id'].tolist() == [i for i in range(1, 29)])
        self.assertTrue(charges['element'].tolist() == self.log.atoms)

        charges = self.log.get_mulliken_charges(heavy_atoms=True)
        self.assertTrue(charges['element'].tolist() == self.log.heavyatoms)
        self.assertTrue(charges['partial_charge'].tolist() == [0.360925, 0.321741, -0.335369, 0.413119, -0.26947,
                                                               -0.505062, -0.511759, 0.258057, 0.267784, -0.335384,
                                                               0.321773, 0.360884, -0.26943, 0.413122, -0.511813,
                                                               -0.505008, 0.268105, 0.257783])

    def test_coordinates(self):
        """ Various tests for getting coordinates of the log file """
        coordinates = self.log.get_coordinates()
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

        coordinates = self.log.get_coordinates_through_scf()
        self.assertTrue([i for i in coordinates['iteration'].unique()] == [i for i in range(0, 11)])

        coordinates = self.log.get_coordinates(scf_iteration=0).round(4)
        self.assertTrue(coordinates['x'].to_list() == [-2.9799, -1.9997, -0.6323, -0.6799, -2.4332, -4.3287, -2.3374, 
                                                       -4.6132, -3.7137, -3.0014, -5.6637, -4.481, -3.817, -3.9767, 0.634, 
                                                       2.0013, 2.9849, 2.4372, 0.6819, 2.3293, 4.3294, 3.6881, 4.633, 
                                                       3.9619, 3.7353, 4.5628, 5.6681, 3.0007])
        self.assertTrue(coordinates['y'].to_list() == [-0.9153, 0.123, -0.3523, -2.135, -2.2084, -0.6969, 1.4271, 0.5998, 
                                                       1.6701, -3.1168, 0.8231, 0.5524, 1.6967, 2.658, 0.3524, -0.127, 0.9038, 
                                                       2.2153, 2.148, -1.4344, 0.6886, -1.5731, -0.7038, -2.6286, -1.3052,
                                                       -0.9441, -0.8594, 3.1337])
        self.assertTrue(coordinates['z'].to_list() == [-0.0515, 0.0122, -0.0714, -0.3392, -0.2825, 0.0092, 0.2637, 0.555, 
                                                       -0.0597, -0.4099, 0.3473, 1.6432, -1.1514, 0.3328, -0.0636, 0.0101, 
                                                       -0.0934, -0.164, -0.1985, 0.2539, 0.0407, 0.6915, -0.1345, 0.5925, 
                                                       1.7537, -1.2031, 0.1842, -0.117])
        
        coordinates = self.log.get_coordinates(heavy_atoms=True)
        self.assertTrue(len(coordinates) == self.log.heavyatomcount)

        data = self.log.get_mulliken_charges(with_coordinates=True)
        self.assertTrue(type(data) == pd.DataFrame)

        data = self.log.get_mulliken_charges(with_coordinates=True, heavy_atoms=True)
        self.assertTrue(type(data) == pd.DataFrame)


class TestBblLog3(unittest.TestCase):

    def setUp(self):
        self.log = GaussianParser("./test_trajectories/bbl/step3.log")

    def test_get_bonds(self):
        """ Test that the parser can extract bond information"""
        result = self.log.get_bonds_from_log()
        self.assertTrue(type(result))
        self.assertTrue(result.iloc[0, 0] == 1)
        self.assertTrue(result.iloc[10, 2] == 1.3935)
        self.assertTrue(result.iloc[30, 3] == "C")
        self.assertTrue(result.iloc[50, 1] == 41)

        result = self.log.get_bonds_from_coordinates()
        self.assertTrue(type(result))
        self.assertTrue(result.iloc[0, 0] == 1)
        self.assertTrue(result.iloc[10, 2] == 1.3935)
        self.assertTrue(result.iloc[30, 3] == "C")
        self.assertTrue(result.iloc[50, 1] == 41)

        result = self.log.get_bonds_from_coordinates(scf_iteration=0)
        self.assertTrue(type(result))
        self.assertTrue(result.iloc[0, 0] == 1)
        self.assertTrue(result.iloc[10, 2] == 1.3935)
        self.assertTrue(result.iloc[30, 3] == "C")
        self.assertTrue(result.iloc[50, 1] == 41)

    def test_get_scf_convergence(self):
        """ Test that the parser can extract the SCF convergence from the log file """
        result = self.log.get_scf_convergence().round(5)
        self.assertTrue(type(result) == pd.DataFrame)
        self.assertTrue(result['energy'].iloc[0] == -0.03369)
        self.assertTrue(result['cycles'].iloc[0] == 1)
        self.assertTrue(result['de'].iloc[0] == 0.01455)

    def test_attributes(self):
        """ test the charge attributes of the parsers """
        self.assertTrue(self.log.charge == -2)
        self.assertTrue(self.log.multiplicity == 3)
        self.assertTrue(self.log.raman is False)
        self.assertTrue(self.log.freq is False)
        self.assertTrue(self.log.atomcount == 140)
        self.assertTrue(self.log.energy == -12531769.403127551)
        self.assertTrue(self.log.functional == "WB97XD")
        self.assertTrue(self.log.basis == "6-311(d,p)")
        self.assertTrue(self.log.heavyatomcount == 110)
        self.assertTrue(self.log.esp is True)
        self.assertTrue(self.log.time_stamp == '2023-09-30 13:04:10')
        self.assertTrue(self.log.n_alpha == 363)
        self.assertTrue(self.log.n_beta == 361)
        self.assertTrue(self.log._n_electrons == 724)
        self.assertTrue(self.log.bandgap == 219.01921)
        self.assertTrue(self.log.homo == -224.979095)
        self.assertTrue(self.log.lumo == -5.959885)
        self.assertTrue(self.log.stable == 'untested')

    def test_mulliken(self):
        """ Test that the parser can extract the mulliken charges from the log file for bbl """
        charges = self.log.get_mulliken_charges()
        self.assertTrue(charges['atom_id'].tolist() == [i for i in range(1, 141)])
        self.assertTrue(charges['element'].tolist() == self.log.atoms)
        self.assertTrue(charges['partial_charge'].iloc[0] == -0.015446)
        self.assertTrue(charges['partial_charge'].iloc[1] == 0.124604)

        spins = self.log.get_mulliken_spin_densities()
        self.assertTrue(spins['atom_id'].tolist() == [i for i in range(1, 141)])
        self.assertTrue(spins['element'].tolist() == self.log.atoms)
        self.assertTrue(spins['spin_density'].iloc[0] == 0.031504)
        self.assertTrue(spins['spin_density'].iloc[1] == -0.001347)

        spins = self.log.get_mulliken_spin_densities(heavy_atoms=True)
        self.assertTrue(spins['element'].tolist() == self.log.heavyatoms)
        self.assertTrue(spins['spin_density'].iloc[0] == 0.030156)
        self.assertTrue(spins['spin_density'].iloc[1] == -0.012499)

    def test_esp(self):
        """ Test that the parser can extract the ESP partial charges from the log file for bbl """
        charges = self.log.get_esp_charges()
        self.assertTrue(type(charges) == pd.DataFrame)
        self.assertTrue(charges['partial_charge'].tolist()[0] == -0.043317)
        self.assertTrue(charges['partial_charge'].tolist()[5] == 0.625322)
        self.assertTrue(charges['partial_charge'].tolist()[10] == -0.112575)
        self.assertTrue(charges['partial_charge'].tolist()[15] == -0.142594)

        charges = self.log.get_esp_charges(heavy_atoms=True)
        self.assertTrue(type(charges) == pd.DataFrame)
        self.assertTrue(charges['element'].tolist() == self.log.heavyatoms)
        self.assertTrue(charges['partial_charge'].tolist()[0] == 0.044489)
        self.assertTrue(charges['partial_charge'].tolist()[5] == -0.284547)
        self.assertTrue(charges['partial_charge'].tolist()[10] == 0.049955)
        self.assertTrue(charges['partial_charge'].tolist()[15] == -0.080793)

    def test_get_coordinates(self):
        """ Test that the parser can extract the coordinates from the log file for heavy atoms for bbl """
        coordinates = self.log.get_coordinates(heavy_atoms=True)
        self.assertTrue(len(coordinates) == self.log.heavyatomcount)

        coordinates = self.log.get_coordinates()
        self.assertTrue(len(coordinates) == self.log.atomcount)
        self.assertTrue(coordinates['x'].iloc[0] == 1.817506)

    def test_get_spin_contamination(self):
        """ Test that the parser can extract the spin contamination from the log file """
        data = self.log.get_spin_contamination()
        self.assertTrue(type(data) == pd.DataFrame)
        self.assertTrue(data['before_annihilation'].iloc[0] == 2.0630)

    def test_get_orbitals(self):
        data = self.log.get_orbitals()
        self.assertTrue(type(data) == pd.DataFrame)
        self.assertTrue(data['energy'].iloc[5] == -50178.424725)
        figure = self.log.get_dos_plot()
        # figure.show()
        self.assertTrue(True)


class TestBblRaman(unittest.TestCase):

    def setUp(self):
        self.log = GaussianParser("./test_trajectories/bbl/raman.log")

    def test_attributes(self):
        """ test the charge attributes of the parsers """
        self.assertTrue(type(self.log.keywords) == list)
        self.assertTrue(self.log._raman is True)

    def test_get_thermochemistry(self):
        """ Test that the parser can extract the thermochemistry numbers from the log file """
        data = self.log.get_thermo_chemistry().round(5)
        self.assertTrue(type(data) == pd.DataFrame)
        self.assertTrue(data['g_corr'].iloc[0] == 3342.10134)
        self.assertTrue(data['e_elec_zp'].iloc[0] == -18887381.23105)
        self.assertTrue(data['g_elec_therm'].iloc[0] == -18887783.57055)


class TestBblLog4Restart(unittest.TestCase):

    def setUp(self):
        path1 = './test_trajectories/bbl/step4.log'
        path2 = './test_trajectories/bbl/step4_restart.log'
        self.log = GaussianParser([path1, path2])

    def test_attributes(self):
        self.assertTrue(type(self.log) == GaussianParser)
        self.assertTrue(self.log.energy == -18889245.002059143)
        self.assertTrue('opt' in self.log.keywords)
        self.assertTrue(self.log.raman is False)
        self.assertTrue(self.log.esp is False)
        self.assertTrue(self.log.complete is True)
        self.assertTrue(self.log.opt is True)
        self.assertTrue(self.log.functional == 'WB97XD')
        self.assertTrue(self.log.basis == '6-311(d,p)')
        self.assertTrue(self.log.restart == True)
        self.assertTrue(self.log.time_stamp == "2024-08-15 18:21:15")

