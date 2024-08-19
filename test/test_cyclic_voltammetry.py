from analytics.experiment_modelling.cyclic_voltammetry import CyclicVoltammogram
import unittest
import pandas as pd
from analytics.materials.electrolytes import Electrolyte
from analytics.materials.ions import Cation, Anion  
from analytics.materials.solvents import Solvent
import plotly.express as px


class TestCyclicVoltammetry(unittest.TestCase):

    na = Cation('Na+')
    cl = Anion('Cl-')
    water = Solvent('H2O')
    electrolyte = Electrolyte(cation=na, anion=cl, solvent=water, pH=7, temperature=298, concentrations={na: 1, cl: 1})

    def test_from_biologic(self):

        cv = CyclicVoltammogram.from_biologic(path = 'test_trajectories/cyclic_voltammetry/biologic1.txt', electrolyte = self.electrolyte)

        self.assertTrue(type(cv) == CyclicVoltammogram)
        self.assertTrue(type(cv.data) == pd.DataFrame)
        self.assertTrue(cv.pH == 7)
        self.assertTrue(cv.temperature == 298)
        self.assertTrue(cv.cation == self.na)
        self.assertTrue(cv.anion == self.cl)
        self.assertTrue(cv.electrolyte == self.electrolyte)
        self.assertTrue('potential' in cv.data.columns) 
        self.assertTrue('current' in cv.data.columns)
        self.assertTrue('cycle' in cv.data.columns)
        self.assertTrue('time' in cv.data.columns)

    def test_plot_1(self):

        cv = CyclicVoltammogram.from_biologic(path = 'test_trajectories/cyclic_voltammetry/biologic1.txt', electrolyte = self.electrolyte)
        data = cv.data
        #px.line(data, x='potential', y='current', color='cycle', markers=True).show()
        self.assertTrue(type(cv.data == pd.DataFrame))

    def test_drop_cycles(self):

        data = (CyclicVoltammogram
                .from_biologic(path = 'test_trajectories/cyclic_voltammetry/biologic1.txt', electrolyte = self.electrolyte)
                .drop_cycles(drop=[1])
                .data
                )

        # px.line(data, x='potential', y='current', color='cycle', markers=True).show()
        self.assertTrue(1 not in data['cycle'].values)

    def test_show_plots_biologic1(self):

        cv = CyclicVoltammogram.from_biologic(path = 'test_trajectories/cyclic_voltammetry/biologic1.txt', electrolyte = self.electrolyte)
        # cv.show_current_potential()
        # cv.show_current_time()
        # cv.show_potential_time()
        self.assertTrue(type(cv.data == pd.DataFrame))

    def test_show_plots_biologic2(self):

        cv = CyclicVoltammogram.from_biologic(path = 'test_trajectories/cyclic_voltammetry/biologic2.txt', electrolyte = self.electrolyte)
        # cv.show_current_potential()
        # cv.show_current_time()
        # cv.show_potential_time()
        self.assertTrue(type(cv.data == pd.DataFrame))

    def test_show_plots_biologic3(self):

        cv = CyclicVoltammogram.from_biologic(path = 'test_trajectories/cyclic_voltammetry/biologic3.txt', electrolyte = self.electrolyte)
        # cv.show_current_potential()
        # cv.show_current_time()
        # cv.show_potential_time()
        self.assertTrue(type(cv.data == pd.DataFrame))

    def test_show_plots_biologic4(self):

        cv = CyclicVoltammogram.from_biologic(path = 'test_trajectories/cyclic_voltammetry/biologic4.txt', electrolyte = self.electrolyte)
        # cv.show_current_potential()
        # cv.show_current_time()
        # cv.show_potential_time()
        self.assertTrue(type(cv.data == pd.DataFrame))

    def test_redox_direction(self):

        cv = CyclicVoltammogram.from_biologic(path = 'test_trajectories/cyclic_voltammetry/biologic1.txt', electrolyte = self.electrolyte)
        data = cv.data
        # px.line(data, x='time', y='current', color='redox', facet_col = 'cycle', markers=True, hover_data=['time']).show()

    def test_get_charge_passed(self):

        cv = CyclicVoltammogram.from_biologic(path = 'test_trajectories/cyclic_voltammetry/biologic1.txt', electrolyte = self.electrolyte)
        integrals = cv.get_charge_passed()
        charges = integrals.assign(anodic_charge = lambda x: x['anodic_charge']*1000).round(4)['anodic_charge'].to_list()
        self.assertTrue(type(integrals) == pd.DataFrame)
        self.assertTrue(charges == [0.004, 2.4893, 0.0040, 2.4983, 0.0040, 2.5017, 0.0039, 2.5035, 0.0039])

    def test_get_charge_passed_biologic2(self):

        cv = CyclicVoltammogram.from_biologic(path = 'test_trajectories/cyclic_voltammetry/biologic2.txt', electrolyte = self.electrolyte)
        integrals = cv.get_charge_passed()
        charges = integrals.assign(anodic_charge = lambda x: x['anodic_charge']*1000).round(4)['anodic_charge'].to_list()
        self.assertTrue(type(integrals) == pd.DataFrame)
        self.assertTrue(charges == [0.0042, 3.7438, 0.0048, 3.7424, 0.0049, 3.7397, 0.0048, 3.7381, 0.0048])

    def test_get_charge_passed_biologic3(self):

        cv = CyclicVoltammogram.from_biologic(path = 'test_trajectories/cyclic_voltammetry/biologic3.txt', electrolyte = self.electrolyte)
        integrals = cv.get_charge_passed()
        charges = integrals.assign(anodic_charge = lambda x: x['anodic_charge']*1000).round(4)['anodic_charge'].to_list()
        self.assertTrue(type(integrals) == pd.DataFrame)
        self.assertTrue(charges == [3.3026, 0.0065, 3.3097, 0.0065, 3.3090, 0.0065, 3.3072, 0.0066])

    def test_get_charge_passed_biologic4(self):

        cv = CyclicVoltammogram.from_biologic(path = 'test_trajectories/cyclic_voltammetry/biologic4.txt', electrolyte = self.electrolyte)
        integrals = cv.get_charge_passed()
        charges = integrals.assign(anodic_charge = lambda x: x['anodic_charge']*1000).round(4)['anodic_charge'].to_list()
        self.assertTrue(type(integrals) == pd.DataFrame)
        self.assertTrue(charges == [48.3689, 5232.7224, 112.5849])

    def test_get_charge_passed_biologic4_av_segments(self):

        cv = CyclicVoltammogram.from_biologic(path = 'test_trajectories/cyclic_voltammetry/biologic1.txt', electrolyte = self.electrolyte)
        integrals = cv.get_charge_passed(average_segments=True)
        charges = integrals.assign(anodic_charge = lambda x: x['anodic_charge']*1000).round(4)['anodic_charge'].to_list()
        self.assertTrue(type(integrals) == pd.DataFrame)
        self.assertTrue(charges == [0.004, 2.4982])

    def test_show_charge_passed(self):

        cv = CyclicVoltammogram.from_biologic(path = 'test_trajectories/cyclic_voltammetry/biologic1.txt', electrolyte = self.electrolyte)
        # cv.show_charge_passed()
        self.assertTrue(type(cv.data == pd.DataFrame))
