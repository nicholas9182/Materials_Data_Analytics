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

    def test_from_benelegic(self):

        cv = CyclicVoltammogram.from_benelogic(path = 'test_trajectories/cyclic_voltammetry/benelogic1.txt', electrolyte = self.electrolyte)

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

        cv = CyclicVoltammogram.from_benelogic(path = 'test_trajectories/cyclic_voltammetry/benelogic1.txt', electrolyte = self.electrolyte)
        data = cv.data
        #px.line(data, x='potential', y='current', color='cycle', markers=True).show()
        self.assertTrue(type(cv.data == pd.DataFrame))

    def test_drop_cycles(self):

        data = (CyclicVoltammogram
                .from_benelogic(path = 'test_trajectories/cyclic_voltammetry/benelogic1.txt', electrolyte = self.electrolyte)
                .drop_cycles(cycles=[1])
                .data
                )

        # px.line(data, x='potential', y='current', color='cycle', markers=True).show()
        self.assertTrue(1 not in data['cycle'].values)

    def test_show_plots(self):

        cv = CyclicVoltammogram.from_benelogic(path = 'test_trajectories/cyclic_voltammetry/benelogic1.txt', electrolyte = self.electrolyte, keep_cycle_1=True)
        # cv.show_current_potential()
        cv.show_current_time()
        # cv.show_potential_time()
        self.assertTrue(type(cv.data == pd.DataFrame))

    def test_redox_direction(self):

        cv = CyclicVoltammogram.from_benelogic(path = 'test_trajectories/cyclic_voltammetry/benelogic1.txt', electrolyte = self.electrolyte)
        data = cv.data
        # px.line(data, x='time', y='current', color='redox', facet_col = 'cycle', markers=True, hover_data=['time']).show()
