from analytics.experiment_modelling.cyclic_voltammetry import CyclicVoltammogram
import unittest
import pandas as pd
from analytics.materials.electrolytes import Electrolyte
from analytics.materials.ions import Cation, Anion  
from analytics.materials.solvents import Solvent
import plotly.express as px
import base64
import mimetypes
from plotly import graph_objects as go


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

    def test_from_biologic_dataframe(self):

        data = pd.read_table('test_trajectories/cyclic_voltammetry/biologic1.txt', sep="\t")
        cv = CyclicVoltammogram.from_biologic(data = data, electrolyte = self.electrolyte)

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

    def test_from_base64_biologic(self):
            
        mime_type = mimetypes.guess_type('test_trajectories/cyclic_voltammetry/biologic1.txt')[0]
        if mime_type is None:
            mime_type = 'text/plain'

        with open('test_trajectories/cyclic_voltammetry/biologic1.txt', 'rb') as file:
            file_content = file.read()
            base64_data = base64.b64encode(file_content).decode('utf-8')
            base64_data = f'data:{mime_type};base64,{base64_data}'
            cv = CyclicVoltammogram.from_html_base64(file_contents = base64_data, electrolyte = self.electrolyte, source='biologic')
            
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

    def test_from_aftermath(self):

        cv = CyclicVoltammogram.from_aftermath(path = 'test_trajectories/cyclic_voltammetry/aftermath1.csv', scan_rate=5)

        # cv.show_current_potential()
        # cv.show_current_time()
        # cv.show_potential_time()
        self.assertTrue(type(cv.data == pd.DataFrame))

    def test_from_base64_aftermath(self):
            
        mime_type = mimetypes.guess_type('test_trajectories/cyclic_voltammetry/aftermath1.csv')[0]
        if mime_type is None:
            mime_type = 'text/plain'

        with open('test_trajectories/cyclic_voltammetry/aftermath1.csv', 'rb') as file:
            file_content = file.read()
            base64_data = base64.b64encode(file_content).decode('utf-8')
            base64_data = f'data:{mime_type};base64,{base64_data}'
            cv = CyclicVoltammogram.from_html_base64(file_contents = base64_data, electrolyte = self.electrolyte, source='aftermath', scan_rate=5)
            # cv.show_current_potential()
            # cv.show_current_time()
            # cv.show_potential_time()
            self.assertTrue(type(cv.data == pd.DataFrame))

    def test_make_cv_with_metadata(self):

        cv = CyclicVoltammogram.from_biologic(path = 'test_trajectories/cyclic_voltammetry/biologic1.txt', metadata = {'scan_rate': 5, 'instrument': 'Biologic'})
        data = cv.data
        self.assertTrue('scan_rate' in data.columns)
        self.assertTrue('instrument' in data.columns)

    def test_get_charge_integration_plot(self):

        cv = CyclicVoltammogram.from_biologic(path = 'test_trajectories/cyclic_voltammetry/biologic1.txt')
        figure = cv.get_charge_integration_plot(cycle=3, direction='reduction')
        # figure.show()
        self.assertTrue(type(cv.data == pd.DataFrame))

    def test_get_charges(self):
        
        cv = CyclicVoltammogram.from_biologic(path = 'test_trajectories/cyclic_voltammetry/biologic1.txt')
        charges = cv.get_charge_passed()
        self.assertTrue('total_charge' in charges.columns)
        self.assertTrue('anodic_charge' in charges.columns)
        self.assertTrue('cathodic_charge' in charges.columns)

    def test_get_charge_passed_biologic1(self):
        cv = CyclicVoltammogram.from_biologic(path='test_trajectories/cyclic_voltammetry/biologic1.txt', electrolyte=self.electrolyte)
        integrals = cv.get_charge_passed()
        self.assertTrue(type(integrals) == pd.DataFrame)
        self.assertTrue('anodic_charge' in integrals.columns)
        self.assertTrue('cathodic_charge' in integrals.columns)
        self.assertTrue('total_charge' in integrals.columns)
        self.assertTrue(all(integrals['anodic_charge'] >= 0))
        self.assertTrue(all(integrals['cathodic_charge'] >= 0))

    def test_get_charge_passed_biologic1_av_segments(self):
        cv = CyclicVoltammogram.from_biologic(path='test_trajectories/cyclic_voltammetry/biologic1.txt', electrolyte=self.electrolyte)
        integrals = cv.get_charge_passed(average_segments=True)
        self.assertTrue(type(integrals) == pd.DataFrame)
        self.assertTrue('anodic_charge' in integrals.columns)
        self.assertTrue('cathodic_charge' in integrals.columns)
        self.assertTrue('anodic_charge_err' in integrals.columns)
        self.assertTrue('cathodic_charge_err' in integrals.columns)
        self.assertTrue(all(integrals['anodic_charge'] >= 0))
        self.assertTrue(all(integrals['cathodic_charge'] >= 0))

    def test_get_maximum_charges_passed_biologic1(self):
        cv = CyclicVoltammogram.from_biologic(path='test_trajectories/cyclic_voltammetry/biologic1.txt', electrolyte=self.electrolyte)
        max_charges = cv.get_maximum_charges_passed()
        self.assertTrue(type(max_charges) == pd.DataFrame)
        self.assertTrue('total_charge' in max_charges.columns)
        self.assertTrue('section' in max_charges.columns)
        self.assertTrue('t_min' in max_charges.columns)
        self.assertTrue('t_max' in max_charges.columns)
        self.assertTrue('type' in max_charges.columns)
        self.assertTrue(all(max_charges['total_charge'] >= 0))
        self.assertTrue(set(max_charges['type']).issubset({'anodic_charge', 'cathodic_charge'}))
        self.assertTrue(max_charges.round(4).total_charge.to_list() == [0.0034, 0.0025, 0.0033, 0.0025, 0.0033, 0.0025, 0.0033, 0.0025])

    def test_get_maximum_charge_integration_plot_anodic(self):
        cv = CyclicVoltammogram.from_biologic(path='test_trajectories/cyclic_voltammetry/biologic5.txt', electrolyte=self.electrolyte)
        figure = cv.get_maximum_charge_integration_plot(section=3)
        # figure.show()
        self.assertTrue(type(figure) == go.Figure)

    def test_downsample(self):
        cv = CyclicVoltammogram.from_biologic(path='test_trajectories/cyclic_voltammetry/biologic5.txt', electrolyte=self.electrolyte)
        cv = cv.downsample(100)
        self.assertTrue(len(cv.data.query('segment == 3')) == 102)
        # cv.get_current_time_plot().show()
        # cv.get_potential_time_plot().show()

    def test_downsample_drop_and_downsample(self):

        cv = (CyclicVoltammogram
              .from_biologic(path='test_trajectories/cyclic_voltammetry/biologic5.txt')
              .drop_cycles(drop=[0, 1, 2])
              .downsample(100)
              )
        
        self.assertTrue(len(cv.data.query('segment == 5')) == 101)
        cv.get_current_time_plot().show()
        cv.get_potential_time_plot().show()
