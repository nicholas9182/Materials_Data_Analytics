from Materials_Data_Analytics.experiment_modelling.cyclic_voltammetry import CyclicVoltammogram
import unittest
import pandas as pd
from Materials_Data_Analytics.materials.electrolytes import Electrolyte
from Materials_Data_Analytics.materials.ions import Cation, Anion  
from Materials_Data_Analytics.materials.solvents import Solvent
import plotly.express as px
import base64
import mimetypes
from plotly import graph_objects as go
from copy import copy


class TestCyclicVoltammetry(unittest.TestCase):

    na = Cation('Na+')
    cl = Anion('Cl-')
    water = Solvent('H2O')
    electrolyte = Electrolyte(cation=na, anion=cl, solvent=water, pH=7, temperature=298, concentrations={na: 1, cl: 1})

    data = pd.read_table('test_trajectories/cyclic_voltammetry/biologic1.txt', sep="\t")

    cv_biologic_1 = CyclicVoltammogram.from_biologic(path = 'test_trajectories/cyclic_voltammetry/biologic1.txt', electrolyte = electrolyte, metadata = {'scan_rate': 5, 'instrument': 'Biologic'})
    cv_biologic_2 = CyclicVoltammogram.from_biologic(data = data, electrolyte = electrolyte)
    cv_biologic_3 = CyclicVoltammogram.from_biologic(path = 'test_trajectories/cyclic_voltammetry/biologic2.txt', electrolyte = electrolyte)
    cv_biologic_4 = CyclicVoltammogram.from_biologic(path = 'test_trajectories/cyclic_voltammetry/biologic3.txt', electrolyte = electrolyte)
    cv_biologic_5 = CyclicVoltammogram.from_biologic(path = 'test_trajectories/cyclic_voltammetry/biologic4.txt', electrolyte = electrolyte)
    cv_biologic_6 = CyclicVoltammogram.from_biologic(path = 'test_trajectories/cyclic_voltammetry/biologic5.txt', electrolyte = electrolyte)
    cv_biologic_7 = CyclicVoltammogram.from_biologic(path = 'test_trajectories/cyclic_voltammetry/biologic6.txt', electrolyte = electrolyte)
    cv_aftermath_1 = CyclicVoltammogram.from_aftermath(path = 'test_trajectories/cyclic_voltammetry/aftermath1.csv', scan_rate=5)
    cv_aftermath_2 = CyclicVoltammogram.from_aftermath(path = 'test_trajectories/cyclic_voltammetry/aftermath2.csv', scan_rate=5)


    mime_type = mimetypes.guess_type('test_trajectories/cyclic_voltammetry/biologic1.txt')[0]
    if mime_type is None:
        mime_type = 'text/plain'

    with open('test_trajectories/cyclic_voltammetry/biologic1.txt', 'rb') as file:
        file_content = file.read()
        base64_data = base64.b64encode(file_content).decode('utf-8')
        base64_data = f'data:{mime_type};base64,{base64_data}'
        cv_html_base64 = CyclicVoltammogram.from_html_base64(file_contents = base64_data, electrolyte = electrolyte, source='biologic')

    
    def test_create_aftermath2(self):
        """ Test the creation from aftermath2 file """
        figure = self.cv_aftermath_2.get_current_time_plot()
        self.assertTrue(type(self.cv_aftermath_2) == CyclicVoltammogram)
        # figure.show()

    def test_from_biologic(self):
        """ Test the from_biologic method """
        self.assertTrue(type(self.cv_biologic_1) == CyclicVoltammogram)
        self.assertTrue(type(self.cv_biologic_1.data) == pd.DataFrame)
        self.assertTrue(self.cv_biologic_1.pH == 7)
        self.assertTrue(self.cv_biologic_1.temperature == 298)
        self.assertTrue(self.cv_biologic_1.cation == self.na)
        self.assertTrue(self.cv_biologic_1.anion == self.cl)
        self.assertTrue(self.cv_biologic_1.electrolyte == self.electrolyte)
        self.assertTrue('potential' in self.cv_biologic_1.data.columns) 
        self.assertTrue('current' in self.cv_biologic_1.data.columns)
        self.assertTrue('cycle' in self.cv_biologic_1.data.columns)
        self.assertTrue('time' in self.cv_biologic_1.data.columns)
        self.assertTrue('segment' in self.cv_biologic_1.data.columns)

    def test_from_biologic_dataframe(self):
        """ Test the from_biologic method with a dataframe """
        self.assertTrue(type(self.cv_biologic_2) == CyclicVoltammogram)
        self.assertTrue(type(self.cv_biologic_2.data) == pd.DataFrame)
        self.assertTrue(self.cv_biologic_2.pH == 7)
        self.assertTrue(self.cv_biologic_2.temperature == 298)
        self.assertTrue(self.cv_biologic_2.cation == self.na)
        self.assertTrue(self.cv_biologic_2.anion == self.cl)
        self.assertTrue(self.cv_biologic_2.electrolyte == self.electrolyte)
        self.assertTrue('potential' in self.cv_biologic_2.data.columns) 
        self.assertTrue('current' in self.cv_biologic_2.data.columns)
        self.assertTrue('cycle' in self.cv_biologic_2.data.columns)
        self.assertTrue('time' in self.cv_biologic_2.data.columns)

    def test_from_base64_biologic(self):
        """ Test the from_html_base64 method """
        self.assertTrue(type(self.cv_html_base64) == CyclicVoltammogram)
        self.assertTrue(type(self.cv_html_base64.data) == pd.DataFrame)
        self.assertTrue(self.cv_html_base64.pH == 7)
        self.assertTrue(self.cv_html_base64.temperature == 298)
        self.assertTrue(self.cv_html_base64.cation == self.na)
        self.assertTrue(self.cv_html_base64.anion == self.cl)
        self.assertTrue(self.cv_html_base64.electrolyte == self.electrolyte)
        self.assertTrue('potential' in self.cv_html_base64.data.columns) 
        self.assertTrue('current' in self.cv_html_base64.data.columns)
        self.assertTrue('cycle' in self.cv_html_base64.data.columns)
        self.assertTrue('time' in self.cv_html_base64.data.columns)

    def test_plot_1(self):
        """ Test the plot method """
        data = self.cv_biologic_1.data
        #px.line(data, x='potential', y='current', color='cycle', markers=True).show()
        self.assertTrue(type(self.cv_biologic_1.data == pd.DataFrame))

    def test_drop_cycles(self):
        """ Test the drop_cycles method """
        data = copy(self.cv_biologic_1).drop_cycles(drop=[1]).data
        # px.line(data, x='potential', y='current', color='cycle', markers=True).show()
        self.assertTrue(1 not in data['cycle'].values)

    def test_show_plots_biologic1(self):
        """ Test the show_current_potential, show_current_time and show_potential_time methods """
        # self.cv_biologic_1.show_current_potential()
        # self.cv_biologic_1.show_current_time()
        # self.cv_biologic_1.show_potential_time()
        self.assertTrue(type(self.cv_biologic_1.data == pd.DataFrame))

    def test_get_charge_passed(self):
        """ Test the get_charge_passed method """
        integrals = self.cv_biologic_1.get_charge_passed()
        charges = integrals.assign(anodic_charge = lambda x: x['anodic_charge']*1000).round(4)['anodic_charge'].to_list()
        self.assertTrue(type(integrals) == pd.DataFrame)
        self.assertTrue(charges == [0.004, 2.4893, 0.0040, 2.4983, 0.0040, 2.5017, 0.0039, 2.5035, 0.0039])

    def test_get_charge_passed_biologic2(self):
        """ Test the get_charge_passed method """
        integrals = self.cv_biologic_3.get_charge_passed()
        charges = integrals.assign(anodic_charge = lambda x: x['anodic_charge']*1000).round(4)['anodic_charge'].to_list()
        self.assertTrue(type(integrals) == pd.DataFrame)
        self.assertTrue(charges == [3.7438, 0.0048, 3.7419, 0.0067, 3.7393, 0.0065, 3.7381, 0.0048])

    def test_get_charge_passed_biologic3(self):
        """ Test the get_charge_passed method for biologic3 """
        integrals = self.cv_biologic_4.get_charge_passed()
        charges = integrals.assign(anodic_charge = lambda x: x['anodic_charge']*1000).round(4)['anodic_charge'].to_list()
        self.assertTrue(type(integrals) == pd.DataFrame)
        self.assertTrue(charges == [3.3026, 0.0065, 3.3097, 0.0065, 3.3090, 0.0065, 3.3072, 0.0066])

    def test_get_charge_passed_biologic4(self):
        """ Test the get_charge_passed method for biologic4 """
        integrals = self.cv_biologic_5.get_charge_passed()
        charges = integrals.assign(anodic_charge = lambda x: x['anodic_charge']*1000).round(4)['anodic_charge'].to_list()
        self.assertTrue(type(integrals) == pd.DataFrame)
        self.assertTrue(charges == [48.3689, 5232.7224, 112.5849])

    def test_get_charge_passed_biologic7(self):
        """ Test the get_charge_passed method for biologic7 """
        integrals = self.cv_biologic_7.get_charge_passed()
        charges = integrals.assign(anodic_charge = lambda x: x['anodic_charge']*1000).round(4)['anodic_charge'].to_list()
        self.assertTrue(type(integrals) == pd.DataFrame)
        self.assertTrue(charges == [0.0038, 6.0156, 0.003, 5.4851])

    def test_get_charge_passed_biologic4_av_segments(self):
        """ Test the get_charge_passed method for biologic4 with average segments """
        integrals = self.cv_biologic_1.get_charge_passed(average_segments=True)
        charges = integrals.assign(anodic_charge = lambda x: x['anodic_charge']*1000).round(4)['anodic_charge'].to_list()
        self.assertTrue(type(integrals) == pd.DataFrame)
        self.assertTrue(charges == [0.004, 2.4982])

    def test_show_charge_passed(self):
        """ Test the show_charge_passed method """
        # self.cv_biologic_1.show_charge_passed()
        self.assertTrue(type(self.cv_biologic_1.data == pd.DataFrame))

    def test_from_aftermath(self):
        """ Test the from_aftermath method """
        # self.cv_aftermath_1.show_current_potential()
        # self.cv_aftermath_1.show_current_time()
        # self.cv_aftermath_1.show_potential_time()
        self.assertTrue(type(self.cv_aftermath_1.data == pd.DataFrame))

    def test_make_cv_with_metadata(self):
        """ Test the make_cv_with_metadata method """
        data = self.cv_biologic_1.data
        self.assertTrue('scan_rate' in data.columns)
        self.assertTrue('instrument' in data.columns)

    def test_get_charge_integration_plot(self):
        """ Test the get_charge_integration_plot method """
        figure = self.cv_biologic_1.get_charge_integration_plot(cycle=3, direction='reduction')
        # figure.show()
        self.assertTrue(type(self.cv_biologic_1.data == pd.DataFrame))

    def test_get_charges(self):
        """ Test the get_charges method """
        charges = self.cv_biologic_1.get_charge_passed()
        self.assertTrue('total_charge' in charges.columns)
        self.assertTrue('anodic_charge' in charges.columns)
        self.assertTrue('cathodic_charge' in charges.columns)
        self.assertTrue(type(charges) == pd.DataFrame)
        self.assertTrue(all(charges['anodic_charge'] >= 0))
        self.assertTrue(all(charges['cathodic_charge'] >= 0))

    def test_get_charge_passed_biologic1_av_segments(self):
        """ Test the get_charge_passed method for biologic1 with average segments """
        integrals = self.cv_biologic_1.get_charge_passed(average_segments=True)
        self.assertTrue(type(integrals) == pd.DataFrame)
        self.assertTrue('anodic_charge' in integrals.columns)
        self.assertTrue('cathodic_charge' in integrals.columns)
        self.assertTrue('anodic_charge_err' in integrals.columns)
        self.assertTrue('cathodic_charge_err' in integrals.columns)
        self.assertTrue(all(integrals['anodic_charge'] >= 0))
        self.assertTrue(all(integrals['cathodic_charge'] >= 0))

    def test_get_maximum_charges_passed_biologic1(self):
        """ Test the get_maximum_charges_passed method for biologic1 """
        max_charges = self.cv_biologic_1.get_maximum_charges_passed()
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
        """ Test the get_maximum_charge_integration_plot method for anodic """
        figure = self.cv_biologic_6.get_maximum_charge_integration_plot(section=3)
        # figure.show()
        self.assertTrue(type(figure) == go.Figure)

    def test_downsample(self):
        """ Test the downsample method with biologic6 """
        cv = copy(self.cv_biologic_6).downsample(100)
        self.assertTrue(len(cv.data.query('segment == 3')) == 102)
        # cv.get_current_time_plot().show()
        # cv.get_potential_time_plot().show()

    def test_downsample_drop_and_downsample(self):
        """ Test the downsample and drop methods together with biologic6 """
        cv = copy(self.cv_biologic_6).drop_cycles(drop=[0, 1, 2]).downsample(100)
        self.assertTrue(len(cv.data.query('segment == 5')) == 101)
        # cv.get_current_time_plot().show()
        # cv.get_potential_time_plot().show()

    def test_get_peaks_biologic5(self):
        """ Test the get_peaks method for biologic5 """
        peaks = self.cv_biologic_6.get_peaks()
        self.assertTrue(type(peaks) == pd.DataFrame)
        self.assertTrue('current_peak' in peaks.columns)
        self.assertTrue('fit_current' in peaks.columns)
        self.assertTrue(all(peaks['current_peak'].notnull()))
        self.assertTrue(all(peaks['fit_current'].notnull()))

    def test_get_peak_plot_anodic(self):
        """ Test the get_peak_plot method for anodic """
        figure = self.cv_biologic_6.get_peak_plot(direction='oxidation', window = 0.02, width=700, height=500)
        self.assertTrue(type(figure) == go.Figure)
        self.assertTrue(len(figure.data) > 0)
        self.assertTrue(any(trace.name.startswith("Fitted") for trace in figure.data))
        self.assertTrue(any(trace.name.startswith("Peak") for trace in figure.data))
        # figure.show()

    def test_get_peak_plot_cathodic(self):
        """ Test the get_peak_plot method for cathodic for biologic6 """
        figure = self.cv_biologic_6.get_peak_plot(direction='reduction', window = 0.02, width=700, height=500)
        self.assertTrue(type(figure) == go.Figure)
        self.assertTrue(len(figure.data) > 0)
        self.assertTrue(any(trace.name.startswith("Fitted") for trace in figure.data))
        self.assertTrue(any(trace.name.startswith("Peak") for trace in figure.data))
        # figure.show()

    def test_get_peak_plot_cathodic_poly(self):
        """ Test the get_peak_plot method for cathodic for biologic6 with polynomial order """
        figure1 = self.cv_biologic_6.get_peak_plot(direction='oxidation', window = 0.02, polynomial_order=2, width=700, height=500)
        figure2 = self.cv_biologic_6.get_peak_plot(direction='oxidation', window = 0.02, polynomial_order=6, width=700, height=500)
        self.assertTrue(type(figure1) == go.Figure)
        self.assertTrue(len(figure1.data) > 0)
        self.assertTrue(len(figure2.data) > 0)
        self.assertTrue(any(trace.name.startswith("Fitted") for trace in figure1.data))
        self.assertTrue(any(trace.name.startswith("Peak") for trace in figure1.data))
        # figure1.show()
        # figure2.show()

    def test_get_plots_peaks_with_cycle(self):
        """ Test the get_plots_peaks_with_cycle method for biologic6 """
        current_figure, potential_figure = self.cv_biologic_6.get_plots_peaks_with_cycle(polynomial_order=4, window=0.02, width=700, height=500)
        self.assertTrue(type(current_figure) == go.Figure)
        self.assertTrue(type(potential_figure) == go.Figure)
        self.assertTrue(len(current_figure.data) > 0)
        self.assertTrue(len(potential_figure.data) > 0)
        self.assertTrue(any(trace.name.startswith("anodic peak") or trace.name.startswith("cathodic peak") for trace in current_figure.data))
        self.assertTrue(any(trace.name.startswith("anodic peak") or trace.name.startswith("cathodic peak") for trace in potential_figure.data))
        # current_figure.show()
        # potential_figure.show()

