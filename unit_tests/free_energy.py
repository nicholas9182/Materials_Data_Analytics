import unittest
import tracemalloc
import os
import plotly.graph_objects as go
from glob import glob
import pandas as pd
import plumed as pl
from analytics.metadynamics.free_energy import FreeEnergySpace, MetaTrajectory, FreeEnergyLine
tracemalloc.start()


class TestMetaTrajectory(unittest.TestCase):

    def test_colvar_read(self):
        """
        checking that the MetaTrajectory is reading in and processing colvar files correctly. Comparing with a direct plumed read in
        """
        file = "../test_trajectories/ndi_na_binding/COLVAR.0"
        cv_traj = MetaTrajectory(file)
        self.assertEqual(cv_traj.data.columns.to_list(), ['time', 'D1', 'CM1', 'bias', 'reweight_bias', 'reweight_factor'])
        self.assertEqual(cv_traj.walker, 0)
        self.assertEqual(cv_traj.cvs, ['D1', 'CM1'])


class TestFreeEnergyLine(unittest.TestCase):

    def test_fes_read(self):
        """
        checking that the 1d fes file is being read in correctly and the cv extracted correctly. Comparing with a direct plumed read in
        """
        file = "../test_trajectories/ndi_na_binding/FES_CM1.dat"
        line = FreeEnergyLine(file)
        compare = pd.DataFrame(pl.read_as_pandas(file))
        self.assertEqual(line.data.loc[:, 'energy'].to_list(), compare.loc[:, 'projection'].to_list())
        self.assertEqual(line.cv, 'CM1')

    def test_fes_read_with_time_data(self):
        """
        checking that alternate constructor works for reading in fes data with strides to get time_data dictionary
        """
        folder = "../test_trajectories/ndi_na_binding/FES_CM1/"
        pattern = "FES*dat"
        all_fes_files = [file for folder, subdir, files in os.walk(folder) for file in glob(os.path.join(folder, pattern))]
        line = FreeEnergyLine(all_fes_files)
        compare1 = pd.DataFrame(pl.read_as_pandas("../test_trajectories/ndi_na_binding/FES_CM1/FES20.dat"))
        compare2 = pd.DataFrame(pl.read_as_pandas("../test_trajectories/ndi_na_binding/FES_CM1/FES23.dat"))
        self.assertEqual(line.time_data[20].loc[:, 'energy'].to_list(), compare1.loc[:, 'projection'].to_list())
        self.assertEqual(line.time_data[23].loc[:, 'energy'].to_list(), compare2.loc[:, 'projection'].to_list())
        self.assertEqual(line.data.loc[:, 'energy'].to_list(), compare2.loc[:, 'projection'].to_list())

    def test_normalise_with_float(self):
        """
        testing that the normalise function works with a single value
        :return:
        """
        file = "../test_trajectories/ndi_na_binding/FES_CM1.dat"
        line = FreeEnergyLine(file)
        line.set_datum(0)
        figure = go.Figure()
        trace = go.Scatter(x=line.data[line.cv], y=line.data['energy'])
        figure.add_trace(trace)
        self.assertTrue(0 in line.data['energy'])
        # figure.show()

    def test_normalise_with_tuple(self):
        """
        testing that the normalise function works with a range
        :return:
        """
        file = "../test_trajectories/ndi_na_binding/FES_CM1.dat"
        line = FreeEnergyLine(file)
        line.set_datum(datum=(6, 8))
        self.assertAlmostEqual(line.data.loc[line.data['CM1'] > 6].loc[line.data['CM1'] < 8]['energy'].mean(), 0)
        figure = go.Figure()
        trace = go.Scatter(x=line.data[line.cv], y=line.data['energy'])
        figure.add_trace(trace)
        # figure.show()

    def test_normalise_with_tuple_on_time_data(self):
        """
        testing that the normalise function works with a range
        :return:
        """
        folder = "../test_trajectories/ndi_na_binding/FES_CM1/"
        pattern = "FES*dat"
        all_fes_files = [file for folder, subdir, files in os.walk(folder) for file in glob(os.path.join(folder, pattern))]
        line = FreeEnergyLine(all_fes_files)
        line.set_datum(datum=(6, 8))
        self.assertAlmostEqual(line.data.loc[line.data['CM1'] > 6].loc[line.data['CM1'] < 8]['energy'].mean(), 0)
        figure = go.Figure()
        trace = go.Scatter(x=line.data[line.cv], y=line.data['energy'])
        figure.add_trace(trace)
        # figure.show()

    def test_get_change_over_time(self):
        """
        testing that the normalise function works with a range
        :return:
        """
        folder = "../test_trajectories/ndi_na_binding/FES_CM1/"
        pattern = "FES*dat"
        all_fes_files = [file for folder, subdir, files in os.walk(folder) for file in glob(os.path.join(folder, pattern))]
        line = FreeEnergyLine(all_fes_files)
        change_data = line.get_time_difference(1, 3)
        figure = go.Figure()
        trace = go.Scatter(x=change_data['time_stamp'], y=change_data['energy_difference'])
        figure.add_trace(trace)
        # figure.show()

    def test_get_change_over_time_using_tuples(self):
        """
        testing that the normalise function works with a range
        :return:
        """
        folder = "../test_trajectories/ndi_na_binding/FES_CM1/"
        pattern = "FES*dat"
        all_fes_files = [file for folder, subdir, files in os.walk(folder) for file in glob(os.path.join(folder, pattern))]
        line = FreeEnergyLine(all_fes_files)
        change_data = line.get_time_difference(region_1=(0.8, 1.2), region_2=(2.8, 3.2))
        figure = go.Figure()
        trace = go.Scatter(x=change_data['time_stamp'], y=change_data['energy_difference'])
        figure.add_trace(trace)
        # figure.show()

    def test_set_datum_twice(self):
        """
        testing that the normalise function works with a single value
        :return:
        """
        folder = "../test_trajectories/ndi_na_binding/FES_CM1/"
        pattern = "FES*dat"
        all_fes_files = [file for folder, subdir, files in os.walk(folder) for file in glob(os.path.join(folder, pattern))]
        line = FreeEnergyLine(all_fes_files)
        data1 = line.set_datum(3).data
        data2 = line.set_datum(3).data
        pd.testing.assert_frame_equal(data2, data1)
        figure = go.Figure()
        trace = go.Scatter(x=line.data[line.cv], y=line.data['energy'])
        figure.add_trace(trace)
        self.assertTrue(0 in line.data['energy'])
        # figure.show()


class TestFreeEnergyLandscape(unittest.TestCase):

    def test_make_landscape(self):
        """
        check that landscape constructor works
        :return:
        """
        landscape = FreeEnergySpace("../test_trajectories/ndi_na_binding/HILLS")
        self.assertTrue('height' in landscape.hills.columns.to_list())
        self.assertTrue('time' in landscape.hills.columns.to_list())
        self.assertEqual(type(landscape), FreeEnergySpace)
        self.assertEqual(landscape.cvs, ['D1', 'CM1'])
        self.assertEqual(landscape.n_walker, 8)
        self.assertEqual(landscape.n_timesteps, 2979)

    def test_hills_plotter_default_values(self):

        landscape = FreeEnergySpace("../test_trajectories/ndi_na_binding/HILLS")
        figures = landscape.get_hills_figures()
        self.assertEqual(len(figures), 8)
        self.assertTrue(figures[0]._validate)
        self.assertTrue(figures[1]._validate)
        self.assertTrue(figures[2]._validate)
        self.assertTrue(figures[3]._validate)

    def test_fes_adder_checks_work(self):

        landscape = FreeEnergySpace("../test_trajectories/ndi_na_binding/HILLS")
        fes = FreeEnergyLine("../test_trajectories/ndi_na_binding/FES_CM1.dat")
        landscape.add_line(fes)
        landscape.add_line(fes)
        self.assertEqual(landscape.lines['CM1'], fes)

    def test_traj_adder_checks_work(self):

        landscape = FreeEnergySpace("../test_trajectories/ndi_na_binding/HILLS")
        traj = MetaTrajectory("../test_trajectories/ndi_na_binding/COLVAR.0")
        landscape.add_metad_trajectory(traj)
        landscape.add_metad_trajectory(traj)
        self.assertEqual(landscape.trajectories[0], traj)
