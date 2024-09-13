import unittest
import tracemalloc
import os
import plotly.graph_objects as go
from glob import glob
import pandas as pd
import plumed as pl
import matplotlib.pyplot as plt
from analytics.metadynamics.free_energy import FreeEnergySpace, MetaTrajectory, FreeEnergyLine, FreeEnergySurface
tracemalloc.start()


class TestMetaTrajectory(unittest.TestCase):

    def test_colvar_read(self):
        """
        checking that the MetaTrajectory is reading in and processing colvar files correctly. Comparing with a direct
        plumed read in
        """
        file = "./test_trajectories/ndi_na_binding/COLVAR.0"
        cv_traj = MetaTrajectory(file)
        self.assertEqual(cv_traj._data.columns.to_list(), ['time', 'D1', 'CM1', 'bias', 'reweight_bias', 'reweight_factor', 'weight'])
        self.assertEqual(cv_traj.walker, 0)
        self.assertEqual(cv_traj.cvs, ['D1', 'CM1'])
        self.assertTrue(cv_traj._opes is False)

    def test_colvar_read_opes(self):
        """
        checking that the MetaTrajectory is reading in and processing colvar files correctly.
        Comparing with a direct plumed read in
        """
        file = "./test_trajectories/ndi_single_opes/COLVAR.0"
        cv_traj = MetaTrajectory(file)
        self.assertEqual(cv_traj._data.columns.to_list(), ['time', 'D1', 'CM1', 'reweight_bias','reweight_factor', 'zed', 'neff', 'nker', 'weight'])
        self.assertEqual(cv_traj.walker, 0)
        self.assertEqual(cv_traj.cvs, ['D1', 'CM1'])
        self.assertTrue(cv_traj._opes is True)


class TestFreeEnergyLine(unittest.TestCase):

    def test_fes_read(self):
        """
        checking that the 1d fes file is being read in correctly and the cv extracted correctly.
        Comparing with a direct plumed read in
        """
        file = "./test_trajectories/ndi_na_binding/FES_CM1.dat"
        pan_file = pl.read_as_pandas(file)
        file = pd.DataFrame(pan_file).rename(columns={'projection': 'energy'})
        line = FreeEnergyLine(file)
        self.assertEqual(line._data.loc[:, 'energy'].to_list(), file.loc[:, 'energy'].to_list())
        self.assertEqual(line.cvs[0], 'CM1')

    def test_fes_read_from_plumed(self):
        """
        checking that the 1d fes file is being read in correctly and the cv extracted correctly.
        Comparing with a direct plumed read in
        """
        file = "./test_trajectories/ndi_na_binding/FES_CM1.dat"
        line = FreeEnergyLine.from_plumed(file)
        self.assertEqual(line.cvs[0], 'CM1')

    def test_fes_read_with_time_data(self):
        """
        checking that alternate constructor works for reading in fes _data with strides to get _time_data dictionary
        """
        folder = "./test_trajectories/ndi_na_binding/FES_CM1/"
        pattern = "FES*dat"

        all_fes_files = [file for folder, subdir, files in os.walk(folder)
                         for file in glob(os.path.join(folder, pattern))]
        
        individual_files = [f.split("/")[-1] for f in all_fes_files]
        time_stamps = [int(''.join(x for x in f if x.isdigit())) for f in individual_files]
        data_frames = [FreeEnergyLine._read_file(f) for f in all_fes_files]
        data = {time_stamps[i]: data_frames[i] for i in range(0, len(time_stamps))}
        line = FreeEnergyLine(data)

        compare = pd.DataFrame(pl.read_as_pandas("./test_trajectories/ndi_na_binding/FES_CM1/FES2.dat"))
        energy_diff = compare.loc[1, 'projection'] - compare.loc[2, 'projection']
        my_diff = line._time_data[2].loc[1, 'energy'] - line._time_data[2].loc[2, 'energy']
        self.assertEqual(energy_diff, my_diff)

    def test_fes_read_with_time_data_from_plumed(self):
        """
        checking that alternate constructor works for reading in fes _data with strides to get _time_data dictionary
        """
        folder = "./test_trajectories/ndi_na_binding/FES_CM1/"
        pattern = "FES*dat"
        all_fes_files = [file for folder, subdir, files in os.walk(folder) for file in glob(os.path.join(folder, pattern))]
        line = FreeEnergyLine.from_plumed(all_fes_files)
        compare = pd.DataFrame(pl.read_as_pandas("./test_trajectories/ndi_na_binding/FES_CM1/FES2.dat"))
        energy_diff = compare.loc[1, 'projection'] - compare.loc[2, 'projection']
        my_diff = line._time_data[2].loc[1, 'energy'] - line._time_data[2].loc[2, 'energy']
        self.assertEqual(energy_diff, my_diff)

    def test_normalise_with_float(self):
        """
        testing that the normalise function works with a single value
        :return:
        """
        file = "./test_trajectories/ndi_na_binding/FES_CM1.dat"
        line = FreeEnergyLine.from_plumed(file)
        line.set_datum({'CM1': 0})
        figure = go.Figure()
        trace = go.Scatter(x=line._data[line.cvs[0]], y=line._data['energy'])
        figure.add_trace(trace)
        self.assertTrue(0 in line._data['energy'])
        # figure.show()

    def test_normalise_with_float_on_time_data(self):
        """
        testing that the normalise function works with a range
        :return:
        """
        folder = "./test_trajectories/ndi_na_binding/FES_CM1/"
        pattern = "FES*dat"
        all_fes_files = [file for folder, subdir, files in os.walk(folder)
                         for file in glob(os.path.join(folder, pattern))]
        line = FreeEnergyLine.from_plumed(all_fes_files)
        line.set_datum(datum={'CM1': 7})
        self.assertTrue(0 in line._time_data[0]['energy'])
        self.assertTrue(0 in line._time_data[1]['energy'])
        self.assertTrue(0 in line._time_data[2]['energy'])
        figure = go.Figure()
        trace = go.Scatter(x=line._data[line.cvs[0]], y=line._data['energy'])
        figure.add_trace(trace)
        # figure.show()

    def test_normalise_with_tuple(self):
        """
        testing that the normalise function works with a range
        :return:
        """
        file = "./test_trajectories/ndi_na_binding/FES_CM1.dat"
        line = FreeEnergyLine.from_plumed(file)
        line.set_datum(datum={'CM1': (6, 8)})
        self.assertAlmostEqual(line._data.loc[line._data['CM1'] > 6].loc[line._data['CM1'] < 8]['energy'].mean(), 0)
        figure = go.Figure()
        trace = go.Scatter(x=line._data[line.cvs[0]], y=line._data['energy'])
        figure.add_trace(trace)
        # figure.show()

    def test_normalise_with_tuple_on_time_data(self):
        """
        testing that the normalise function works with a range
        :return:
        """
        folder = "./test_trajectories/ndi_na_binding/FES_CM1/"
        pattern = "FES*dat"
        all_fes_files = [file for folder, subdir, files in os.walk(folder)
                         for file in glob(os.path.join(folder, pattern))]
        line = FreeEnergyLine.from_plumed(all_fes_files)
        line.set_datum(datum={'CM1': (6, 8)})
        self.assertAlmostEqual(line._data.loc[line._data['CM1'] > 6].loc[line._data['CM1'] < 8]['energy'].mean(), 0)
        figure = go.Figure()
        trace = go.Scatter(x=line._data[line.cvs[0]], y=line._data['energy'])
        figure.add_trace(trace)
        # figure.show()

    def test_get_change_over_time(self):
        """
        testing that the normalise function works with a range
        :return:
        """
        folder = "./test_trajectories/ndi_na_binding/FES_CM1/"
        pattern = "FES*dat"
        all_fes_files = [file for folder, subdir, files in os.walk(folder)
                         for file in glob(os.path.join(folder, pattern))]
        line = FreeEnergyLine.from_plumed(all_fes_files)
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
        folder = "./test_trajectories/ndi_na_binding/FES_CM1/"
        pattern = "FES*dat"
        all_fes_files = [file for folder, subdir, files in os.walk(folder)
                         for file in glob(os.path.join(folder, pattern))]
        line = FreeEnergyLine.from_plumed(all_fes_files)
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
        folder = "./test_trajectories/ndi_na_binding/FES_CM1/"
        pattern = "FES*dat"
        all_fes_files = [file for folder, subdir, files in os.walk(folder)
                         for file in glob(os.path.join(folder, pattern))]
        line = FreeEnergyLine.from_plumed(all_fes_files)
        data1 = line.set_datum({'CM1': 3})._data
        data2 = line.set_datum({'CM1': 3})._data
        pd.testing.assert_frame_equal(data2, data1)
        figure = go.Figure()
        trace = go.Scatter(x=line._data[line.cvs[0]], y=line._data['energy'])
        figure.add_trace(trace)
        self.assertTrue(0 in line._data['energy'])
        # figure.show()

    def test_get_error_from_time_dynamics(self):
        """
        testing that the normalise function works with a range
        :return:
        """
        folder = "./test_trajectories/ndi_na_binding/FES_CM1/"
        pattern = "FES*dat"
        all_fes_files = [file for folder, subdir, files in os.walk(folder)
                         for file in glob(os.path.join(folder, pattern))]
        line = FreeEnergyLine.from_plumed(all_fes_files)
        line = line.set_errors_from_time_dynamics(5, bins=100)
        self.assertTrue('energy_err' in line._data.columns.to_list())
        self.assertTrue('population_err' in line._data.columns.to_list())

    def test_bulk_add_lines_alternate_constructor(self):
        """
        testing bulk adding lines to a free energy line
        :return:
        """
        here_dir = "./test_trajectories/ndi_na_binding/"

        shape = FreeEnergySpace.from_standard_directory(here_dir)
        self.assertTrue(type(shape) == FreeEnergySpace)
        self.assertTrue(shape.lines['CM1'].cvs == ['CM1'])
        self.assertTrue(len(shape.trajectories) == 2)

    def test_make_line_from_plumed_bias_exchange(self):
        """
        Testing bias exchange fes read in
        :return:
        """
        file_list = [f for f in glob("./test_trajectories/ndi_bias_exchange/FES_D1/*")]
        fes = FreeEnergyLine.from_plumed(file_list)
        self.assertTrue(type(fes) == FreeEnergyLine)


class TestFreeEnergySurface(unittest.TestCase):

    def test_surface_reader(self):

        file = "./test_trajectories/ndi_na_binding/FES_CM1_D1/FES"
        surface = FreeEnergySurface.from_plumed(file)
        self.assertTrue("energy" in surface._data.columns.to_list())
        self.assertTrue("D1" in surface._data.columns.to_list())
        self.assertTrue("CM1" in surface._data.columns.to_list())
        self.assertTrue("population" in surface._data.columns.to_list())

    def test_surface_datum_with_floats(self):

        file = "./test_trajectories/ndi_na_binding/FES_CM1_D1/FES"
        surface = FreeEnergySurface.from_plumed(file)
        surface = surface.set_datum({'CM1': 0.03, 'D1': 5})
        figure = go.Figure()
        (figure
         .add_trace(go.Contour(
            x=surface._data['CM1'],
            y=surface._data['D1'],
            z=surface._data['energy'],
            colorscale='Jet'))
         )
        self.assertTrue(0 in surface._data['energy'].values.tolist())
        # figure.show()

    def test_surface_reweighting(self):

        surface = FreeEnergySpace.from_standard_directory("./test_trajectories/ndi_na_binding/")
        data = surface.get_reweighted_surface(cvs=["CM2", "CM3"], bins=[-0.5, 0.5, 1.5, 2.5, 3.5]).get_data()
        figure = go.Figure()
        figure.add_trace(go.Heatmap(x=data["CM2"], y=data["CM3"], z=data['energy']))
        figure.update_layout(template='plotly_dark')
        # figure.show()
        self.assertTrue(type(data) == pd.DataFrame)

    def test_surface_reweight_with_symmetry(self):

        space = FreeEnergySpace.from_standard_directory("./test_trajectories/ndi_na_binding/")
        data = (space
                .get_reweighted_surface(cvs=["CM2", "CM3"], bins=[-0.5, 0.5, 1.5, 2.5, 3.5])
                .set_as_symmetric('y=x')
                .get_data()
                )
        figure = go.Figure()
        figure.add_trace(go.Heatmap(x=data["CM2"], y=data["CM3"], z=data['energy']))
        figure.update_layout(template='plotly_dark')
        # figure.show()
        self.assertTrue(type(data) == pd.DataFrame)

    def test_surface_reweight_with_symmetry_err(self):

        space = FreeEnergySpace.from_standard_directory("./test_trajectories/ndi_na_binding/")
        data = (space
                .get_reweighted_surface(cvs=["CM2", "CM3"], bins=[-0.5, 0.5, 1.5, 2.5, 3.5])
                .set_as_symmetric('y=x')
                .get_data()
                )
        figure = go.Figure()
        figure.add_trace(go.Heatmap(x=data["CM2"], y=data["CM3"], z=data['symmetry_error']))
        figure.update_layout(template='plotly_dark')
        # figure.show()
        self.assertTrue(type(data) == pd.DataFrame)

    def test_get_forces(self):

        space = FreeEnergySpace.from_standard_directory("./test_trajectories/ndi_na_binding/")

        force = (space
                 .get_reweighted_surface(cvs=["CM2", "CM3"], bins=[-0.5, 0.5, 1.5, 2.5, 3.5])
                 .set_as_symmetric('y=x')
                 .get_mean_force()
                 .assign(CM2_grad=lambda x: x['CM2_grad']/30)
                 .assign(CM3_grad=lambda x: x['CM3_grad']/30)
                 )

        plt.figure(figsize=(10, 6))
        plt.quiver(force['CM2'], force['CM3'], force['CM2_grad'], force['CM3_grad'], scale=5)

        # plt.show()
        self.assertTrue(type(force) == pd.DataFrame)


class TestFreeEnergySpace(unittest.TestCase):

    landscape = FreeEnergySpace("./test_trajectories/ndi_na_binding/HILLS")

    traj0 = MetaTrajectory("./test_trajectories/ndi_na_binding/COLVAR_REWEIGHT.0")
    traj1 = MetaTrajectory("./test_trajectories/ndi_na_binding/COLVAR_REWEIGHT.1")

    landscape.add_metad_trajectory(traj0)
    landscape.add_metad_trajectory(traj1)

    landscape_opes = FreeEnergySpace("./test_trajectories/ndi_single_opes/Kernels.data")
    traj0_opes = MetaTrajectory("./test_trajectories/ndi_single_opes/COLVAR.0")
    landscape_opes.add_metad_trajectory(traj0_opes)

    def test_make_landscape(self):
        """
        check that landscape constructor works
        :return:
        """
        landscape = FreeEnergySpace("./test_trajectories/ndi_na_binding/HILLS")
        self.assertTrue('height' in landscape._hills.columns.to_list())
        self.assertTrue('time' in landscape._hills.columns.to_list())
        self.assertEqual(type(landscape), FreeEnergySpace)
        self.assertEqual(landscape.cvs, ['D1', 'CM1'])
        self.assertEqual(landscape.n_walker, 2)
        self.assertEqual(landscape.n_timesteps, 50)

    def test_make_landscape_opes(self):
        """
        check that landscape constructor works
        :return:
        """
        landscape = FreeEnergySpace("./test_trajectories/ndi_single_opes/Kernels.data")
        self.assertTrue('height' in landscape._hills.columns.to_list())
        self.assertTrue('time' in landscape._hills.columns.to_list())
        self.assertEqual(type(landscape), FreeEnergySpace)
        self.assertEqual(landscape.cvs, ['D1', 'CM1'])
        self.assertEqual(landscape.n_walker, 1)
        self.assertEqual(landscape.n_timesteps, 40)

    def test_hills_plotter_default_values(self):

        figures = self.landscape.get_hills_figures()
        self.assertTrue(self.landscape._opes is False)
        self.assertEqual(len(figures), 2)
        self.assertTrue(figures[0]._validate)
        self.assertTrue(figures[1]._validate)

    def test_average_hills_figure(self):

        figure = self.landscape.get_average_hills_figure()
        self.assertTrue(figure._validate)

    def test_kernels_plotter_default_values(self):

        figures = self.landscape_opes.get_hills_figures()
        self.assertTrue(self.landscape_opes._opes is True)
        self.assertEqual(len(figures), 1)
        self.assertTrue(figures[0]._validate)

    def test_hills_plotter_default_values_opes(self):

        figures = self.landscape_opes.get_hills_figures()
        self.assertEqual(len(figures), 1)
        self.assertTrue(figures[0]._validate)

    def test_fes_adder_checks_work(self):

        landscape = FreeEnergySpace("./test_trajectories/ndi_na_binding/HILLS")
        fes = FreeEnergyLine.from_plumed("./test_trajectories/ndi_na_binding/FES_CM1.dat")
        landscape.add_line(fes)
        landscape.add_line(fes)
        self.assertEqual(landscape.lines['CM1'], fes)

    def test_traj_adder_checks_work(self):

        landscape = FreeEnergySpace("./test_trajectories/ndi_na_binding/HILLS")
        traj = MetaTrajectory("./test_trajectories/ndi_na_binding/COLVAR.0")
        landscape.add_metad_trajectory(traj)
        landscape.add_metad_trajectory(traj)
        self.assertEqual(landscape.trajectories[0], traj)

    def test_traj_adder_checks_work_opes(self):

        landscape = FreeEnergySpace("./test_trajectories/ndi_single_opes/Kernels.data")
        traj = MetaTrajectory("./test_trajectories/ndi_single_opes/COLVAR.0")
        landscape.add_metad_trajectory(traj)
        landscape.add_metad_trajectory(traj)
        self.assertEqual(landscape.trajectories[0], traj)

    def test_two_bin_reweighted_cv(self):
        """
        Function to test that it is normalising properly when using two bins
        :return:
        """
        fes = self.landscape.get_reweighted_line('D1', bins=[6, 6.5, 7]).set_datum({'D1': 6})
        self.assertEqual(fes._data[fes._data['D1'] == 6.25]['energy'].values[0], 0)

    def test_reweighted_line_cv(self):
        """
        Function to test that it is normalising properly when using two bins
        :return:
        """
        fes = self.landscape.get_reweighted_line('D1', bins=10).set_datum({'D1': 0})
        self.assertEqual(fes._data['energy'].values[0], 0)
        self.assertEqual(round(fes._data['energy'].values[2], 5), 0.88307)
        self.assertEqual(round(fes._data['energy'].values[4], 5), 5.70183)
        self.assertEqual(round(fes._data['energy'].values[6], 5), 11.0178)
        self.assertEqual(round(fes._data['energy'].values[8], 5), 14.10185)

    def test_reweighted_line_cv_adaptive(self):
        """
        Function to test that it is normalising properly when using two bins
        :return:
        """
        fes = self.landscape.get_reweighted_line('D1', bins=10, adaptive_bins=True).set_datum({'D1': 0})
        self.assertEqual(fes._data['energy'].values[0], 0)
        self.assertEqual(round(fes._data['energy'].values[2], 3), 0.435)
        self.assertEqual(round(fes._data['energy'].values[4], 3), 1.863)
        self.assertEqual(round(fes._data['energy'].values[6], 3), 7.596)
        self.assertEqual(round(fes._data['energy'].values[8], 3), 12.964)

    def test_reweighted_line_cv_adaptive_with_walker_err(self):
        """
        Function to test that it is normalising properly when using two bins
        :return:
        """
        fes = self.landscape.get_reweighted_line_with_walker_error('D1', bins=10, adaptive_bins=True)
        self.assertEqual(round(fes._data['energy'].values[0], 3), -6.535)
        self.assertEqual(round(fes._data['energy'].values[2], 3), -6.659)
        self.assertEqual(round(fes._data['energy'].values[4], 3), -4.764)
        self.assertEqual(round(fes._data['energy'].values[6], 3), 0.668)
        self.assertEqual(round(fes._data['energy'].values[8], 3), 6.001)

    def test_two_bin_reweighted_cv_opes(self):
        """
        Function to test that it is normalising properly when using two bins
        :return:
        """
        fes = self.landscape_opes.get_reweighted_line('D1', bins=[0, 0.9, 7]).set_datum({'D1': 0})
        self.assertEqual(fes._data[fes._data['D1'] == 0.45]['energy'].values[0], 0)

    def test_two_bin_reweighted_cv_one_condition(self):
        """
        Function to test that it is normalising properly when using two bins
        :return:
        """
        fes = self.landscape.get_reweighted_line('D1', bins=[6, 6.4, 7], conditions='D1 < 7').set_datum({'D1': 6})
        self.assertEqual(fes._data[fes._data['D1'] == 6.2]['energy'].values[0], 0)

    def test_two_bin_reweighted_cv_opes_one_condition(self):
        """
        Function to test that it is normalising properly when using two bins
        :return:
        """
        fes = self.landscape_opes.get_reweighted_line('D1', bins=[0, 0.9, 7], conditions='D1 < 5').set_datum({'D1': 0})
        self.assertEqual(fes._data[fes._data['D1'] == 0.45]['energy'].values[0], 0)

    def test_two_bin_reweighted_cv_two_condition(self):
        """
        Function to test that it is normalising properly when using two bins
        :return:
        """
        fes = (self
               .landscape
               .get_reweighted_line('D1', bins=[6, 6.4, 7], conditions=['D1 < 7', 'D1 < 6.8'])
               .set_datum({'D1': 6})
               )
        self.assertEqual(fes._data[fes._data['D1'] == 6.2]['energy'].values[0], 0)

    def test_two_bin_reweighted_cv_with_time_stamps(self):
        """
        Function to test that it is normalising properly when using two bins
        :return:
        """
        fes = self.landscape.get_reweighted_line('D1', bins=[6, 6.4, 7], n_timestamps=5).set_datum({'D1': 6})
        self.assertEqual(fes._data[fes._data['D1'] == 6.2]['energy'].values[0], 0)
        self.assertTrue(type(fes._time_data) == dict)
        self.assertTrue(type(fes._time_data[1]) == pd.DataFrame)
        self.assertTrue(type(fes._time_data[3]) == pd.DataFrame)
        self.assertTrue(type(fes._time_data[5]) == pd.DataFrame)

    def test_two_bin_reweighted_cv_with_time_stamps_opes(self):
        """
        Function to test that it is normalising properly when using two bins
        :return:
        """
        fes = self.landscape_opes.get_reweighted_line('D1', bins=[0, 0.9, 7], n_timestamps=5).set_datum({'D1': 0})
        self.assertEqual(fes._data[fes._data['D1'] == 0.45]['energy'].values[0], 0)
        self.assertTrue(type(fes._time_data) == dict)
        self.assertTrue(type(fes._time_data[1]) == pd.DataFrame)
        self.assertTrue(type(fes._time_data[3]) == pd.DataFrame)
        self.assertTrue(type(fes._time_data[5]) == pd.DataFrame)

    def test_temperature_parsed_to_traj(self):
        """
        Function to test that it is normalising properly when using two bins
        :return:
        """
        landscape = FreeEnergySpace("./test_trajectories/ndi_na_binding/HILLS", temperature=320)
        traj0 = MetaTrajectory("./test_trajectories/ndi_na_binding/COLVAR_REWEIGHT.0", temperature=320)
        landscape.add_metad_trajectory(traj0)
        self.assertEqual(landscape.trajectories[0].temperature, 320)

    def test_temperature_parsed_to_traj_opes(self):
        """
        Function to test that it is normalising properly when using two bins
        :return:
        """
        landscape = FreeEnergySpace("./test_trajectories/ndi_single_opes/Kernels.data", temperature=320)
        traj0 = MetaTrajectory("./test_trajectories/ndi_single_opes/COLVAR.0", temperature=320)
        landscape.add_metad_trajectory(traj0)
        self.assertEqual(landscape.trajectories[0].temperature, 320)

    def test_metadata_parsed_to_traj(self):
        """
        Function to test that it is normalising properly when using two bins
        :return:
        """
        landscape = FreeEnergySpace("./test_trajectories/ndi_na_binding/HILLS", metadata=dict(oligomer='NDI'))
        traj0 = MetaTrajectory("./test_trajectories/ndi_na_binding/COLVAR_REWEIGHT.0")
        landscape.add_metad_trajectory(traj0)
        self.assertEqual(landscape.trajectories[0]._metadata['oligomer'], 'NDI')

    def test_metadata_parsed_to_traj_opes(self):
        """
        Function to test that it is normalising properly when using two bins
        :return:
        """
        landscape = FreeEnergySpace("./test_trajectories/ndi_single_opes/Kernels.data", metadata=dict(oligomer='NDI'))
        traj0 = MetaTrajectory("./test_trajectories/ndi_single_opes/COLVAR.0")
        landscape.add_metad_trajectory(traj0)
        self.assertEqual(landscape.trajectories[0]._metadata['oligomer'], 'NDI')

    def test_temperature_parsed_to_line(self):
        """
        Function to test that it is normalising properly when using two bins
        :return:
        """
        landscape = FreeEnergySpace("./test_trajectories/ndi_na_binding/HILLS", temperature=320)
        traj0 = FreeEnergyLine.from_plumed("./test_trajectories/ndi_na_binding/FES_CM1.dat")
        landscape.add_line(traj0)
        self.assertEqual(landscape.lines['CM1'].temperature, 320)

    def test_metad_parsed_to_line(self):
        """
        Function to test that it is normalising properly when using two bins
        :return:
        """
        landscape = FreeEnergySpace("./test_trajectories/ndi_na_binding/HILLS", metadata=dict(oligomer='NDI'))
        traj0 = FreeEnergyLine.from_plumed("./test_trajectories/ndi_na_binding/FES_CM1.dat")
        landscape.add_line(traj0)
        self.assertEqual(landscape.lines['CM1']._metadata['oligomer'], 'NDI')

    def test_bulk_add_trajectories_alternate_constructor_opes(self):
        """
        testing bulk adding trajectories to a free energy line
        :return:
        """
        here_dir = "./test_trajectories/ndi_single_opes/"
        shape = FreeEnergySpace.from_standard_directory(here_dir, colvar_string_matcher="COLVAR.")
        self.assertTrue(type(shape) == FreeEnergySpace)
        self.assertTrue(len(shape.trajectories) == 1)
        self.assertTrue(shape.n_walker == 1)

    def test_bulk_add_trajectories_alternate_constructor(self):
        """
        testing bulk adding trajectories to a free energy line
        :return:
        """
        here_dir = "./test_trajectories/ndi_na_binding/"
        shape = FreeEnergySpace.from_standard_directory(here_dir, colvar_string_matcher="COLVAR.")
        self.assertTrue(type(shape) == FreeEnergySpace)
        self.assertTrue(len(shape.trajectories) == 2)
        self.assertTrue(shape.n_walker == 2)

    def test_bulk_add_trajectories_alternate_constructor_temp_check(self):
        """
        testing bulk adding trajectories to a free energy line
        :return:
        """
        here_dir = "./test_trajectories/ndi_na_binding/"
        shape = FreeEnergySpace.from_standard_directory(here_dir, colvar_string_matcher="COLVAR.", temperature=350)
        line = shape.get_reweighted_line('D1', bins=80)
        self.assertTrue(type(line) == FreeEnergyLine)
        self.assertTrue(type(shape) == FreeEnergySpace)
        self.assertTrue(len(shape.trajectories) == 2)
        self.assertTrue(shape.n_walker == 2)

    def test_bulk_add_trajectories_alternate_constructor_temp(self):
        """
        testing bulk adding trajectories to a free energy line
        :return:
        """
        here_dir = "./test_trajectories/ndi_na_binding/"
        shape = FreeEnergySpace.from_standard_directory(here_dir, colvar_string_matcher="COLVAR.", temperature=350)
        self.assertTrue(type(shape) == FreeEnergySpace)
        self.assertTrue(len(shape.trajectories) == 2)
        self.assertTrue(shape.n_walker == 2)
        self.assertTrue(shape.trajectories[0].temperature == 350)
        self.assertTrue(shape.lines['D1'].temperature == 350)
        self.assertTrue(shape.surfaces[0].temperature == 350)

    def test_bulk_add_trajectories_alternate_constructor_temp_with_metadata(self):
        """
        testing bulk adding trajectories to a free energy line
        :return:
        """
        here_dir = "./test_trajectories/ndi_na_binding/"
        shape = FreeEnergySpace.from_standard_directory(here_dir, colvar_string_matcher="COLVAR.", temperature=350,
                                                        metadata={'ion': 'Na'})
        self.assertTrue(type(shape) == FreeEnergySpace)
        self.assertTrue(len(shape.trajectories) == 2)
        self.assertTrue(shape.n_walker == 2)
        self.assertTrue(shape.trajectories[0].temperature == 350)
        self.assertTrue(shape.lines['D1'].temperature == 350)
        self.assertTrue(shape.surfaces[0].temperature == 350)

    def test_bulk_add_trajectories_alternate_constructor_opes_traj_metad(self):
        """
        testing bulk adding trajectories to a free energy line
        :return:
        """
        here_dir = "./test_trajectories/ndi_single_opes/"
        shape = FreeEnergySpace.from_standard_directory(
            here_dir, colvar_string_matcher="COLVAR.", metadata={'unit': 'NDI'}
        ).trajectories[0].get_data(with_metadata=True)
        self.assertTrue('unit' in shape.columns)

    def test_one_walker_reweighted_with_walker_error(self):
        """
        Function to test that it returns error when only one walker is present.
        :return:
        """
        with self.assertRaises(ValueError):
            self.landscape_opes.get_reweighted_line_with_walker_error("D1", bins=[0, 4, 7]).set_datum({"D1": 0})
    
    def test_two_bin_reweighted_with_walker_error(self):
        """
        Function to test that it weights correctly with errors when using two bins.
        :return:
        """
        fes = self.landscape.get_reweighted_line_with_walker_error("D1", bins=[6, 6.4, 7]).set_datum({"D1": 0})
        self.assertEqual(fes._data[fes._data["D1"] == 6.2]["energy"].values[0], 0)

    def test_two_bin_reweighted_with_walker_one_condition(self):
        """
        Function to test that it reweights correctly with errors when using two bins and one condition.
        :return:
        """
        fes = self.landscape.get_reweighted_line_with_walker_error(
            "D1", bins=[6, 6.4, 7], conditions="D1 < 7").set_datum({"D1": 0}
                                                                 )
        self.assertEqual(fes._data[fes._data["D1"] == 6.2]["energy"].values[0], 0)
    
    def test_two_bin_reweighted_with_walker_two_conditions(self):
        """
        Function to test that it reweights correctly with errors when using two bins and one condition.
        :return:
        """
        fes = self.landscape.get_reweighted_line_with_walker_error(
            "D1", bins=[6, 6.4, 7.0], conditions=["D1 < 7", "D1 < 8"]).set_datum({"D1": 0}
                                                                             )
        self.assertEqual(fes._data[fes._data["D1"] == 6.2]["energy"].values[0], 0)

    def test_reweighted_with_walker_two_conditions(self):
        """
        Function to test that it reweights correctly with errors when using two bins and one condition.
        :return:
        """
        fes = self.landscape.get_reweighted_line_with_walker_error(
            "D1", bins=10, conditions=["D1 < 7", "D1 < 8"]).set_datum({"D1": 0})
        data = fes.get_data().round(4)
        self.assertTrue(type(data) == pd.DataFrame)
        self.assertTrue(data['D1'].iloc[2] == 6.3497)
        self.assertTrue(data['energy'].iloc[6] == 12.1793)
        self.assertTrue(data['energy_err'].iloc[1] == 1.1086)

    def test_bulk_add_trajectories_alternate_constructor_opes_walker_err(self):
        """
        testing bulk adding trajectories to a free energy line
        :return:
        """
        here_dir = "./test_trajectories/ndi_na_binding/"
        shape = FreeEnergySpace.from_standard_directory(here_dir)
        shape.get_reweighted_line_with_walker_error('CM1', bins=200)
        self.assertTrue(shape.n_walker == 2)


class TestFreeEnergySpaceBiasExchange(unittest.TestCase):

    hills = ['./test_trajectories/ndi_bias_exchange/HILLS.0',
             './test_trajectories/ndi_bias_exchange/HILLS.1',
             './test_trajectories/ndi_bias_exchange/HILLS.2',
             './test_trajectories/ndi_bias_exchange/HILLS.3']

    landscape = FreeEnergySpace.from_standard_directory('./test_trajectories/ndi_bias_exchange/', hills_file=hills)

    def test_attributes(self):
        self.assertTrue(self.landscape.cvs == ['D1', 'CM1', 'CM2', 'CM3'])
        self.assertTrue(self.landscape.dt == 0.0004)
        self.assertTrue(self.landscape.max_time == 0.0596)
        self.assertTrue(self.landscape.n_timesteps == 149)
        self.assertTrue(self.landscape.opes is False)
        self.assertTrue(self.landscape.temperature == 298)

    def test_get_hills_figures(self):
        figures = self.landscape.get_hills_figures()
        self.assertTrue(self.landscape._biasexchange is True)
        self.assertEqual(len(figures), 4)

    def test_get_average_hills_figure(self):
        figure = self.landscape.get_average_hills_figure()
        self.assertTrue(type(self.landscape) == FreeEnergySpace)

    def test_get_max_hills_figure(self):
        figure = self.landscape.get_max_hills_figure()
        self.assertTrue(type(self.landscape) == FreeEnergySpace)

    def test_get_hills_figures_hp(self):
        figures = self.landscape.get_hills_figures(height_power=0.5)
        self.assertTrue(self.landscape._biasexchange is True)
        self.assertEqual(len(figures), 4)

    def test_get_CM1_data(self):
        data = self.landscape.lines['CM1'].get_data()
        self.assertTrue(type(data) == pd.DataFrame)

