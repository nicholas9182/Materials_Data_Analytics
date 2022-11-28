import unittest
import tracemalloc
import os
from glob import glob
import pandas as pd
import plumed as pl
from analytics.metadynamics.free_energy import FreeEnergyLandscape, MetaTrajectory, FreeEnergyLine
tracemalloc.start()


class TestMetaTrajectory(unittest.TestCase):

    def test_colvar_read(self):
        """
        checking that the MetaTrajectory is reading in and processing colvar files correctly. Comparing with a direct plumed read in
        """
        file = "../test_trajectories/ndi_na_binding/COLVAR.0"
        cv_traj = MetaTrajectory(file)
        compare = pd.DataFrame(pl.read_as_pandas("../test_trajectories/ndi_na_binding/COLVAR.0"))

        pd.testing.assert_frame_equal(cv_traj.data, compare)
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

        pd.testing.assert_frame_equal(line.data, compare)
        self.assertEqual(line.cv, 'CM1')

    def test_fes_read_with_time_data(self):
        """
        checking that alternate constructor works for reading in fes data with strides to get time_data dictionary
        """
        folder = "../test_trajectories/ndi_na_binding/FES_CM1/"
        pattern = "FES*dat"
        all_fes_files = [file for folder, subdir, files in os.walk(folder) for file in glob(os.path.join(folder, pattern))]
        line = FreeEnergyLine.with_strides(all_fes_files)
        compare = pd.DataFrame(pl.read_as_pandas("../test_trajectories/ndi_na_binding/FES_CM1/FES23.dat"))
        pd.testing.assert_frame_equal(line.time_data[23], compare)

    def test_fes_plot_line_with_normalise(self):
        """
        checking that the plotter normalises properly
        :return:
        """
        line = FreeEnergyLine("../test_trajectories/ndi_na_binding/FES_CM1.dat")
        plot = line.plot_line(ymax=-40, normalise=0)
        self.assertTrue(plot._validate)


class TestFreeEnergyLandscape(unittest.TestCase):

    def test_make_landscape(self):
        """
        check that landscape constructor works
        :return:
        """
        landscape = FreeEnergyLandscape("../test_trajectories/ndi_na_binding/HILLS")
        self.assertTrue('height' in landscape.hills.columns.to_list())
        self.assertTrue('time' in landscape.hills.columns.to_list())
        self.assertEqual(type(landscape), FreeEnergyLandscape)
        self.assertEqual(landscape.cvs, ['D1', 'CM1'])
        self.assertEqual(landscape.n_walker, 8)
        self.assertEqual(landscape.n_timesteps, 2978)

    def test_hills_plotter_default_values(self):

        landscape = FreeEnergyLandscape("../test_trajectories/ndi_na_binding/HILLS")
        figures = landscape.get_hills_figures()
        self.assertEqual(len(figures), 8)
        self.assertTrue(figures[0]._validate)
        self.assertTrue(figures[1]._validate)
        self.assertTrue(figures[2]._validate)
        self.assertTrue(figures[3]._validate)

    def test_fes_adder_checks_work(self):

        landscape = FreeEnergyLandscape("../test_trajectories/ndi_na_binding/HILLS")
        fes = FreeEnergyLine("../test_trajectories/ndi_na_binding/FES_CM1.dat")
        landscape.add_fes_line(fes)
        landscape.add_fes_line(fes)
        self.assertEqual(landscape.lines[0], fes)

    def test_traj_adder_checks_work(self):

        landscape = FreeEnergyLandscape("../test_trajectories/ndi_na_binding/HILLS")
        traj = MetaTrajectory("../test_trajectories/ndi_na_binding/COLVAR.0")
        landscape.add_metad_trajectory(traj)
        landscape.add_metad_trajectory(traj)
        self.assertEqual(landscape.trajectories[0], traj)
