import unittest
import tracemalloc
import os
from glob import glob
import pandas as pd
import plumed as pl
from analytics.metadynamics.free_energy import FreeEnergyLandscape, MetaTrajectory, FreeEnergyLine
os.environ["TF_CPP_MIN_LOG_LEVEL"] = "2"
tracemalloc.start()


class TestMetaTrajectory(unittest.TestCase):

    def test_colvar_read(self):
        """
        checking that the MetaTrajectory is reading in and processing colvar files correctly. Comparing with a direct plumed read in
        """
        file = "../test_trajectories/ndi_na_binding/COLVAR.0"
        walker = int(file.split(".")[-1])
        cv_traj = MetaTrajectory(file, walker_num=walker)
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


class TestFreeEnergyLandscape(unittest.TestCase):

    def test_make_landscape(self):

        landscape = FreeEnergyLandscape("../test_trajectories/ndi_na_binding/HILLS")
        self.assertTrue('height' in landscape.hills.columns.to_list())
        self.assertTrue('time' in landscape.hills.columns.to_list())
        self.assertEqual(type(landscape), FreeEnergyLandscape)
        self.assertEqual(landscape.cvs, ['D1', 'CM1'])
        self.assertEqual(landscape.n_walker, 8)
        self.assertEqual(landscape.n_timesteps, 2978)
