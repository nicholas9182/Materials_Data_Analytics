import unittest
import numpy as np
import pandas as pd
from analytics.metadynamics.path_analysis import Path
from analytics.metadynamics.free_energy import FreeEnergySpace
from analytics.metadynamics.path_analysis import SurfacePath


class TestPath(unittest.TestCase):

    mep = pd.DataFrame({
        "CM6": [f for f in np.linspace(0, 8, 10, endpoint=False)],
        "CM7": [f for f in np.linspace(8, 0, 10, endpoint=True)]
    })

    def test_path(self):

        path = Path(self.mep)
        self.assertTrue(type(path) == Path)

    def test_from_points(self):

        points = [[8, 0], [0, 0], [0, 8]]
        cvs = ['CM6', 'CM7']
        path = Path.from_points(points, n_steps=21, cvs=cvs).get_data()
        self.assertTrue(type(path) == pd.DataFrame)
        self.assertTrue(path['CM6'].iloc[4] == 4.8)
        self.assertTrue(path['CM7'].iloc[4] == 0.0)

    def test_get_surface_force(self):

        points = [[3, 0], [0, 0], [0, 3]]
        cvs = ['CM2', 'CM3']

        surface = (FreeEnergySpace
                   .from_standard_directory("./test_trajectories/ndi_na_binding/")
                   .get_reweighted_surface(cvs=["CM2", "CM3"], bins=[-0.5, 0.5, 1.5, 2.5, 3.5])
                   .set_as_symmetric('y=x')
                   )

        path = SurfacePath.from_points(points, n_steps=21, cvs=cvs, shape=surface)

        forces = path._get_surface_forces(index=5)

        self.assertTrue(type(forces) == np.ndarray)
