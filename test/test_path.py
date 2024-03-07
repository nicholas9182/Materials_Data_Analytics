import unittest
import numpy as np
import pandas as pd
from analytics.metadynamics.path_analysis import Path


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

