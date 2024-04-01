import pandas as pd
import numpy as np


class Path:

    def __init__(self, path: pd.DataFrame):
        """
        initialisation function for a path
        :param path: pandas data frame with the path information
        """
        self._path = path
        self._cvs = path.columns.to_list()
        self._dimensions = len(self._cvs)
        self._time_data = None

    @classmethod
    def from_points(cls, points: list, n_steps: int, cvs: list):
        """
        alternate constructor to make the path from a list of points
        :param points: the points the path should go through
        :param n_steps: the number of steps for the path
        :param cvs: the collective variables the path is defined in
        :return: a Path object
        """
        dimensions = len(cvs)
        for i in points:
            if len(i) != dimensions:
                raise ValueError("All your points need to have the right dimension!")

        if n_steps % len(points) != 0:
            raise ValueError("Make sure your n_steps is a multiple of the number of points")

        segment_steps = int(n_steps / (len(points)-1))

        path_list = []
        for p in range(0, len(points)-1):
            segment_path = pd.DataFrame()
            for v in range(0, len(cvs)):
                segment_path[cvs[v]] = [s for s in np.linspace(points[p][v], points[p+1][v], segment_steps, endpoint=False)]
            path_list.append(segment_path)

        end_points = pd.DataFrame([points[-1]], columns=cvs)
        path_list.append(end_points)
        path = pd.concat(path_list).reset_index(drop=True)
        return cls(path)

    def get_data(self):
        """
        function to get the path data
        :return:
        """
        return self._path.round(3)



