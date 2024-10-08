import pandas as pd
import numpy as np
from Materials_Data_Analytics.metadynamics.free_energy import FreeEnergySurface
from Materials_Data_Analytics.metadynamics.free_energy import FreeEnergyShape


class Path:

    def __init__(self, path: pd.DataFrame, shape: FreeEnergyShape = None):
        """
        initialisation function for a path
        :param path: pandas data frame with the path information
        :param shape: free energy shape
        """
        self._path = path
        self._cvs = path.columns.to_list()
        self._dimensions = len(self._cvs)
        self._time_data = None
        self._shape = shape
        if shape is not None:
            if shape.dimension != self._dimensions:
                raise ValueError("Check your path and shape are the same dimension")

    @classmethod
    def from_points(cls, points: list, n_steps: int, cvs: list, shape: FreeEnergyShape = None):
        """
        alternate constructor to make the path from a list of points
        :param points: the points the path should go through
        :param n_steps: the number of steps for the path
        :param cvs: the collective variables the path is defined in
        :param shape: the shape the path is defined on
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
        return cls(path=path, shape=shape)

    def get_data(self):
        """
        function to get the path data
        :return:
        """
        return self._path.round(3)


class SurfacePath(Path):

    def __init__(self, path: pd.DataFrame, shape: FreeEnergySurface):

        super().__init__(path=path)
        self._shape = shape

        if self._dimensions != 2 or shape.dimension != 2:
            raise ValueError("Check the path and surface both have dimension of 2!")
        if self._cvs[0] not in shape.cvs or self._cvs[1] not in shape.cvs:
            raise ValueError("Check the path and surface cvs are the same!")

    def _get_surface_forces(self, index: int):
        """
        Function to get the force from a free energy surface acting on the i'th point of the path
        :param index: The point to add to get the forces for
        :return: the forces acting on that point
        """
        if index < 0 or index > self._path.index.max():
            raise ValueError("The index needs to be between 0 and the max index")

        surface_forces = self._shape.get_mean_force()
        cv1 = self._cvs[0]
        cv2 = self._cvs[1]

        if index == self._path.index.max() or index == 0:
            fx = 0
            fy = 0
        else:
            x = self._path[cv1].iloc[index]
            y = self._path[cv2].iloc[index]
            fx = FreeEnergyShape.get_nearest_value(surface_forces, {cv1: x}, val_col=f"{cv1}_grad")
            fy = FreeEnergyShape.get_nearest_value(surface_forces, {cv2: y}, val_col=f"{cv2}_grad")

        return np.array([fx, fy])
