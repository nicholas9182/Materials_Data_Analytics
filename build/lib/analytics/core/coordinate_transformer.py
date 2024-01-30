import numpy as np
import pandas as pd


class CoordinateTransformer:

    def __init__(self, data: pd.DataFrame, x: str = 'x', y: str = 'y', z: str = 'z'):
        """
        class to take coordinate data and perform transformations
        :param data: pandas dataframe with the coordinates
        :param x: column name with the x coordinates
        :param y: column name with the y coordinates
        :param z: column name with the z coordinates
        :return: self
        """

        if x not in data.columns:
            raise ValueError("Your x column isn't in your coordinates")
        elif y not in data.columns:
            raise ValueError("Your y column isn't in your coordinates")
        elif z not in data.columns:
            raise ValueError("Your z column isn't in your coordinates")
        else:
            pass

        if data[x].dtype not in ['int64', 'int32', 'float64', 'float32']:
            raise ValueError("Check your x coordinate dtypes!")
        elif data[y].dtype not in ['int64', 'int32', 'float64', 'float32']:
            raise ValueError("Check your x coordinate dtypes!")
        elif data[z].dtype not in ['int64', 'int32', 'float64', 'float32']:
            raise ValueError("Check your x coordinate dtypes!")
        else:
            pass

        if len(data[x]) == len(data[y]) == len(data[z]):
            pass
        else:
            raise ValueError("Check your column lengths")

        self._x_name = x
        self._y_name = y
        self._z_name = z
        self._n = range(0, len(data[x]))
        self._data = data
        self._coordinates = [np.array([
            [self._data[x].iloc[i]],
            [self._data[y].iloc[i]],
            [self._data[z].iloc[i]]]) for i in self._n]

    @property
    def data(self):
        return self._data

    def rotate(self, theta_x: float = 0, theta_y: float = 0, theta_z: float = 0):
        """
        function to rotate coordinates
        :param theta_x: rotation around the x-axis, in degrees
        :param theta_y: rotation around the y-axis, in degrees
        :param theta_z: rotation around the z-axis, in degrees
        :return:
        """
        # convert degrees to radians
        theta_x = theta_x * 0.0174533
        theta_y = theta_y * 0.0174533
        theta_z = theta_z * 0.0174533

        # define the rotation matrices
        rx = np.array([
            [1, 0, 0],
            [0, np.cos(theta_x), -np.sin(theta_x)],
            [0, np.sin(theta_x), np.cos(theta_x)]
        ])

        ry = np.array([
            [np.cos(theta_y), 0, np.sin(theta_y)],
            [0, 1, 0],
            [-np.sin(theta_y), 0, np.cos(theta_y)]
        ])

        rz = np.array([
            [np.cos(theta_z), -np.sin(theta_z), 0],
            [np.sin(theta_z), np.cos(theta_z), 0],
            [0, 0, 1]
        ])

        r = np.matmul(rz, np.matmul(ry, rx))

        # matrix multiply to get the new coordinates
        new_coords = [np.matmul(r, self._coordinates[i]) for i in self._n]

        # write out to self
        self._data[self._x_name] = [c.flat[0] for c in new_coords]
        self._data[self._y_name] = [c.flat[1] for c in new_coords]
        self._data[self._z_name] = [c.flat[2] for c in new_coords]

        return self
