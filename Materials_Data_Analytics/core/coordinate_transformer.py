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
    
    def translate(self, x: float = 0, y: float = 0, z: float = 0):
        """
        function to translate coordinates
        :param x: translation in the x direction
        :param y: translation in the y direction
        :param z: translation in the z direction
        :return:
        """
        self._data[self._x_name] = self._data[self._x_name] + x
        self._data[self._y_name] = self._data[self._y_name] + y
        self._data[self._z_name] = self._data[self._z_name] + z

        return self
    
    def translation_fit(self, reference_data_frame: pd.DataFrame):
        """
        Function to translate the data to the reference data frame to reduce the mean squared error between the two dataframe's coordinates
        """
        if self._data.shape[0] != reference_data_frame.shape[0]:
            raise ValueError("The number of data points in the two dataframes must be equal")
        
        if self._x_name not in reference_data_frame.columns or self._y_name not in reference_data_frame.columns or self._z_name not in reference_data_frame.columns:
            raise ValueError("The reference data frame must have the same column names for the coordinates as the data of the coordinate transformer")
        
        coords = self._data[[self._x_name, self._y_name, self._z_name]].values
        reference = reference_data_frame[[self._x_name, self._y_name, self._z_name]].values
        centroid_coords = np.mean(coords, axis=0)
        centroid_reference = np.mean(reference, axis=0)
        translation_vector = centroid_reference - centroid_coords
        self.translate(translation_vector[0], translation_vector[1], translation_vector[2])
        return self
    
    def rotation_fit(self, reference_data_frame: pd.DataFrame):
        """
        Function to rotate the data to the reference data frame to reduce the mean squared error between the two dataframe's coordinates
        """
        if self._data.shape[0] != reference_data_frame.shape[0]:
            raise ValueError("The number of data points in the two dataframes must be equal")
        
        if self._x_name not in reference_data_frame.columns or self._y_name not in reference_data_frame.columns or self._z_name not in reference_data_frame.columns:
            raise ValueError("The reference data frame must have the same column names for the coordinates as the data of the coordinate transformer")
        
        coords = self._data[[self._x_name, self._y_name, self._z_name]].values
        reference = reference_data_frame[[self._x_name, self._y_name, self._z_name]].values

        centroid_coords = np.mean(coords, axis=0)
        centroid_reference = np.mean(reference, axis=0)

        coords_centered = coords - centroid_coords
        reference_centered = reference - centroid_reference

        h = np.matmul(coords_centered.T, reference_centered)
        u, s, vh = np.linalg.svd(h)
        r = np.matmul(vh.T, u.T)

        theta_y = np.degrees(np.arcsin(-r[2, 0]))
        if np.cos(theta_y) > 1e-6:
            theta_x = np.degrees(np.arctan2(r[2, 1] / np.cos(theta_y), r[2, 2] / np.cos(theta_y)))
            theta_z = np.degrees(np.arctan2(r[1, 0] / np.cos(theta_y), r[0, 0] / np.cos(theta_y)))
        else:
            theta_x = np.degrees(np.arctan2(r[1, 2], r[1, 1]))
            theta_z = 0 

        self.rotate(theta_x=theta_x, theta_y=theta_y, theta_z=theta_z).translate(centroid_reference[0], centroid_reference[1], centroid_reference[2])

        return self
    

class PdbParser:
    """
    Class to handle parsing of pdb files, pring pandas dataframes to pdb files, both single and trajectories
    """
    def __init__(self):
        pass

    @staticmethod
    def _write_pdb(data: pd.DataFrame, filename: str, path: str = '.', 
                   x_col: str = 'x', y_col: str = 'y', z_col: str = 'z', 
                   element_col: str = 'element', time: float = 0.0, model: int = 1, step: int = 1,
                   cell_x: float = 1000.0, cell_y: float = 1000.0, cell_z: float = 1000.0,
                   overwrite: bool = True):
        """
        function to write a pandas dataframe to a pdb file
        :param data: pandas dataframe with the coordinates
        :param filename: name of the pdb file to write
        :param path: path to write the pdb file
        :param x_col: column name with the x coordinates
        :param y_col: column name with the y coordinates
        :param z_col: column name with the z coordinates
        :param element_col: column name with the element type
        :param time: time of the pdb file for when outputting a pdb trajectory
        :param model: model number for when outputting a pdb trajectory
        :param step: step number for when outputting a pdb trajectory
        :param cell_x: x dimension of the simulation box
        :param cell_y: y dimension of the simulation box
        :param cell_z: z dimension of the simulation box
        :param overwrite: boolean to overwrite the file if it exists or to append to a file if it exists
        :return: None
        """

        if filename[-4:] == '.pdb':
            filename = filename[:-4]

        if all(col in data.columns for col in [x_col, y_col, z_col, element_col]) is False:
            raise ValueError("Check your values of x_col, y_col, z_col, and element_col")

        mode = 'w' if overwrite else 'a'

        with open(f"{path}/{filename}.pdb", mode) as f:
            f.write(f"REMARK    GENERATED BY PDBParser in Materials_Data_Analytics \n")
            f.write(f"TITLE     GENERATED BY PDBParser in Materials_Data_Analytics  t=   {time} step= {step} \n")
            f.write(f"REMARK    THIS IS A SIMULATION BOX\n")
            f.write(f"CRYST1  {cell_x:>8.3f}  {cell_y:>8.3f}  {cell_z:>8.3f}  90.00  90.00  90.00 P 1           1\n")
            f.write(f"MODEL        {model}\n")
            for i in range(0, len(data)):
                f.write(f"ATOM  {i+1:5}  {data[element_col].iloc[i]:<2}  UNK     1    {data[x_col].iloc[i]:>8.3f}{data[y_col].iloc[i]:>8.3f}{data[z_col].iloc[i]:>8.3f}  1.00  0.00          {data[element_col].iloc[i]:<2}\n")
            f.write("TER\n")
            f.write("ENDMDL\n")

        return None

    @staticmethod
    def pandas_to_pdb(data: pd.DataFrame, grouping_variables: list = [], **kwargs):
        """ 
        Function to write a pandas dataframe to a pdb file 
        :param data: pandas dataframe with the coordinates
        :grouping_variables: list of columns to group by
        :return: None
        """
        
        if all(col in data.columns for col in grouping_variables) is False:
            raise ValueError("Check your values of grouping_variables")

        if len(grouping_variables)  == 0:
            PdbParser._write_pdb(data, **kwargs)

        basename = kwargs['filename'] + "_"

        if len(grouping_variables) > 0:
            for name, group in data.groupby(grouping_variables):
                for n in name:
                    tag = str(n) + '_'
                tag = tag[:-1]
                kwargs['filename'] = basename + tag
                PdbParser._write_pdb(group, **kwargs)

        return None
    
    @staticmethod
    def pandas_to_pdb_trajectory(data: pd.DataFrame, time_col: str, fit_t0 = False, **kwargs):
        """
        Function to write a pandas dataframe to a pdb trajectory file
        :param data: pandas dataframe with the coordinates
        :param time_col: column name with the time data
        :return: None
        """
        data = data.reset_index().sort_values(by=[time_col, 'index']).drop(columns=['index'])
        t_0 = data[time_col].iloc[0]
        model_0 = 1
        step_0 = 1
        data_0 = data.query(f"{time_col} == {t_0}")

        if fit_t0 is False:
            PdbParser._write_pdb(data_0, time=t_0, model=model_0, step=step_0, overwrite=True, **kwargs)
            for name, group in data.query(f"{time_col} != {t_0}").groupby(time_col):
                model_0 += 1
                step_0 += 1
                PdbParser._write_pdb(group, time=name, model=model_0, step=step_0, overwrite=False, **kwargs)
        else:
            PdbParser._write_pdb(data_0, time=t_0, model=model_0, step=step_0, overwrite=True, **kwargs)
            for name, group in data.query(f"{time_col} != {t_0}").groupby(time_col):
                model_0 += 1
                step_0 += 1
                group = CoordinateTransformer(group).rotation_fit(data_0).data
                PdbParser._write_pdb(group, time=name, model=model_0, step=step_0, overwrite=False, **kwargs)
        
        return None
