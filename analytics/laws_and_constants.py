import pandas as pd
import numpy as np

Kb = 0.008314463  # in kJ/mol


def boltzmann_invert_energy_to_population(data: pd.DataFrame, x_col: str, temperature: float = 298, y_col: str = 'energy',
                                          y_col_out: str = 'population') -> pd.DataFrame:
    """
    function to take a dataframe and perform a Boltzmann inversion on a column
    :param data: dataframe
    :param temperature: temperature
    :param x_col: column containing the x values of the function
    :param y_col: column containing the y values of the function
    :param y_col_out: the name of the new column with the population distribution
    :return: dataframe with a new column
    """
    data[y_col_out] = np.exp((-data[y_col])/(Kb * temperature))
    data[y_col_out] = data[y_col_out]/np.abs(np.trapz(x=data[x_col], y=data[y_col_out]))

    return data
