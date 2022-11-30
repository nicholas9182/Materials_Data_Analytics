import pandas as pd
import numpy as np

Kb = 0.008314463  # in kJ/mol


def boltzmann_invert_energy_to_population(data: pd.DataFrame,  temperature, x_col: str, y_col: str = 'energy', y_col_out: str = 'population'):

    data[y_col_out] = np.exp((-data[y_col])/(Kb * temperature))
    data[y_col_out] = data[y_col_out]/np.abs(np.trapz(x=data[x_col], y=data[y_col_out]))
    return data
