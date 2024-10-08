import pandas as pd
import numpy as np

# fundamentals
Kb = 0.008314463  # in kJ/mol
R = 8.314         # (J)/(K mol)
F = 96485.332     # (J)/(V mol)

#standard reduction potentials vs SHE
potential_standard_O2_H2O2 = 0.695

# diffusion rates
diff_rate_O2  = 1.9e-5 # cm2/s

# viscosity parameters
viscosity_aq = 0.01 # cm2/s

def boltzmann_energy_to_population(data: pd.DataFrame, x_col: str, temperature: float = 298, y_col: str = 'energy',
                                   y_col_out: str = 'population', discrete_bins: bool = False) -> pd.DataFrame:
    """
    function to take a dataframe and perform a Boltzmann inversion on a column
    :param data: dataframe
    :param temperature: temperature
    :param x_col: column containing the x values of the function
    :param y_col: column containing the y values of the function
    :param y_col_out: the name of the new column with the population distribution
    :param discrete_bins: set to true if x is a non-continuous function
    :return: dataframe with a new column
    """
    data[y_col_out] = np.exp((-data[y_col])/(Kb * temperature))
    if not discrete_bins:
        area = np.trapz(x=data[x_col], y=data[y_col_out])
    else:
        area = np.sum(data[y_col_out])
    data[y_col_out] = data[y_col_out]/np.abs(area)
    return data


def boltzmann_population_to_energy(data: pd.DataFrame, temperature: float = 298, y_col: str = 'population',
                                   y_col_out: str = 'energy') -> pd.DataFrame:
    """
    function to take a dataframe and perform a Boltzmann inversion on a column
    :param data: dataframe
    :param temperature: temperature
    :param y_col: column containing the y values of the function
    :param y_col_out: the name of the new column with the population distribution
    :return: dataframe with a new column
    """
    data[y_col_out] = -np.log(data[y_col]) * Kb * temperature
    return data


def lorentzian(x: list, x0: float, w: float, h: float) -> list[float]:
    """
    lorentzian function
    :param x: a list of x values
    :param x0: the center of the lorentzian
    :param w: width of the lorentzian peak
    :param h: height of the lorentzian peak
    :return: list of y values
    """
    y = [(h/np.pi)*((w/2)/((x-x0)**2 + (w/2)**2)) for x in x]
    return y
