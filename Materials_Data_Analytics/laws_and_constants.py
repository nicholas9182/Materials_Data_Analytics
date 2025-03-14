import pandas as pd
import numpy as np

# fundamentals
KB = 1.380649e-23         # Boltzmann constant in J/K
R = 8.31446261815324      # Universal gas constant in J/(mol K)
F = 96485.3321233100184   # Faraday constant in C/mol
E = 1.602176634e-19       # Elementary charge in C
NA = 6.02214076e23        # Avogadro's number in mol^-1

#standard reduction potentials vs SHE
POTENTIAL_STANDARD_O2_H2O2 = 0.695

# diffusion rates
DIFF_RATE_O2  = 1.9e-5 # cm2/s

# viscosity parameters
VISCOSITY_AQ = 0.01 # cm2/s

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
    data[y_col_out] = np.exp((-data[y_col])/(KB*NA/1000 * temperature))
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
    data[y_col_out] = -np.log(data[y_col]) * KB * NA / 1000 * temperature
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
