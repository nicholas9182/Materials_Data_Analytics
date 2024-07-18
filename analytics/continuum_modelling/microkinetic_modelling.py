import pandas as pd
import numpy as np
from analytics.laws_and_constants import R, F, pka_H2O2, pka_HO2, potential_formal_O2_supO2, potential_standard_O2_H2O2, diff_rate_O2, viscosity_aq


class MicroKineticModel():

    """
    Microkinetic model containing high level functions which apply to all microkinetic modelling
    """

    def __init__(self, pH: float, rotation_rate: float, bulk_concentration: float, temperature: float = 298) -> None:
        self._pH = pH
        self._rotation_rate = rotation_rate
        self._bulk_concentration = bulk_concentration
        self._temperature = temperature
        self._f = F / (R * self._temperature)

    @property
    def temperature(self):
        return self._temperature
    
    @property
    def pH(self):
        return self._pH
    
    @property
    def rotation_rate(self):
        return self._rotation_rate
    
    @property
    def bulk_concentration(self):
        return self._bulk_concentration
    
    @property
    def f(self): 
        return round(self._f)
    
    def _calculate_rate_constant(self, E1: float, E2: float, n1: int = 1, n2: int = 1) -> float:
        """
        Function to calculate a rate constant from two formal reduction potentials of species
        :param E1: the formal reduction potential of species 1
        :param E2: the formal reduction potential of species 2
        :param n1: stoichiometric number of species 1
        :param n2: stoichiometric number of species 2 
        :return rate_constant: the rate_constant of the reaction
        """
        rate_constant = np.exp(self._f*(n2*E2 - n1*E1))
        return rate_constant
    
    def _calculate_electrochemical_rate_constant(self, E1: float, E: float, beta: float, rate_constant_zero_overvoltage: float, forward: bool):
        """
        Function to calculate an electrochemical rate constant at a given potential
        :param E1: The formal reduction voltage of the reduced species based on pH and reaction conditions
        :param E: The applied potential
        :param beta: the symmetry coefficient
        :param rate_constant_zero_overvoltage: rate constant at zero overvoltage for reaction
        """
        if forward is True:
            alpha = -beta
        elif forward is False:
            alpha = 1-beta
        else: 
            return ValueError("Forward needs to be boolean!")

        return rate_constant_zero_overvoltage * np.exp(alpha * self._f * (E - E1))
    
    def _calculate_diffusion_layer_thickness(self, diff_rate: float, viscosity: float) -> float:
        """
        Function to calculate the depth of the diffusion layer according to Levich model
        :param diff_rate: The diffusion rate of the reacting species in cm2/s
        :param viscosity: the viscosity of the electrolyte in cm2/s
        """
        return 1.61 * diff_rate**(1/3) * (self._rotation_rate*2*np.pi/60)**(-1/2) * viscosity**(1/6)
    
    def _calculate_mass_transfer_coefficient(self, diff_rate: float, diff_layer_thickness: float) -> float: 
        """
        Function to calculate the mass transfer coefficient
        :param diff_rate: The diffusion rate of the reacting species in cm2/s
        :param diff_layer_thickness: the diffusion layer thickness
        """
        return diff_rate/diff_layer_thickness
    