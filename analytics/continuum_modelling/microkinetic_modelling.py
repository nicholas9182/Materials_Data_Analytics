import pandas as pd
import numpy as np
from analytics.laws_and_constants import R, F
from analytics.materials.electrolytes import Electrolyte
from analytics.materials.polymers import Polymer
from analytics.materials.solutes import Solute


class MicroKineticModel():
    """
    Top class for a microkinetic model with general functions and attributes that apply to all microkinetic models
    """
    def __init__(self, electrolyte: Electrolyte, polymer: Polymer, rotation_rate: float, temperature: float = 298, pH: float = None) -> None:
        self._electrolyte = electrolyte
        self._polymer = polymer
        self._rotation_rate = rotation_rate
        self._temperature = temperature
        self._f = F / (R * self._temperature)
        self._pH = pH if pH is not None else self._electrolyte.pH
        self._check_tempertures()

    def _check_tempertures(self):
        electrolyte_temp = self._electrolyte.temperature
        model_temp = self._temperature
        if electrolyte_temp != model_temp:
            raise ValueError("Electrolyte and model temperatures are not the same")

    @property
    def temperature(self):
        return self._temperature
    
    @property
    def rotation_rate(self):
        return self._rotation_rate
    
    @property
    def f(self):
        return round(self._f)
    
    @property
    def electrolyte(self):
        return self._electrolyte
    
    @property
    def polymer(self):
        return self._polymer
    
    @property
    def pH(self):
        return self._electrolyte.pH
    
    @property
    def cation(self):
        return self.electrolyte.cation
    
    @property
    def anion(self):
        return self.electrolyte.anion

    @property
    def solvent(self):
        return self.electrolyte.solvent
    
    @property
    def solute(self):
        return self.electrolyte.solute
    
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
            raise ValueError("Forward needs to be boolean!")

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


class OxygenReductionModel(MicroKineticModel):
    
    """
    Class for microkinetic modelling of oxygen reduction reaction in an aqueous electrolyte
    """

    def __init__(self, electrolyte: Electrolyte, polymer: Polymer, rotation_rate: float, temperature: float = 298, pH: float = 7) -> None:
        """
        Initiation function to read in the general properties of an ORR microkinetic model
        """
        super().__init__(electrolyte = electrolyte, polymer = polymer, rotation_rate = rotation_rate, temperature = temperature, pH = pH)

        if self.solvent.name != 'Water':
            raise ValueError("This model is only for aqueous electrolytes")
        if self.solute.name != 'Oxygen':
            raise ValueError("This model is only for oxygen reduction, you need oxygen as a solute")
        
        self._h202 = Solute('H2O2')
        self._ho2 = Solute('HO2')
        self._x = self._calculate_x()

    def _calculate_x(self) -> int:
        """
        Function to determine x in the reaction 2 O2- + xH20 -> HxO2^(x-2) + xOH-.  Depends on the pH and the pka of H2O2 and HO2
        :return x: integer value  
        """

        if self._pH > self._h202.pka and self._pH < 15:
            x = 1
        elif self._pH > self._ho2.pka and self._pH < self._h202.pka:
            x = 2
        elif self._pH < self._h202.pka and self._pH > -1:
            x = 0
        else: 
            raise ValueError("Check your pH value")
        
        return x
    
    @property
    def x(self):
        return self._x
    
    @property
    def h202(self):
        return self._h202
    
    @property
    def ho2(self):
        return self._ho2

