import pandas as pd
import numpy as np
import torch
from Materials_Data_Analytics.laws_and_constants import R, F
from Materials_Data_Analytics.materials.electrolytes import Electrolyte
from Materials_Data_Analytics.materials.polymers import Polymer, NType
from Materials_Data_Analytics.materials.solutes import Solute
from scipy.optimize import fsolve


class MicroKineticModel():
    """
    Top class for a microkinetic model with general functions and attributes that apply to all microkinetic models
    """
    def __init__(self, 
                 electrolyte: Electrolyte, 
                 polymer: Polymer | list[Polymer], 
                 rotation_rate: float
                 ) -> None:

        self._electrolyte = electrolyte
        self._polymer = polymer if type(polymer) == list else [polymer]
        self._rotation_rate = rotation_rate
        self._pH = self.electrolyte._pH
        self._temperature = self.electrolyte._temperature
        self._solute = self._electrolyte._solute
        self._f = F / (R * self._temperature)
        self._diffusivities = self._electrolyte._diffusivities
        self._viscosity = self._electrolyte._viscosity

    @property
    def viscosity(self):
        return self._viscosity

    @property  
    def diffusivities(self):
        return self._diffusivities

    @property
    def temperature(self):
        return self._electrolyte.temperature
    
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
        return self._polymer if len(self._polymer) > 1 else self._polymer[0]
    
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
    
    def _calculate_electrochemical_rate_coefficient(self, E1: float, E: float, beta: float, rate_constant_zero_overvoltage: float, forward: bool):
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
    
    def _calculate_diffusion_layer_thickness(self, diff_rate: float) -> float:
        """
        Function to calculate the depth of the diffusion layer according to Levich model
        :param diff_rate: The diffusion rate of the reacting species in cm2/s
        """
        return 1.61 * diff_rate**(1/3) * (self._rotation_rate*2*np.pi/60)**(-1/2) * self.electrolyte._viscosity**(1/6)
    
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

    def __init__(self, 
                 electrolyte: Electrolyte, 
                 polymer: Polymer, 
                 rotation_rate: float
                 ) -> None:

        super().__init__(electrolyte = electrolyte, polymer = polymer, rotation_rate = rotation_rate)

        if self.solvent.name != 'Water':
            raise ValueError("This model is only for aqueous electrolytes")
        if self.solute.name != 'Oxygen':
            raise ValueError("This model is only for oxygen reduction, you need oxygen as a solute")
        
        if self._electrolyte._diffusivities is None:
            raise ValueError("You need to provide the diffusivities of the solutes in the electrolyte for an ORR model")
        
        self._h202 = Solute('H2O2')
        self._ho2 = Solute('HO2')
        self._x = self._calculate_x()
        self._o2 = [i for i in self._solute if i.name == 'Oxygen'][0]
        self._o2._diffusion_layer_thickness = self._calculate_diffusion_layer_thickness(self._electrolyte._diffusivities[self._o2])
        self._o2._mass_transfer_coefficient = self._electrolyte._diffusivities[self._o2]/self._o2._diffusion_layer_thickness

    @property
    def mass_transfer_coefficient(self):
        return self._o2._mass_transfer_coefficient

    @property
    def diffusion_layer_thickness(self):
        return self._o2._diffusion_layer_thickness

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
    
    @property
    def o2(self):
        return self._o2


class ECpD(OxygenReductionModel):

    """
    Class to model an ORR reaction with a polymer in the electrolyte following an ECpD mechanism
    """

    def __init__(self, 
                 electrolyte: Electrolyte, 
                 polymer: NType, 
                 rotation_rate: float
                 ) -> None:
        
        super().__init__(electrolyte = electrolyte, polymer = polymer, rotation_rate = rotation_rate)
        

    def calculate_k2(self) -> float:
        """
        Function to calculate the equilibrium constant for the Cp step of the reaction pathway. 
        Assumes the polymer with the most positive formal reduction potential is the one that is reduced
        """
        reduction_potentials = [i._formal_reduction_potential for i in self._polymer]
        E1 = max(reduction_potentials)
        E2 = self.o2._formal_reduction_potentials['O2_superoxide']

        return self._calculate_rate_constant(E1 = E1, E2 = E2, n1 = 1, n2 = 1)
    
    def calculate_k3(self) -> float:
        """
        Function to calculate the equilibrium constant for the D step of the reaction pathway
        """
        E1 = self.o2._formal_reduction_potentials['O2_superoxide']
        E2 = self.o2.calculate_formal_potential_O2_HXO2(T = self.temperature, x = self._x, pH = self.pH)

        return self._calculate_rate_constant(E1 = E1, E2 = E2, n1 = 2, n2 = 1)

    def calculate_ksf1(self, E: float, k01: float, beta: float) -> float:
        """
        Function to calculate the electrochemical rate constant for the forward reaction of the first step
        Assumes the polymer with the most positive formal reduction potential is the one that is reduced
        :param E: the applied potential
        :param k01: the rate constant at zero overvoltage
        :param beta: the symmetry coefficient
        :return ksf1: the electrochemical rate constant for the forward reaction of the first step
        """
        reduction_potentials = [i._formal_reduction_potential for i in self._polymer]
        E1 = max(reduction_potentials)
        return self._calculate_electrochemical_rate_coefficient(E1 = E1, E = E, beta = beta, rate_constant_zero_overvoltage = k01, forward = True)

    def calculate_ksb1(self, E: float ,k01: float, beta: float) -> float:
        """
        Function to calculate the electrochemical rate constant for the backward reaction of the first step
        :param E: the applied potential
        :param k01: the rate constant at zero overvoltage
        :param beta: the symmetry coefficient
        :return ksb1: the electrochemical rate constant for the backward reaction of the first step
        """
        reduction_potentials = [i._formal_reduction_potential for i in self._polymer]
        E1 = max(reduction_potentials)
        return self._calculate_electrochemical_rate_coefficient(E1 = E1, E = E, beta = beta, rate_constant_zero_overvoltage = k01, forward = False)
    
    def calculate_v1(self, E: float, k01: float, beta: float, thetaN: float, thetaP: float) -> float:
        """
        Function to calculate the rate expression of the first step of the reaction
        :param E: the applied potential
        :param rate_const_zero_overvoltage: the rate constant at zero overvoltage
        :param beta: the symmetry coefficient
        :param thetaN: the coverage of the polymer in its neutral state
        :param thetaP: the coverage of the polymer in its polaron state
        :return v1: the rate expression of the first step
        """
        first_term = self.calculate_ksf1(E = E, k01 = k01, beta = beta) * thetaN
        second_term = self.calculate_ksb1(E = E, k01 = k01, beta = beta) * thetaP
        return first_term - second_term
    
    def calculate_v2(self, kf2: float, thetaP: float, thetaN: float, CS02: float, CS02_superoxide: float) -> float:
        """
        Function to calculate the rate expression of the second step of the reaction
        :param kf2: the rate constant of the second step
        :param thetaP: the coverage of the polymer in its polaron state
        :param thetaN: the coverage of the polymer in its neutral state
        :param CS02: the concentration of O2 at the surface
        :param CS02_superoxide: the concentration of O2 in the superoxide state at the surface
        :return v2: the rate expression of the second step
        """
        kb2 = kf2/self.calculate_k2()
        first_term = kf2 * thetaP * CS02 
        second_term = kb2 * thetaN * CS02_superoxide
        return first_term - second_term
    
    def calculate_v3(self, kf3: float, CS02_superoxide: float):
        """
        Function to calculate the rate expression of the third step of the reaction
        :param kf3: the rate constant of the third step
        :param CS02_superoxide: the concentration of O2 in the superoxide state at the surface
        :return v3: the rate expression of the third step
        """
        return kf3 * CS02_superoxide**2
    
    def get_disk_current_density(self, v1: float):
        """
        Function to get the current density at the disk electrode
        :param v1: the rate expression of the first step
        """
        return -F * v1
    
    def get_ring_current_density(self, v3: float, Neff = 0.25):
        """
        Function to get the current density at the ring electrode
        :param v3: the rate expression of the third step
        :param Neff: the collection efficiency of the ring electrode
        """
        return 2 * F * v3 * Neff
    
    def _calculate_rate_equations(self, vars, E: float, k01: float, beta: float, kf2: float, kf3: float, thetaN: float, thetaP: float, CS02: float, CS02_superoxide: float) -> tuple:
        """
        Function to calculate the rate equations for the ECpD model
        :param vars: the variables to solve for
        :param E: the applied potential
        :param k01: the rate constant at zero overvoltage
        :param beta: the symmetry coefficient
        :param kf2: the rate constant of the second step
        :param kf3: the rate constant of the third step
        :param thetaN: the coverage of the polymer in its neutral state
        :param thetaP: the coverage of the polymer in its polaron state
        :param CS02: the concentration of O2 at the surface
        :param CS02_superoxide: the concentration of O2 in the superoxide state at the surface
        """
        v1 = self.calculate_v1(E = E, k01 = k01, beta = beta, thetaN = thetaN, thetaP = thetaP)
        v2 = self.calculate_v2(kf2 = kf2, thetaP = thetaP, thetaN = thetaN, CS02 = CS02, CS02_superoxide = CS02_superoxide)
        v3 = self.calculate_v3(kf3 = kf3, CS02_superoxide = CS02_superoxide)
        return v1, v2, v3
    
    def solve_parameters(self, E: float, k01: float, kf2: float, kf3: float, beta: float):
        """
        Function to solve for the parameters of the ECpD model
        :param k01: the rate constant of the first step
        :param kf2: the rate constant of the second step
        :param kf3: the rate constant of the third step
        :param thetaN: the coverage of the polymer in its neutral state
        :param thetaP: the coverage of the polymer in its polaron state
        :param CS02: the concentration of O2 at the surface
        :param CS02_superoxide: the concentration of O2 in the superoxide state at the surface
        """
        def equations(vars):

            thetaN, thetaP, CS02, CS02_superoxide = vars
            
            v1, v2, v3 = self._calculate_rate_equations(vars = vars, 
                                                        E = E, 
                                                        k01 = k01, 
                                                        beta = beta, 
                                                        kf2 = kf2, 
                                                        kf3 = kf3, 
                                                        thetaN = thetaN, 
                                                        thetaP = thetaP, 
                                                        CS02 = CS02, 
                                                        CS02_superoxide = CS02_superoxide)

            eq16 = v2 - v3 - self.mass_transfer_coefficient * (self.electrolyte._concentrations[self._o2] - CS02_superoxide)
            eq17 = 2*v3 - v2 + self.mass_transfer_coefficient * CS02_superoxide**2
            eq18 = v1 - v2
            eq19 = thetaN + thetaP - 1

            return eq16, eq17, eq18, eq19
        
        guess = [0.5, 0.5, 0, 0]
        solution = fsolve(equations, guess)

        return {'thetaN': solution[0], 'thetaP': solution[1], 'CS02': solution[2], 'CS02_superoxide': solution[3]}
    
    def solve_parameters2(self, E: float, k01: float, kf2: float, kf3: float, beta: float):
        """
        Function to solve for the parameters of the ECpD model
        :param k01: the rate constant of the first step
        :param kf2: the rate constant of the second step
        :param kf3: the rate constant of the third step
        :param thetaN: the coverage of the polymer in its neutral state
        :param thetaP: the coverage of the polymer in its polaron state
        :param CS02: the concentration of O2 at the surface
        :param CS02_superoxide: the concentration of O2 in the superoxide state at the surface
        """
        def equations(vars):

            thetaN, thetaP, CS02, CS02_superoxide = vars
            
            v1, v2, v3 = self._calculate_rate_equations(vars = vars, 
                                                        E = E, 
                                                        k01 = k01, 
                                                        beta = beta, 
                                                        kf2 = kf2, 
                                                        kf3 = kf3, 
                                                        thetaN = thetaN, 
                                                        thetaP = thetaP, 
                                                        CS02 = CS02, 
                                                        CS02_superoxide = CS02_superoxide)

            eq16 = v2 - v3 - self.mass_transfer_coefficient * (self.electrolyte._concentrations[self._o2] - CS02_superoxide)
            eq17 = 2*v3 - v2 + self.mass_transfer_coefficient * CS02_superoxide**2
            eq18 = v1 - v2
            eq19 = thetaN + thetaP - 1

            return eq16, eq17, eq18, eq19
        
        guess = [0.5, 0.5, 0, 0]
        solution = fsolve(equations, guess)

        return {'thetaN': solution[0], 'thetaP': solution[1], 'CS02': solution[2], 'CS02_superoxide': solution[3]}
    
    def get_e_sweep(self, E_min: float, E_max: float, E_n: 20, k01: float, kf2: float, kf3: float, beta: float):
        """
        Function to get the parameters of the ECpD model for a range of potentials
        :param E_min: the minimum potential
        :param E_max: the maximum potential
        :param E_n: the number of potentials to calculate
        :param k01: the rate constant of the first step
        :param kf2: the rate constant of the second step
        :param kf3: the rate constant of the third step
        :param beta: the symmetry coefficient
        """
        data = (pd
                .DataFrame({'potential': np.linspace(E_min, E_max, E_n)})
                .assign(
                    params = lambda x: [self.solve_parameters2(E = E, k01 = k01, kf2 = kf2, kf3 = kf3, beta = beta) for E in x['potential']],
                    thetaN = lambda x: [i['thetaN'] for i in x['params']],
                    thetaP = lambda x: [i['thetaP'] for i in x['params']],
                    CS02 = lambda x: [i['CS02'] for i in x['params']],
                    CS02_superoxide = lambda x: [i['CS02_superoxide'] for i in x['params']]
                    )
                .assign(
                    v1 = lambda x: [self.calculate_v1(E = x['potential'].iloc[i], k01 = k01, beta = beta, 
                                                      thetaN = x['thetaN'].iloc[i], thetaP = x['thetaP'].iloc[i]) for i in range(len(x))],
                    v2 = lambda x: [self.calculate_v2(kf2 = kf2, thetaP = x['thetaP'].iloc[i], thetaN = x['thetaN'].iloc[i],
                                                        CS02 = x['CS02'].iloc[i], CS02_superoxide = x['CS02_superoxide'].iloc[i]) for i in range(len(x))],
                    v3 = lambda x: [self.calculate_v3(kf3 = kf3, CS02_superoxide = x['CS02_superoxide'].iloc[i]) for i in range(len(x))],
                    disk_current_density = lambda x: [self.get_disk_current_density(v1 = x['v1'].iloc[i]) for i in range(len(x))],
                    ring_current_density = lambda x: [self.get_ring_current_density(v3 = x['v3'].iloc[i]) for i in range(len(x))]
                    )
                )
                
        return data

