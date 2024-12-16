# Continuum Modelling

The continuum_modelling module provides frameworks for modelling materials properties and comparing modells to experiments.

## Microkinetic Modelling

Microkinetic modelling is used to model ring and disk currents during a rotating ring-disk electrode experiment. Generally the model provides the following information; the expected currents at the ring and disk, the proportion of each species present at the electrode surface, and the rate of each reaction.

### Authors

 - Dr. Nicholas Siemons ([GitHub](https://github.com/nicholas9182)) (nsiemons@stanford.edu)

### Getting Started with ECpD

Currently the only supported electrochemical model is the ECpD model, which is defined in publication "Origins of hydrogen peroxide selectivity
during oxygen reduction on organic mixed ionic–electronic conducting polymers", DOI: 10.1039/d3ee02102e.

Initially an electrolyte must be created using the Electrolyte class,

```python
from Materials_Data_Analytics.materials.solutes import MolecularOxygen
from Materials_Data_Analytics.materials.solvents import Solvent 
from Materials_Data_Analytics.materials.ions import Cation, Anion 
from Materials_Data_Analytics.materials.electrolytes import Electrolyte
from Materials_Data_Analytics.continuum_modelling.microkinetic_modelling import ECpD 

na_cation = Cation(name='Na+')
cl_anion = Anion(name='Cl-')
water_solvent = Solvent('water')
oxygen_solute = MolecularOxygen()

my_electrolyte = Electrolyte(solvent=water_solvent, 
                            cation=na_cation, 
                            anion=cl_anion, 
                            concentrations={na_cation: 0.1, cl_anion: 0.1, oxygen_solute: 0.0008}, 
                            solute=oxygen_solute, 
                            pH=14.2, 
                            temperature=298,
                            diffusivities={oxygen_solute: 0.000019},
                            viscosity=0.01
                            )
```

followed by the creation of the polymer object,

```python
from Materials_Data_Analytics.materials.polymers import NType

my_polymer = NType('NDI', formal_reduction_potential = -0.3159)
```

followed by creation of the model

```python
my_model = ECpD(electrolyte=my_electrolyte, polymer=my_polymer, rotation_rate=1600)
```

where a rotation rate for the electrode is required.

### Attributes

The ECpD class contains various attributes which can easily be accessed:
 - ```viscosity``` - The viscosity of the electrolyte
 - ```diffusivities``` - The diffusivities of the solutes in the electrolyte
 - ```temperature``` - The temperature of the electrolyte
 - ```pH``` - The pH of the electrolyte
 - ```rotation_rate``` - The rotation rate of the electrode
 - ```electrolyte``` - The electrolyte object
 - ```polymer``` - The polymer object
 - ```f``` - F/RT
 - ```cation``` - The cation in the electrolyte
 - ```anion``` - The anion in the electrolyte
 - ```solute``` - The solute in the electrolyte
 - ```solvent``` - The solvent in the electrolyte
 - ```mass_transfer_coefficient``` - The mass transfer coefficient of the electrolyte
 - ```diffusion_layer_thickness``` - The thickness of the diffusion layer

### Methods

The ECpD model can be used to calculate the expected currents at the ring and disk, the proportion of each species present at the electrode surface, and the rate of each reaction.

```python
my_model.calculate_k2() # calculate the rate constant for the second reaction
my_model.calculate_k3() # calculate the rate constant for the third reaction
my_model.calculate_ksf1(E) # calculate the forward electrochemical rate constant for the first reaction
my_model.calculate_ksb1(E) # calculate the backward electrochemical rate constant for the first reaction
my_model.calculate_v1(E, k01, beta, thetaN, thetaP) # calculate the rate of the first reaction
my_model.calculate_v2(kf2, thetaP, thetaN, CS02, CS02_superoxide) # calculate the rate of the second reaction
my_model.calculate_v3(kf3, CS02_superoxide) # calculate the rate of the third reaction
my_model.get_disk_current_density(v1) # get the disk current density from the rate constant
my_model.get_ring_current_density(v3, Neff) # get the ring current density from the rate constant
my_model.solve_parameters(E, k01, kf2, kf3, beta) # solve the reaction rate parameters for the model
my_model.get_e_sweep(E_min, E_max, E_n, k01, kf2, kf3, beta) # get the current densities and coverages for a potential range
```


