# Quantum Chemical Simulations

The quantum_chemistry module contains classes for parsing and analysing quantum chemical simulations.

## Gaussian

### Authors

 - Dr. Nicholas Siemons ([GitHub](https://github.com/nicholas9182)) (nsiemons@stanford.edu)

### Getting Started

A GaussianParser object can be created either from a gaussian log output file, or from a list of output files (if, for example calculation restarts are rquired). At the moment it is assumed that the flag #p was included in the keywords line of the gaussian input file, indicating the log file to be more verbose in its output.  

All energies are returned in eV, and distances in Angstroms.

```python
from Materials_Data_Analytics.quantum_chemistry.gaussian import GaussianParser

my_gaussian = GaussianParser('path/to/logfile.log')
my_gaussian = GaussianParser(['path/to/logfile1.log', 'path/to/logfile2.log'])
```

### Attributes
The GausianParser class contains various attributes which can easily be accessed:
 - ```complete``` - A boolean indicating whether the calculation completed successfully
 - ```stable``` - A string indicating the stability of the system
 - ```restart``` - A boolean indicating whether the calculation was restarted
 - ```keywords``` - A list of the keywords used in the calculation
 - ```basis_set``` - The basis set used in the calculation
 - ```functional``` - The functional used in the calculation
 - ```charge``` - The charge of the system
 - ```multiplicity``` - The multiplicity of the system
 - ```pop``` - A boolean indicating whether population analysis was performed
 - ```solvent``` - A boolean indicating whether a solvent was used
 - ```esp``` - A boolean indicating whether an ESP calculation was performed
 - ```freq``` - A boolean indicating whether a frequency calculation was performed
 - ```raman``` - A boolean indicating whether a raman calculation was performed
 - ```scf_iterations``` - The number of scf iterations
 - ```energy``` - The final DFT energy of the system (in eV)
 - ```atom_count``` - The number of atoms in the system
 - ```heavyatomcount``` - The number of heavy atoms in the system
 - ```atoms``` - A list with the elements in the system
 - ```heavy_atoms``` - A list with the heavy atoms in the system
 - ```time_stamp``` - The time stamp of the calculation
 - ```n_alpha``` - The number of alpha electrons
 - ```n_beta``` - The number of beta electrons
 - ```n_electrons``` - The number of electrons
 - ```homo``` - The HOMO energy level with reference to the vacuum level
 - ```lumo``` - The LUMO energy level with reference to the vacuum level
 - ```bandgap``` - The bandgap of the material

 If a vibrational analysis has been done, the parser will contain the additional keywords:
 - ```thermal_energy_corrections``` - A dictionary of energy corrections calculated during the thermochemical analysis
 - ```free_energy``` - the free energy of the molecule as calculated via a thermochemical analysis


### Methods
Once created, various operations can be performed on the object. These checking the SCF energy during an optimization or checking for spin contamination:

```python
scf_energies = my_gaussian.get_scf_convergence() # extract the energies during the scf 
scf_coordinates = my_gaussian.get_coordinates_through_scf() # get the atom coordinates during the scf
spin_contamination = my_gaussian.get_spin_contamination() # get spin contamination data for the final optimized density
```

thermochemical analysis using:
    
```python
thermochemical_data = my_gaussian.get_thermo_chemistry() # get the thermochemical data from the log file
```

coordinate and bond analysis:

```python
coordinate_data = my_gaussian.get_coordinates() # get the coordinates of the atoms
bond_data = my_gaussian.get_bonds_from_log() # get the bonds and their lengths from bond information in the log file
bond_data = my_gaussian.get_bonds_from_coordinates(cutoff = 1.8) # get the bonds from the coordinates using a cuttoff
my_gaussian.get_optimisation_trajectory('opt_traj.pdb') # write the optimisation trajectory to a pdb file
```

charge and spin analysis:

```python
charge_data = my_gaussian.get_mulliken_charges() # get the mulliken charges for each atom
spin_data = my_gaussian.get_mulliken_spin_densities() # get the spin density for each atom
charge_data = my_gaussian.get_esp_charges() # get the esp charges for each atom
```

frequency analysis:

```python
frequency_data = my_gaussian.get_raman_frequencies() # get the raman frequencies
raman_spectra = my_gaussian.get_raman_spectra() # get the raman spectra for the system
```

or orbital analysis:

```python 
dos_plot = my_gaussian.get_dos_plot(height=800, width=800) # get the density of states plot
```
