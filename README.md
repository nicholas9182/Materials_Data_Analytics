# Materials_Data_Analytics

A python package for the handling and analysis of a wide range of synthetic and experimental data for the development of next generation energy materials.

## Authors

 - Dr. Nicholas Siemons (nsiemons@stanford.edu)
 - Srikant Sagireddy (srikant@stanford.edu)

## Philosophy

The package is designed to be as user-friendly as possible, with a focus on ease of use, readability and distributility. The package is designed to be as modular as possible, with the aim of allowing users to easily extend the package to suit their own needs. Wherever suitable, the code has been written so as to be method chainable, allowing for more concise and readable code. 

Generally any class methods will do one of four things - 
 - modify the self of the object in place
 - return a pandas dataframe, which can then be method chained with the usual pandas methods
 - return a plotly.express figure, which can then be modified with the usual plotly methods. In these cases, arguments can be passed to the method to modify the figure according with the plotly documentation through the use of **kwargs.
 - display a plotly figure

 Additionally, the creation of almost all objects can be done by parsing a metadata dictionary to the object. This means that, for example, measurements corresponding to different systems can easily be compared by calculating their properties along with their metadata into a long-format dataframe, and then comparing those properties using the usual pandas and plotly methods.

## Key Analysis Types

- **Cyclic Voltammetry**

- **Gaussian Quantum Chemistry**
- **Gromacs+Plumed Metadynamics** 

## Installation

To install the package, clone the repository and run the following command:

```sh
pip install ./path/to/Materials_Data_Analytics
```

or alternatively to install the most recent version from PyPi, run 

```sh
pip install Materials_Data_Analytics
```

<br><br>

# Functionality 

## Experiment Modelling

### Cyclic Voltammetry 

Currently there are three supported ways of initiating a CyclicVoltammogram object:
 - From three iterables of equal length, representing the voltage, current and time data. 
    ```python
    from Materials_Data_Analytics.experiment_modelling.cyclic_voltammetry import CyclicVoltammogram
    my_cv = CyclicVoltammogram(potential=voltage_data, current=current_data, time=time_data)
    ```
 - From a .txt file exported from biologic software
    ```python
    my_cv = CyclicVoltammogram.from_biologic('path/to/file.txt')
    ```

 - From a .csv file exported from aftermath software
    ```python
    my_cv = CyclicVoltammogram.from_aftermath('path/to/file.csv')
    ```

Note that the CycleVoltammogram object contains current in the units of mA, and potential in V against the reference used in the experiment. 

The CyclicVoltammogram has various important attributes, including 
 - ```data``` - a pandas dataframe containing the potential, current, time, cycle, and direction data
 - ```pH``` - the pH of the solution
 - ```temperature``` - the temperature of the solution
 - ```cation``` - the cation used in the solution
 - ```anion``` - the anion used in the solution
 - ```steps_per_cycle``` - the number of voltage steps per cycle


Once created, various self-modifying operations can be performed on it. These include dowsampling and cycle selecting.

```python
my_cv.downsample(200) # Downsample the data to 200 points per cycle
my_cv.drop_cycles(drop=[1, 2, 3]) # Drop the first three cycles
my_cv.drop_cycles(keep=[1, 2, 3]) # Keep only the first three cycles
```

Once created, basic plots can be generated using 
    
```python
figure = mv_cv.get_current_potential_plot(height=800, width=800)
figure = mv_cv.get_current_time_plot(height=800, width=800)
figure = mv_cv.get_potential_time_plot(height=800, width=800)
```

These figures can be modified using the usual plotly methods.

Analysis can be performed on the amount of charge passing in and out of the system, as well as the peak current and potential values. 

```python
charge_passed = my_cv.get_charge_passed() # get the charges passed per cycle
max_charge_passed = my_cv.get_max_charge_passed() # get the maximum charge passed in a cycle
peak_data = my_cv.get_peaks() # get the peak current and potential data
```

with plots alongside, such as 

```python
figure = my_cv.get_charge_passed_plot(height=800, width=800)
figure = my_cv.get_max_charge_passed_plot(height=800, width=800)
figure = my_cv.get_charge_integration_plot(cycle=2, direction='reduction', height=800, width=800)
figure = my_cv.get_maximum_charge_integration_plot(section=2, height=800, width=800)
figure = my_cv.get_peak_plot(height=800, width=800)
```

<br><br>

### X-ray diffraction

To be written


<br><br>

## Quantum Chemical Simulations

### Gaussian log file analysis

A GaussianParser object can be created either from a gaussian log output file, or from a list of output files (if, for example calculation restarts are rquired). At the moment it is assumed that the flag #p was included in the keywords line of the gaussian input file, indicating the log file to be more verbose in its output. 

```python
from Materials_Data_Analytics.quantum_chemistry.gaussian import GaussianParser
my_gaussian = GaussianParser('path/to/logfile.log')
my_gaussian = GaussianParser(['path/to/logfile1.log', 'path/to/logfile2.log'])
```

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
orbital_data = my_gaussian.get_orbitals() # get the orbital energies as a pandas dataframe
dos_plot = my_gaussian.get_dos_plot(height=800, width=800) # get the density of states plot
```

<br><br>

## Atomistic Simulations

### MetaDynamics Analysis

Metadynamcs analysis comprises four main classes - 

 - MetaTrajectory.  The MetaTrajectory represents the trajectory of your simulation in a space defined by your collective variables.  It can be created from a colvar file produced from a metadynamics simulation in gromacs with plumed. 

```python
from Materials_Data_Analytics.metadynamics.free_energy import MetaTrajectory
my_traj = MetaTrajectory(colvar_file='path/to/colvar')
```

The MetaTrajectory will automatically compute the weighting for each from from the biased ensemble to perform reweighting into the unbiased ensemble. The trajectory data can then be obtained using 

```python
data = my_traj.get_data()
```

and has the following attributes:
- ```walker``` - the walker number if doing multiwalker calculations
- ```cvs``` - the collective variables used in the simulation
- ```temperature``` - the temperature of the simulation
- ```opes``` - A boolean indicating whether OPES was used

The next two classes are both part of the more general FreeEnergyShape class. The FreeEnergyLine and FreeEnergySurface represent free energy profiles over either one or two variables. The FreeEnergyShape has the following attributes:

- ```cvs``` - the collective variables used in the simulation
- ```temperature``` - the temperature of the simulation
- ```dimension``` - the dimension of the free energy profile (1 for a line or 2 for a surface). 

They can be created either from a pandas dataframe with the cv and free energy data, or from a plumed file - 

```python
from Materials_Data_Analytics.metadynamics.free_energy import FreeEnergyLine, FreeEnergySurface
my_line = FreeEnergyLine(data=data) # create a line from a pandas dataframe
my_line = FreeEnergyLine.from_plumed('path/to/1DFES.dat') # create a line from a plumed file
my_surface = FreeEnergySurface(data=data) # create a surface from a pandas dataframe
my_surface = FreeEnergySurface.from_plumed('path/to/2DFES.dat') # create a surface from a plumed file
```

Finally, the FreeEnergySpace is a class which holds all the trajectories and shapes so as to perform free energy analysis on the system. The free energy space has the following attributes:

- ```trajectories``` - a list of MetaTrajectory objects
- ```lines``` - a list of FreeEnergyLine objects
- ```surfaces``` - a list of FreeEnergySurface objects
- ```n_walker``` - the number of walkers in the simulation
- ```sigmas``` - the hill widths used in the simulation
- ```hills``` - a dataframe with the hills deposited in the simulation
- ```n_timesteps``` - the number of timesteps in the simulation
- ```max_time``` - the maximum time in the simulation
- ```dt``` - the timestep of the simulation
- ```cvs``` - the collective variables used in the simulation
- ```opes``` - a boolean indicating whether OPES was used
- ```biasexchange``` - a boolean indicating whether bias exchange was used
- ```temperature``` - the temperature of the simulation.

To create a FreeEnergySpace object, you can either create an empty one and add trajectories, lines and surfaces, or you can create one with a HILLS file obtained from a metadynamics trajectory:

```python
from Materials_Data_Analytics.metadynamics.free_energy import FreeEnergySpace
my_space = FreeEnergySpace(hills_file='path/to/HILLS') # create a space from a HILLS file
my_space = FreeEnergySpace() # create an empty space
```

It can then be populated with trajectories, lines and surfaces:

```python
my_space.add_trajectory(my_traj)
my_space.add_line(my_line)
my_space.add_surface(my_surface)
```

Alternatively, you can create a FreeEnergySpace object from a directory where a metadynamics simulation was run. If you do this, then the FES's must be computed with the plumed sum_hills tool and the fes for each CV or pair of CVs must be in its own subdirectory. Furthermore, the MD files for each walker must be in their own subdirectory.  

```python
my_space = FreeEnergySpace.from_standard_directory('path/to/directory')
```

Once created, the FreeEnergySpace object can be used to perform various analysis on the system. This includes reweighting the trajectories to the unbiased ensemble, and calculating the free energy profiles. For example. visualising the hills deposited in a system:

```python
my_space = FreeEnergySpace(hills_file='path/to/HILLS')
data = my_space.get_long_hills() # get the hills data in long format
data = my_space.get_hills_average_across_walkers() # get the hills data averaged across the walkers
data = my_space.get_hills_max_across_walkers() # get the hills data with the maximum height
figures = my_space.get_hills_figures(height=800, width=800) # get a hills figure for each walker
figure = my_space.get_average_hills_figure(height=800, width=800) # get the hills figure averaged across the walkers
figure = my_space.get_max_hills_figure(height=800, width=800) # get the hills figure with the maximum height
```

Alternatively, new FreeEnergyLines and FreeEnergySurfaces can be created from the space through reweighting the trajectories:

```python
new_line = my_space.get_reweighted_line(cv='cv', bins=100) # get a reweighted line
new_surface = my_space.get_reweighted_surface(cvs=['cv1','cv2'], bins=100) # get a reweighted surface
new_line = my_space.get_reweighted_line_with_walker_error(cv='cv', bins=100) # get the reweighted line with errors as deviation across the walkers
```


Enjoy!

If you find this code useful, please cite it according to the CITATION.cff file. 