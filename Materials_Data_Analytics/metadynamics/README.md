# Atomistic Simulations

The atomistic simulations module contains classes for the analysis of atomistic metadynamics simulations. Currently it is able to analyse data generated from gromacs simulations with plumed.

## Metatrajectory

### Authors

 - Dr. Nicholas Siemons ([GitHub](https://github.com/nicholas9182)) (nsiemons@stanford.edu)

### Getting Started

The Metatrajectory class is used to analyse the output of a metadynamics simulation. It can be created from a COLVAR file output from plumed and represents the trajectory in a space defined by the collective variables used in the simulation. 

```python
from Materials_Data_Analytics.metadynamics.free_energy import MetaTrajectory

my_traj = MetaTrajectory(colvar_file='path/to/colvar')
```

The MetaTrajectory will automatically compute the weighting for each from from the biased ensemble to perform reweighting into the unbiased ensemble. The trajectory data can then be obtained using 

```python
data = my_traj.get_data()
```

### Attributes
A metatrajectory has the following attributes:
- ```walker``` - the walker number if doing multiwalker calculations
- ```cvs``` - the collective variables used in the simulation
- ```temperature``` - the temperature of the simulation
- ```opes``` - A boolean indicating whether OPES was used

### Methods
The metatrajectory has the following methods:
- ```get_data()``` - returns the trajectory data

<br>

## FreeEnergyLine, FreeEnergySurface, FreeEnergyShape

### Authors

 - Dr. Nicholas Siemons ([GitHub](https://github.com/nicholas9182)) (nsiemons@stanford.edu)

### Getting Started

The FreeEnergyLine and FreeEnergySurface class are both child classes of the more general FreeEnergyShape class. The FreeEnergyLine and FreeEnergySurface represent free energy profiles over either one or two variables. 

They can be created either from a pandas dataframe with the cv and free energy data, or from a plumed file - 

```python
from Materials_Data_Analytics.metadynamics.free_energy import FreeEnergyLine, FreeEnergySurface

my_line = FreeEnergyLine(data=data) # create a line from a pandas dataframe
my_line = FreeEnergyLine.from_plumed('path/to/1DFES.dat') # create a line from a plumed file
my_surface = FreeEnergySurface(data=data) # create a surface from a pandas dataframe
my_surface = FreeEnergySurface.from_plumed('path/to/2DFES.dat') # create a surface from a plumed file
```

### Attributes
FreeEnergyLine and FreeEnergySurface have the following attributes:
- ```cvs``` - the collective variables used in the simulation
- ```temperature``` - the temperature of the simulation
- ```dimension``` - the dimension of the free energy profile (1 for a line or 2 for a surface). 
- ```metadata``` - a dictionary containing any additional metadata

### Methods
The FreeEnergyLine and FreeEnergySurface have the following methods:
- ```get_data()``` - returns the free energy data
- ```get_nearest_value()``` - returns the free energy value at the nearest point to a given value
- ```set_datum()``` - sets the free energy value to zero at a given point
- ```get_time_difference()``` - returns the energy difference between two points on a FES as a function of simulation time 
- ```set_errors_from_time_dynamics()``` - sets the errors on the FES from the time dynamics of the simulation
- ```set_as_symmetric()``` - sets the FES as symmetric about a given symmetry line

<br>

## FreeEnergySpace

### Authors

 - Dr. Nicholas Siemons ([GitHub](https://github.com/nicholas9182)) (nsiemons@stanford.edu)

### Getting Started

Finally, the FreeEnergySpace is a class which holds all the trajectories and shapes so as to perform free energy analysis on the system. The free energy space has the following attributes:

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

### Attributes
The FreeEnergySpace has the following attributes:
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

### Methods
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
