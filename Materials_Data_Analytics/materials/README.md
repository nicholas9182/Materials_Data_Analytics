# Materials

The materials module contains classes for the construction of materials objects. These objects can be used to store information about materials, such as their composition, structure, properties, redox potentials etc. At the moment there are six main types of materials:

- **ions**
- **solutes**
- **solvents**
- **electrolytes**

These objects can be created and parsed into other classes for the construction of, for example, a microkinetic model. 

## Ions

### Authors

 - Dr. Nicholas Siemons ([GitHub](https://github.com/nicholas9182)) (nsiemons@stanford.edu)

### Getting Started

Ion objects can be created by passing in the name of the ion. For example, to create a sodium ion object:

```python
from Materials_Data_Analytics.materials.ions import Cation
from Materials_Data_Analytics.materials.ions import Anion

na_cation = Cation(name='Na+')
cl_anion = Anion(name='Cl-')
```

Most common ions are currently supported for construction in this way. Alternatively, the ion object can be created by passing in the name of the ion and the charge of the ion using the custom constructor class. For example, to create a sodium ion object:

```python
custom_anion = Anion.from_custom_inputs(name='Pt-', charge=-1, formula='Pt')
```

### Attributes
Ions have various attributes, including:
 - ```name``` - the name of the ion
 - ```charge``` - the charge of the ion
 - ```formula``` - the chemical formula of the ion

### Methods
Currently ions have no public methods. 

<br>

## Solutes

### Authors

 - Dr. Nicholas Siemons ([GitHub](https://github.com/nicholas9182)) (nsiemons@stanford.edu)

### Getting Started

Solute objects can be created by passing in the name of the solute. For example, to represent molecular oxygen solute:

```python
from Materials_Data_Analytics.materials.solutes import Solute

my_solute = Solute(name='H2O2')
```

Some common solutes for aqeous electrochemistry are currently supported for construction in this way. Oxygen is represented by its own class, which can be created by passing in the name of the solute. For example, to create a molecular oxygen solute object:

```python
from Materials_Data_Analytics.materials.solutes import MolecularOxygen

o2_solute = MolecularOxygen()
```

### Attributes
Solutes have the following attributes:
 - ```name``` - the name of the solute
 - ```formal_reduction_potentials``` - the formal reduction potentials of the solute
 - ```formula``` - the chemical formula of the solute

### Methods
Currently solutes have no public methods.

<br>

## Solvents

### Authors

 - Dr. Nicholas Siemons ([GitHub](https://github.com/nicholas9182)) (nsiemons@stanford.edu)

### Getting Started

Solvent objects can be created by passing in the name of the solvent. For example, to represent water solvent:

```python
from Materials_Data_Analytics.materials.solvents import Solvent

water_solvent = Solvent('water')
```

Currently most common solvents are supported for construction in this way, including water, methanol, ethanol, acetone, acetonitrile, dimethylsulfoxide, and dimethylformamide. Equally a custom solvent can be created, 

```python
custom_solvent = Solvent.from_custom_inputs(name='EMIM-TFSI', formula='C8H15F6N2O4S2', pH=6)
```

### Attributes
Solvents have the following attributes:
 - ```name``` - the name of the solvent
 - ```formula``` - the chemical formula of the solvent
 - ```pH``` - the pH of the solvent

### Methods
Currently solvents have no public methods.

<br>

## Electrolytes

### Authors

 - Dr. Nicholas Siemons ([GitHub](https://github.com/nicholas9182)) (nsiemons@stanford.edu)

### Getting Started

Electrolyte objects can be created by passing ions, solutes, and solvents into the constructor along with their concentrations as a dictionary. For example:

```python
from Materials_Data_Analytics.materials.electrolytes import Electrolyte

na_cation = Cation(name='Na+')
cl_anion = Anion(name='Cl-')
water_solvent = Solvent('water')
oxygen_solute = MolecularOxygen()

my_electrolyte = Electrolyte(solvent=water_solvent, 
                             cation=na_cation, 
                             anion=cl_anion, 
                             concentrations={na_cation: 0.1, cl_anion: 0.1, oxygen_solute: 0.0008}, 
                             solute=oxygen_solute, 
                             pH=14.2
                             )

```

### Attributes

Electrolytes have the following attributes:
 - ```solvent``` - the solvent of the electrolyte
 - ```cation``` - the cation of the electrolyte
 - ```anion``` - the anion of the electrolyte
 - ```concentrations``` - a dictionary of the concentrations of the ions and solutes in the electrolyte
 - ```solute``` - the solute of the electrolyte
 - ```pH``` - the pH of the electrolyte

### Methods

Currently electrolytes have no public methods.
