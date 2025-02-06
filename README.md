# Materials_Data_Analytics

A python package for the handling and analysis of a wide range of synthetic and experimental data for the development of next-generation energy materials. 
Its modular, method-chainable design makes it easy for researchers to analyze complex datasets efficiently, promoting reproducibility and extensibility.

## Authors

 - Dr. Nicholas Siemons ([GitHub](https://github.com/nicholas9182)) (nsiemons@stanford.edu)
 - Dr. Arianna Magni ([GitHub](https://github.com/magaris)) (amagni@stanford.edu)
 - Srikant Sagireddy ([GitHub](https://github.com/sbsagireddy)) (srikant@stanford.edu)

## Why

The package was developed to streamline the analysis of a wide range of materials data, with a focus on the development of next-generation energy materials. The package is designed to be as user-friendly as possible, with a focus on ease of use, readability and distributility. The package is designed to be as modular as possible, with the aim of allowing users to easily extend the package to suit their own needs. 

This package looks to keep many aspects of materials data analysis in one place. In doing so it allows for complex analysis to be done in a single environment, and for the results of that analysis to be easily compared. Furthermore, it will allow for analysis involving data from a variety of sources. 

## Philosophy

The package is designed to be as user-friendly as possible, with a focus on ease of use, readability and distributility. The package is designed to be as modular as possible, with the aim of allowing users to easily extend the package to suit their own needs. Wherever suitable, the code has been written so as to be method chainable, allowing for more concise and readable code. 

Generally any class methods will do one of four things - 
 - modify the self of the object in place
 - return a pandas dataframe, which can then be method chained with the usual pandas methods
 - return a plotly.express figure, which can then be modified with the usual plotly methods. In these cases, arguments can be passed to the method to modify the figure according with the plotly documentation through the use of **kwargs.
 - display a plotly figure

 Additionally, the creation of almost all objects can be done by parsing a metadata dictionary to the object. This means that, for example, measurements corresponding to different systems can easily be compared by calculating their properties along with their metadata into a long-format dataframe, and then comparing those properties using the usual pandas and plotly methods. Finally, internally this package leverages the power of pandas and plotly for the handling and visualization of data.

## Key Analysis Types

- **Cyclic Voltammetry**: Analyze and visualize electrochemical measurements for battery and fuel cell development.
- **GIWAXS Analysis**: Interpret grazing-incidence wide-angle X-ray scattering data to understand material structures.
- **Gaussian Quantum Chemistry**: Analyze and visualize results from quantum chemistry calculations.
- **Gromacs+Plumed Metadynamics**: Analyze and visualize results from metadynamics simulations. 

## Installation

To install the package, clone the repository and run the following command:

```sh
pip install ./path/to/Materials_Data_Analytics
```

or alternatively to install the most recent version from PyPi, run 

```sh
pip install Materials_Data_Analytics
```

For usage instructions, see the README.md files in the module folders.

#### **Dependencies**

- scipy
- pandas
- plotly
- matplotlib
- typer
- click
- numpy
- networkx
- MDAnalysis
- dash
- kaleido
- pyFAI
- pygix
- Datetime

## License

This project is licensed under the MIT License - see the LICENSE.md file for details.
