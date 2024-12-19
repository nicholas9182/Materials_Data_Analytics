# Experiment Modelling

The experiment_modelling module provides frameworks for analysing experimental data. Currently, the module supports cyclic voltammetry data, and X-ray diffraction data. The general structure of the module is to create a class for each type of measurement you may take. The measurements can optionally be given metadata. The metadata could, for example, contain an independant variable of an experiment, allowing for the construction of an experiment from multiple measurements.

## Cyclic Voltammetry 

### Authors

 - Dr. Nicholas Siemons ([GitHub](https://github.com/nicholas9182)) (nsiemons@stanford.edu)

### Getting Started

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
    my_cv = CyclicVoltammogram.from_aftermath('path/to/file.csv', scan_rate=1500)
    ```

Note that the CycleVoltammogram object contains current in the units of mA, and potential in V against the reference used in the experiment. 

### Attributes

The CyclicVoltammogram has various important attributes, including 
 - ```data``` - a pandas dataframe containing the potential, current, time, cycle, and direction data
 - ```pH``` - the pH of the solution
 - ```temperature``` - the temperature of the solution
 - ```cation``` - the cation used in the solution
 - ```anion``` - the anion used in the solution
 - ```steps_per_cycle``` - the number of voltage steps per cycle

### Methods

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

## X-ray diffraction

### Authors

 - Dr. Arianna Magni ([GitHub](https://github.com/magaris)) (amagni@stanford.edu)
 - Dr. Nicholas Siemons ([GitHub](https://github.com/nicholas9182)) (nsiemons@stanford.edu)

To be written

