# HydrEns_eval
The evaluation tools package provides a set of methods to evaluate the performance and reliability of rainfall/runoff forecasts
against radar/gauge observations.
It also works with ensemble datasets, and provides visual comparisons and contingency table based metrics for either grid or array data. 
# Dependencies
        * Weather data Harmonizer package developed by Michael Wagner at TUDresden michael.wagner@tu-dresden.de
        * xarray package            V 2023.4.2 
        * xskillscore package       V 0.00.24
        * geopandas package         V 0.12.02
        * pyshp package             V 2.03.01
        * pyproj package            V 3.05.00
        * odc-geo package           V 0.03.03
        * rioxarray package         V 0.14.01
        * pandas                    V 1.04.00
![HydrEns_eval](https://user-images.githubusercontent.com/97175071/235663636-2d083277-5818-40f8-a775-cf79f64ec4b1.png)
# Structure
# fr_to_netcdf_tools.py
This module includes all the functions to convert between the .grib files by the DWD to netCDF format, to be later used as input for the other modules.
# fr_entities_tools.py
This module includes all the forecast and observation objects and methods to be used for evaluation.
# fr_Conttools.py
This module includes the contingency object that utilizes the forecast and observation objects to perform the comparison. 
# fr_ROCtools.py
This module includes the relative operation characteristics curve object and plotting method.
