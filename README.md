# HydrEns_eval
The evaluation tools package provides a set of methods to evaluate the performance and reliability of rainfall/runoff forecasts
against radar/gauge observations.
It also works with ensemble datasets, and provides visual comparisons and contingency table based metrics for either grid or array data. 
# Dependencies
        * Weather data Harmonizer package developed by Michael Wagner at TUDresden michael.wagner@tu-dresden.de
			---- can be retrieved at https://gitlab.hrz.tu-chemnitz.de/ihmtud/weatherdataharmonizer 
        * xarray package            V 2023.4.2 
        * xskillscore package       V 0.00.24
        * geopandas package         V 0.12.02
        * pyshp package             V 2.03.01
        * pyproj package            V 3.05.00
        * odc-geo package           V 0.03.03
        * rioxarray package         V 0.14.01
        * pandas                    V 1.04.00
# Structure
![HydrEns_eval (4)](https://github.com/MohamedElGhorab95/HydrEns_eval/assets/97175071/da684461-5e5b-4bbf-8155-47abb7a2c7c4)





# fr_to_netcdf_tools.py
This module includes all the functions to convert from the .grib files by the DWD to netCDF format, to be later used as input for the other modules.
# fr_entities_tools.py
This module includes all the forecast and observation objects and methods to be used for evaluation.
# fr_Conttools.py
This module includes the contingency object that utilizes the forecast and observation objects to perform the comparison. 
# fr_ROCtools.py
This module includes the relative operation characteristics curve object and plotting method.
# License
HydrEns_eval is licensed under Apache-2.0. You may obtain a copy
of the License here in LICENSE or at http://www.apache.org/licenses/LICENSE-2.0.