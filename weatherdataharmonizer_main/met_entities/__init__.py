#   @author: Michael Wagner
#   @organization: TU Dresden
#   @contact: michael.wagner@tu-dresden.de
"""
The package met_entities consists of different classes to handle different weather data formats of the german weather
service (DWD). It supports RadolanRW, RadvorRQ, IconD2, IconD2EPS datasets.

MetEntities.py module contains MetEntities class as an abstract class that defines methods for handling
meteorological data.

RadolanRW.py module contains RadolanRW class that provides all relevant data of RadolanRW data from DWD.

RadvorRQ.py: dito for RadvorRQ data

RadolanRV.py: dito for RadolanRV data

IconD2.py: dito for Icon-D2 data

IconD2EPS: dito for Icon-D2-EPS data

IconEU.py: dito for Icon-EU data

IconEUEPS.py: dito for Icon-EU-EPS data

CosmoD2.py: dito for Cosmo-D2 data

CosmoD2EPS.py: dito for Cosmo-D2-EPS data

GeoReferencedData.py module contains GeoReferencedData class that serves as a container for spatial data including
coordinates, data and their description; further elementary functions are implemented (regridding, cropping).

LonLatTime.py module contains LonLatTime class that is used as a container for coordinate- or time-related data to
facilitate the export to netCDF.

VariableDescription.py module contains classes DataDescription (further description of climate data),
TimeDescription (description of time aspects of climate data), RegridDescription (attributes necessary for regridding),
CropDescription (attribute necessary for cropping), and NcDimVarDescription (describes the dimension and variable names
of a netcdf file in detail).

Exceptions.py module contains typical exceptions raised in met_entities.

WeatherData.py module contains WeatherData class that provides collection and harmonizing routines for all types of
above-mentioned supported data.
"""
