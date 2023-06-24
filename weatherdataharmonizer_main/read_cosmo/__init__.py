#   @author: Michael Wagner
#   @organization: TU Dresden
#   @contact: michael.wagner@tu-dresden.de
"""
The read_cosmo package provides access to Cosmo formatted data from the german weather service (DWD).

The module readCosmo.py provides a function to read in typical Cosmo data, either as bz2 compressed or uncompressed grib
data. In the case of compressed data a temporary file is written (in OS specific temporary folder) and deleted
afterwards. Further a function is incorporated to return coordinates (lon/lat).

Metadata module has Metadata class with relevant metadata of Cosmo formatted data.

CosmoData module has CosmoData class with relevant data of Cosmo formatted data.
"""