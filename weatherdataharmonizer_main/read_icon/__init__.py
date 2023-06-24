#   @author: Michael Wagner
#   @organization: TU Dresden
#   @contact: michael.wagner@tu-dresden.de
"""
The read_icon package provides access to Icon formatted data from the german weather service (DWD).

The module readIcon.py provides a function to read in typical Icon data, either as bz2 compressed or uncompressed grib
data. In the case of compressed data a temporary file is written (in OS specific temporary folder) and deleted
afterwards. Further a function is incorporated to return coordinates (lon/lat).

Metadata module has Metadata class with relevant metadata of Icon formatted data.

IconData module has IconData class with relevant data of Icon formatted data.
"""
