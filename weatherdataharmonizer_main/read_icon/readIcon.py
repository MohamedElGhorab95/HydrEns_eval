#   @author: Michael Wagner
#   @organization: TU Dresden
#   @contact: michael.wagner@tu-dresden.de
import importlib.resources
import os
import datetime as dt

import numpy as np
import pygrib

from met_entities.LonLatTime import LonLatTime
from met_entities.VariableDescription import DataDescription
from read_icon.IconData import IconData
from read_icon.Metadata import Metadata
from aux_tools import grib2_tools as gt
from met_entities.Exceptions import *


def read_icon_d2(filename, scale_factor=1, fill_value=-999, variant='d2', short=None):
    """
    Read in IconD2 or IconEU data from file, either as grib2 file or as bz2 compressed grib2 file. Currently, IconD2,
    IconD2EPS, IconEU, and IconEUEPS are supported.

    :param filename: name of file
    :type filename: str
    :param scale_factor: the final data has to be multiplied with this value
    :type scale_factor: float, optional
    :param fill_value: missing data is filled with that value
    :type fill_value: float, optional
    :param variant: specify variant, either 'd2' (default), 'd2eps', 'eu', or 'eueps'
    :type variant: str, optional
    :param short: if int16 or int32, the large data variables are cast to numpy.int16/int32 to minimize memory
        usage; the user must pay attention to the scale_factor to include the necessary precision
    :type short: str, optional
    :return: Icon data of four 15 min values (for some older data also one 1 h value possible) as in the content of the
        grib2 file
    :rtype: list
    """
    filename_orig = filename
    if short and scale_factor >= 1:
        print(f'the attribute short casts the data to 16 bit signed integer for memory saving; the scale_factor = '
              f'{scale_factor} imposes a low precision')
    if short and np.isnan(fill_value):
        raise Exception(f'if short data type is requested, the fill_value must not be nan')
    if short and (short != 'int16' and short != 'int32'):
        raise Exception(f'datatype {short} not supported; please change to int16 or int32')

    # read data from icon grib2 file
    icon_data_realis = []  # list of realisations with lists of IconData() objects - only used at variant='eps'
    icon_data = []  # list of IconData() objects

    idx_bz2 = False
    if filename[-4:] == '.bz2':
        # file is bz2 compressed
        idx_bz2 = True
        filename = gt.write_temporary_grib_file(filename)

    grbs = pygrib.open(filename)
    ct = 0
    for grb in grbs:
        ct = ct + 1

        # get metadata
        metadata = Metadata()
        metadata.centre = grb['centre']
        metadata.centre_description = grb['centreDescription']
        metadata.datum = dt.datetime(grb['year'], grb['month'], grb['day'], grb['hour'], grb['minute'])
        if 'yearOfEndOfOverallTimeInterval' in grb.keys():
            # 15 min products do not contain useful content in grib message key 'forecastTIme' -> build information
            # from 'yearOfEndOfOverallTimeInterval' and similar
            metadata.datum_end_of_overall_time_interval = dt.datetime(grb['yearOfEndOfOverallTimeInterval'],
                                                                      grb['monthOfEndOfOverallTimeInterval'],
                                                                      grb['dayOfEndOfOverallTimeInterval'],
                                                                      grb['hourOfEndOfOverallTimeInterval'],
                                                                      grb['minuteOfEndOfOverallTimeInterval'])
        else:
            # currently only hourly products contain useful content in grib message key 'forecastTime'
            metadata.datum_end_of_overall_time_interval = metadata.datum + dt.timedelta(minutes=grb['forecastTime'])
        metadata.datum_iso = metadata.datum_end_of_overall_time_interval.replace(tzinfo=dt.timezone.utc).isoformat()
        metadata.prediction_time = metadata.datum_end_of_overall_time_interval - metadata.datum
        metadata.step_type = grb['stepType']
        metadata.grid_type = grb['gridType']
        metadata.name = grb['name']
        metadata.units = grb['units']
        metadata.missing_value = grb['missingValue']
        metadata.number_of_missing = grb['numberOfMissing']
        metadata.number_of_data_points = grb['numberOfDataPoints']
        metadata.fill_value = fill_value

        # get data and replace missing values with fill_value
        data = grb.values
        if variant.lower() == 'eu':
            data = np.flipud(data)
        idx_missing = np.where(data == metadata.missing_value)
        data = data / scale_factor
        data[idx_missing] = fill_value

        if short == 'int16':
            icon_data.append(IconData(data=np.short(np.around(data)), metadata=metadata))
        elif short == 'int32':
            icon_data.append(IconData(data=np.int32(np.around(data)), metadata=metadata))
        else:
            icon_data.append(IconData(data, metadata=metadata))
        if variant.lower() == 'd2eps' and grbs.messages == 80 and ct == 4:
            icon_data_realis.append(icon_data)
            icon_data = []
            ct = 0
        elif variant.lower() == 'd2eps' and grbs.messages == 20:
            icon_data_realis.append(icon_data)
            icon_data = []
        elif variant.lower() == 'd2eps' and grbs.messages not in [20, 80]:
            raise ForecastFileFormatError(f'for IconD2EPS the number of grib messages must be either 20 or 80\n'
                                          f'the current number of messages is {grbs.messages}')

    grbs.close()

    if idx_bz2:
        os.remove(filename)

    if variant.lower() in ['d2', 'eu', 'eueps']:
        return icon_data
    elif variant.lower() == 'd2eps':
        return icon_data_realis


def get_lonlat(variant, clon_file=None, clat_file=None):
    """
    Read longitudes and latitudes from clon/clat files (either grib2 files or bz2 compressed files).

    :param variant: variant of Icon data (case-insensitive), currently supported: Icon-D2 - 'd2', Icon-D2-EPS - 'd2eps',
        Icon-EU - 'eu', Icon-EU-EPS - 'eueps'
    :type variant: str
    :param clon_file: grib2 or bz2 compressed grib2 file with center longitudes of Icon; if not given the file from the
        resources directory in this package is used
    :type clon_file: str, optional
    :param clat_file: grib2 or bz2 compressed grib2 file with center latitudes of Icon; if not given the file from the
        resources directory in this package is used
    :type clat_file: str, optional
    :return: center longitudes and latitudes
    :rtype: met_entities.LonLatTime.LonLatTime, met_entities.LonLatTime.LonLatTime
    """
    if clon_file is None and clat_file is None:
        # import files from sub folder resources
        if variant.lower() == 'd2':
            clon_file = 'resources/icon-d2_germany_icosahedral_time-invariant_2022110800_000_0_clon.grib2'
            clat_file = 'resources/icon-d2_germany_icosahedral_time-invariant_2022110800_000_0_clat.grib2'
        elif variant.lower() == 'd2eps':
            clon_file = 'resources/icon-d2-eps_germany_icosahedral_time-invariant_2022110900_000_0_clon.grib2'
            clat_file = 'resources/icon-d2-eps_germany_icosahedral_time-invariant_2022110900_000_0_clat.grib2'
        elif variant.lower() == 'eu':
            clon_file = 'resources/icon-eu_europe_regular-lat-lon_time-invariant_2023032300_RLON.grib2.bz2'
            clat_file = 'resources/icon-eu_europe_regular-lat-lon_time-invariant_2023032300_RLAT.grib2.bz2'
        elif variant.lower() == 'eueps':
            clon_file = 'resources/icon-eu-eps_europe_icosahedral_time-invariant_2023032300_clon.grib2.bz2'
            clat_file = 'resources/icon-eu-eps_europe_icosahedral_time-invariant_2023032300_clat.grib2.bz2'
        else:
            raise Exception(f'Icon variant {variant} not supported yet')
        import read_icon
        resources_dir = str(importlib.resources.files(read_icon))
        clon_file_complete = os.path.join(resources_dir, clon_file)
        clat_file_complete = os.path.join(resources_dir, clat_file)
    else:
        # clon and clat files are given externally
        clon_file_complete = clon_file
        clat_file_complete = clat_file

    idx_bz2_clon = False
    if clon_file_complete[-4:] == '.bz2':
        # the given clon file is bz2 compressed
        idx_bz2_clon = True
        clon_file_complete = gt.write_temporary_grib_file(clon_file_complete)
    clon_data = pygrib.open(clon_file_complete).message(1).values
    clon_description = DataDescription(units='degrees_east', long_name='longitude of center',
                                       coordinate_system='WGS 84, EPSG:4326')
    clon = LonLatTime(data=clon_data, data_description=clon_description)

    idx_bz2_clat = False
    if clat_file_complete[-4:] == '.bz2':
        # the given clat file is bz2 compressed
        idx_bz2_clat = True
        clat_file_complete = gt.write_temporary_grib_file(clat_file_complete)
    clat_data = pygrib.open(clat_file_complete).message(1).values
    if variant.lower() == 'eu':
        clat_data = np.flipud(clat_data)
    clat_description = DataDescription(units='degrees_north', long_name='latitude of center',
                                       coordinate_system='WGS 84, EPSG:4326')
    clat = LonLatTime(data=clat_data, data_description=clat_description)

    if idx_bz2_clon:
        os.remove(clon_file_complete)
    if idx_bz2_clat:
        os.remove(clat_file_complete)

    return clon, clat
