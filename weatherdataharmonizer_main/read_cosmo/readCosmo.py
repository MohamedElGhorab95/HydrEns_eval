#   @author: Michael Wagner
#   @organization: TU Dresden
#   @contact: michael.wagner@tu-dresden.de

import os
import numpy as np
import pygrib
import datetime as dt

from met_entities.LonLatTime import LonLatTime
from met_entities.VariableDescription import DataDescription
from read_cosmo.CosmoData import CosmoData
from read_cosmo.Metadata import Metadata
from aux_tools import grib2_tools as gt


def read_cosmo_d2(filename, scale_factor=1, fill_value=-999, variant='d2', short=None):
    """
    Read in CosmoD2 data from file, either as grib2 file or as bz2 compressed grib2 file. Currently, both variants,
    CosmoD2 and CosmoD2EPS are supported.

    :param filename: name of file
    :type filename: str
    :param scale_factor: the final data has to be multiplied with this value
    :type scale_factor: float, optional
    :param fill_value: missing data is filled with that value
    :type fill_value: float, optional
    :param variant: specify variant, either 'd2' (default) or 'd2eps'
    :type variant: str, optional
    :param short: if int16 or int32, the large data variables are cast to numpy.int16/int32 to minimize memory
        usage; the user must pay attention to the scale_factor to include the necessary precision
    :type short: str, optional
    :return: Cosmo data of four 15 min values as in the content of the grib2 file
    :rtype: list
    """
    if short and scale_factor >= 1:
        print(f'the attribute short casts the data to 16 bit signed integer for memory saving; the scale_factor = '
              f'{scale_factor} imposes a low precision')
    if short and np.isnan(fill_value):
        raise Exception(f'if short data type is requested, the fill_value must not be nan')
    if short and (short != 'int16' and short != 'int32'):
        raise Exception(f'datatype {short} not supported; please change to int16 or int32')

    # read data from cosmo d2 grib2 file
    cosmo_data_realis = []  # list of realisations with lists of CosmoData() objects - only used at variant='eps'
    cosmo_data = []  # list of CosmoData() objects

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
            metadata.datum_end_of_overall_time_interval = dt.datetime(grb['yearOfEndOfOverallTimeInterval'],
                                                                      grb['monthOfEndOfOverallTimeInterval'],
                                                                      grb['dayOfEndOfOverallTimeInterval'],
                                                                      grb['hourOfEndOfOverallTimeInterval'],
                                                                      grb['minuteOfEndOfOverallTimeInterval'])
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
            idx_missing = np.where(data == metadata.missing_value)
            data = data / scale_factor
            data[idx_missing] = fill_value

            if short == 'int16':
                cosmo_data.append(CosmoData(data=np.short(np.around(data)), metadata=metadata))
            elif short == 'int32':
                cosmo_data.append(CosmoData(data=np.int32(np.around(data)), metadata=metadata))
            else:
                cosmo_data.append(CosmoData(data=data, metadata=metadata))
            if variant.lower() == 'd2eps' and grbs.messages == 80 and ct == 4:
                    cosmo_data_realis.append(cosmo_data)
                    cosmo_data = []
                    ct = 0
            elif variant.lower() == 'd2eps' and grbs.messages == 20:
                    cosmo_data_realis.append(cosmo_data)
                    cosmo_data = []

    grbs.close()

    if idx_bz2:
            os.remove(filename)

    if variant.lower() == 'd2':
            return cosmo_data
    elif variant.lower() == 'd2eps':
            return cosmo_data_realis


def get_lonlat(corners=False):
        """
        Calculate longitudes and latitudes for centerpoints of the last Cosmo-D2 rasters (716 rows, 651 columns).

        :param corners: if true, the coordinates of all four corner points are returned (default: False)
        :type corners: bool, optional
        :return: longitude and latitude objects
        :rtype: met_entities.LonLatTime.LonLatTime, met_entities.LonLatTime.LonLatTime
        """
        lon_n = -170 * np.pi / 180
        lat_n = 40 * np.pi / 180

        # due to unclear cosmo d2 conversion equation from DWD rely on adapted python code from another source:
        # https://gis.stackexchange.com/questions/10808/manually-transforming-rotated-lat-lon-to-regular-lat-lon
        lon_s = lon_n * 180 /np.pi - 180
        lat_s = -lat_n * 180 / np.pi

        phi = lon_s * np.pi / 180
        theta = (90 + lat_s) * np.pi / 180

        lon_corner = np.array([-7.5, 5.5]) * np.pi / 180
        lat_corner = np.array([-6.3, 8.0]) * np.pi / 180
        resolution = 0.02 * np.pi / 180

        if corners:
            lon_r = lon_corner
            lat_r = lat_corner
        else:
            lon_r = np.arange(lon_corner[0], lon_corner[1] + resolution, resolution)
            lat_r = np.arange(lat_corner[0], lat_corner[1] + resolution, resolution)

        lon_r_mat, lat_r_mat = np.meshgrid(lon_r, lat_r)

        x = np.cos(lon_r_mat) * np.cos(lat_r_mat)
        y = np.sin(lon_r_mat) * np.cos(lat_r_mat)
        z = np.sin(lat_r_mat)
        phi = -phi
        theta = -theta
        x_new = np.cos(theta) * np.cos(phi) * x + np.sin(phi) * y + np.sin(theta) * np.cos(phi) * z
        y_new = -np.cos(theta) * np.sin(phi) * x + np.cos(phi) * y - np.sin(theta) * np.sin(phi) * z
        z_new = -np.sin(theta) * x + np.cos(theta) * z

        lon_g = np.flipud(np.arctan2(y_new, x_new)) * 180 / np.pi
        lat_g = np.flipud(np.arcsin(z_new)) * 180 / np.pi

        lon_description = DataDescription(units='degrees_east',
                                          long_name='longitude of center',
                                          coordinate_system='WGS84, EPSG:4326')
        lon = LonLatTime(data=lon_g, data_description=lon_description)

        lat_description = DataDescription(units='degrees_north',
                                          long_name='latitude of center',
                                          coordinate_system='WGS84, EPSG:4326')
        lat = LonLatTime(data=lat_g, data_description=lat_description)

        return lon, lat
