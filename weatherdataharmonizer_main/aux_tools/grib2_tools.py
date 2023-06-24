#   @author: Michael Wagner
#   @organization: TU Dresden
#   @contact: michael.wagner@tu-dresden.de

import tempfile
import bz2
import pygrib
import os
import numpy as np
import copy

def write_temporary_grib_file(filename):
    """
    Write a temporary grib2 file from a bz2 compressed file. It is saved in OS dependent directory for temporary files
    and can be removed after use with the filename.

    :param filename: bz2 compressed grib file
    :type filename: str
    :return: filename of temporary file
    :rtype: str
    """
    with tempfile.NamedTemporaryFile(prefix='icon', delete=False) as tmp_file, bz2.BZ2File(filename, 'rb') as file:
        for data in iter(lambda: file.read(100 * 1024), b''):
            tmp_file.write(data)
    return tmp_file.name

def get_keys(filename, message_number=1):
    """
    Get keys and corresponding values of the first message in the grib2 file.

    :param filename: name of grib2 or bz2 compressed grib file
    :type filename: str
    :param message_number: number of grib message
    :type message_number: int, optional
    :return: keys and values
    :rtype: dict
    """
    idx_bz2 = False
    if filename[-4:] == '.bz2':
        # file is bz2 compressed
        idx_bz2 = True
        filename = write_temporary_grib_file(filename)

    grbs = pygrib.open(filename)
    grb = grbs.message(message_number)
    grb_dict = {}
    for key in grb.keys():
        try:
            grb_dict[key] = grb[key]
            # print(f'{key}: {grb[key]}')
        except RuntimeError:
            grb_dict[key] = 'na'
            # print(f'{key} not found')

    grbs.close()

    if idx_bz2:
        os.remove(filename)

    return grb_dict


def accum_to_instantaneous_flattened(data_list, fill_value):
    """
    Convert accumulated forecast data series, e.g. from IconD2, to a flattened list with instantaneous values. All 15
    minutes forecasts are treated equally. Occasionally occurring negative precipitation values are corrected to zero.

    :param data_list: list of accumulated data (e.g. IconData or CosmoData class)
    :type data_list: list
    :param fill_value: missing data is filled with that value
    :type fill_value: float
    :return: 1D list of forecast data content from grib files (all 15 minutes datasets in a row, e.g. in IconD2)
    :type: list
    """
    data_list_flattened = [quarter for hour in data_list for quarter in hour]  # flattened for each 15 min
    num_data = len(data_list_flattened)

    for quarter in range(num_data - 1, 0, -1):
        # calculate instantaneous values from accumulated
        # loop backwards through the list as the upper list element itself is changed in each step and must not
        # influence following steps
        if np.isnan(fill_value):
            idx_missing = np.isnan(data_list_flattened[quarter].data)
        else:
            idx_missing = data_list_flattened[quarter].data == fill_value
        # data_list_flattened[quarter].data[idx_missing] = np.nan
        data_list_flattened[quarter].data = data_list_flattened[quarter].data - \
                                                 data_list_flattened[quarter - 1].data

        # correct negative values (due to machine precision) to zero
        idx_negative = data_list_flattened[quarter].data < 0
        if np.any(idx_negative):
            # print(f'negative values were corrected to 0')
            data_list_flattened[quarter].data[idx_negative] = 0
        data_list_flattened[quarter].data[idx_missing] = fill_value

    return data_list_flattened


def average_to_instantaneous_flattened(data_list, fill_value, start_time_step):
    """
    Convert averaged forecast data series, e.g. from IconD2, to a flattened list with instantaneous values. All 15
    minutes forecasts are treated equally. Occasionally occurring negative radiation values are corrected to zero. The
    de-averaging is built like proposed in the Icon description (Reinert et al.: DWD Database Reference for Global and
    Regional ICON and ICON-EPS Forecasting System. p. 47, Version 2.2.0, 2022).

    :param data_list: list of accumulated data (e.g. IconData class)
    :type data_list: list
    :param fill_value: missing data is filled with that value
    :type fill_value: float
    :param start_time_step: defines the time step at the start of the database (zero based); especially necessary, if an
        intermediate dataset with 15 min values is given to the function
    :type start_time_step: float
    :return: 1D list of forecast data content from grib files (all 15 minutes datasets in a row, e.g. in IconD2)
    :type: list
    """
    data_list_flattened = [quarter for hour in data_list for quarter in hour]  # flattened for each 15 min / 1 h
    num_data = len(data_list_flattened)

    for quarter in range(num_data - 1, 0, -1):
        # calculate instantaneous values from averaged
        # loop backwards through the list as the upper list element itself is changed in each step and must not
        # influence following steps
        if np.isnan(fill_value):
            idx_missing = np.isnan(data_list_flattened[quarter].data)
        else:
            idx_missing = data_list_flattened[quarter].data == fill_value
        data_list_flattened[quarter].data = (start_time_step + quarter) * data_list_flattened[
            quarter].data - (start_time_step + quarter - 1) * data_list_flattened[quarter - 1].data

        # correct negative values (due to machine precision) to zero
        idx_negative = data_list_flattened[quarter].data < 0
        if np.any(idx_negative):
            # print(f'negative values were corrected to 0')
            data_list_flattened[quarter].data[idx_negative] = 0
        data_list_flattened[quarter].data[idx_missing] = fill_value

    return data_list_flattened


def harmonize_time_step(data_list, fill_value, forecast_time):
    """
    Harmonize time step (dt) in data_list. The target dt is taken from the first two time steps in forecast_time. A
    changing dt is allowed at an arbitrary number of elements in data_list, as far they are at the end of that list and
    have an integer divider compared to the first dt. The data disaggregation uses a simple block design and fills all
    values evenly with the sum_disaggregated = sum_original/divider. If dt is consistent, the original data is returned.

    :param data_list: list of various data in met entities classes (e.g. IconEU)
    :type data_list: list
    :param fill_value: missing data is filled with that value
    :type fill_value: float
    :param forecast_time: list of values for forecast times
    :type forecast_time: list
    :return: new data_list and forecast_time objects
    :type: list, list
    """
    dt = np.array(np.array(forecast_time[1:]) - np.array(forecast_time[0:-1]))
    if all(dt == dt[0]):
        return data_list, forecast_time
    else:
        # time step harmonization necessary
        # find indexes with different dt
        idxs = np.where(dt != dt[0])[0]
        if np.max(idxs) != len(dt) - 1 and all(np.arange(min(idxs), max(idxs) + 1) == idxs):
            raise Exception('changing dt only allowed for tailing time steps')

        # change values
        data_disag = []

        for idx in idxs:
            divider = dt[idx] / dt[0]
            if not np.isclose(divider%1, 0):
                raise Exception(f'dt = {divider}, different dt must be integer values')
            data = copy.deepcopy(data_list[idx + 1])
            idx_missing = data.data == fill_value
            data.data = data.data / divider
            data.data[idx_missing] = fill_value
            for i in range(int(divider)):
                data_disag.append(data)
        for _ in idxs:
            data_list.pop()

        # replace original data in data_list
        for data in data_disag:
            data_list.append(data)

        # build updated forecast_time vector
        forecast_time_disag = list(np.linspace(forecast_time[0], forecast_time[-1], len(data_list)))

        return data_list, forecast_time_disag
