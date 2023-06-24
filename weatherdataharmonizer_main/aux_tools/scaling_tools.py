#   @author: Michael Wagner
#   @organization: TU Dresden
#   @contact: michael.wagner@tu-dresden.de

import numpy as np


def get_fill_value_unpack(gr_data, scale_factor_nc):
    """
    Include the scale factor for missing values only in the geodata. This is used for packed netcdf export.

    :param gr_data: one set of geodata
    :type gr_data: met_entities.GeoReferencedData.GeoReferencedData
    :param scale_factor_nc: scale_factor for netcdf writing; leads to the variable attribute scale_factor
    :type scale_factor_nc: float
    :return: nd-matrix with data from gr_data and missing values multiplied with scale_factor_nc
    :rtype: numpy.ndarray
    """

    fill_value_unpack = np.ones(gr_data.data.shape)
    fill_value_unpack[gr_data.data == gr_data.data_description.fill_value] = scale_factor_nc

    return fill_value_unpack

def gr_data_scaling(gr_data, scale_undo, scale_factor, scale_factor_nc):
    """
    Prepare geodata for export to netcdf. The internal scaling can be rolled back. If the netcdf variable shall be
    packed the missing values have to be multiplied with the according scale factor first. This is necessary as the
    packing ist done before writing to netcdf and the final values in netcdf raw data are compared to the attribute
    _FillValue to mark missing values.

    :param gr_data: one set of geodata
    :type gr_data: met_entities.GeoReferencedData.GeoReferencedData
    :param scale_undo: if True the internal scale factor is taken back
    :type scale_undo: bool
    :param scale_factor: internal scale factor; e.g. used at reading in large data to save memory; leads to the variable
        attribute scale_factor_internal
    :type scale_factor: float
    :param scale_factor_nc: scale_factor for netcdf writing; leads to the variable attribute scale_factor
    :type scale_factor_nc: float
    :return: nd-matrix with back scaled data from gr_data and missing values multiplied with scale_factor_nc
    :rtype: numpy.ndarray
    """

    fill_value_unpack = get_fill_value_unpack(gr_data, scale_factor_nc)

    if scale_undo:
        value_scale_undo = np.ones(gr_data.data.shape)
        value_scale_undo[gr_data.data != gr_data.data_description.fill_value] = scale_factor
        return gr_data.data * value_scale_undo * fill_value_unpack
    else:
        return gr_data.data * fill_value_unpack
