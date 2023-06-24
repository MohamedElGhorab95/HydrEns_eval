#   @author: Michael Wagner
#   @organization: TU Dresden
#   @contact: michael.wagner@tu-dresden.de

import netCDF4 as nc
import numpy as np

from met_entities.LonLatTime import LonLatTime
import met_entities.VariableDescription as vd
from met_entities.GeoReferencedData import GeoReferencedData
from met_entities.RadolanRW import RadolanRW
from met_entities.RadvorRQ import RadvorRQ
from met_entities.RadolanRV import RadolanRV
from met_entities.IconD2 import IconD2
from met_entities.IconD2EPS import IconD2EPS
from met_entities.IconEU import IconEU
from met_entities.IconEUEPS import IconEUEPS
from met_entities.CosmoD2 import CosmoD2
from met_entities.CosmoD2EPS import CosmoD2EPS


def read_netcdf(filename, scale_factor_result=None, fill_value_result=None, short_result=None):
    """
    Read in data from netcdf file formerly created with export_netcdf in instances of MetEntities subclasses. The type
    is derived automatically from variable naming in the netcdf file. It results in a MetEntities subclass (e.g.
    RadolanRW, IconD2EPS or others) with all variables and descriptions as in the class content if loaded from raw data.
    The scale_factor in netcdf is multiplied with scale_factor_internal to ensure precision preservation.
    Explained: The scale_factor is typically only used for packing the data in netcdf to save hard disk storage.
    Additionally, the datatype from the netcdf file is also taken as internal python datatype. Without multiplication of
    scale_factor the data could lose precision in case of integer datatypes (which are used to save memory).

    Additionally, a scale_factor_result can be given to request a specific scale_factor that shall be valid in the end.

    :param filename: filename of the netcdf file that has the format used in MetEntities.export_netcdf.
    :type filename: str
    :param scale_factor_result: the final data has to be multiplied with this value; if None the applied scale_factor
        is the multiplication of scale_factor in netcdf and scale_factor_internal (see description above)
    :type scale_factor_result: float, optional
    :param fill_value_result: fill value of the final result; if given, the netcdf internally used fill value is
        replaced with this one
    :type fill_value_result: float, optional
    :param short_result: if int16, int32 or None, the resulting netcdf variables are cast to numpy.int16/int32/float no
        matter the data type in netcdf; the user must pay attention to the scale_factor_result to ensure the necessary
        precision
    :type short_result: str, optional
    :return: instance of a MetEntities subclass
    :rtype: met_entities.MetEntities.MetEntities
    """
    print(f'load {filename}')

    if short_result and scale_factor_result >= 1:
        print(f'the attribute short casts the data to 16/32 bit signed integer for memory saving; the scale_factor = '
              f'{scale_factor_result} imposes a low precision')
    if short_result and np.isnan(fill_value_result):
        raise Exception(f'if short data type is requested, the fill_value must not be nan')
    if short_result and (short_result != 'int16' and short_result != 'int32'):
        raise Exception(f'datatype {short_result} not supported; please change to int16, int32, or None')

    data = None

    # open netcdf
    with nc.Dataset(filename, mode='r') as ncfile:
        nc_vars = ncfile.variables

        # define target class by netcdf variable name
        switch_forecast = True
        switch_eps = False
        for key in nc_vars.keys():
            if key.lower().startswith('radolanrw'):
                data = RadolanRW()
                switch_forecast = False
                break
            elif key.lower().startswith('radvorrq'):
                data = RadvorRQ()
                break
            elif key.lower().startswith('radolanrv'):
                data = RadolanRV()
                break
            elif key.lower().startswith('icond2eps'):
                data = IconD2EPS()
                switch_eps = True
                break
            elif key.lower().startswith('icond2') and not key.lower().startswith('icond2eps'):
                data = IconD2()
                break
            elif key.lower().startswith('iconeueps'):
                data = IconEUEPS()
                switch_eps = True
                break
            elif key.lower().startswith('iconeu') and not key.lower().startswith('iconeueps'):
                data = IconEU()
                break
            elif key.lower().startswith('cosmod2eps'):
                data = CosmoD2EPS()
                switch_eps = True
                break
            elif key.lower().startswith('cosmod2') and not key.lower().startswith('cosmod2eps'):
                data = CosmoD2()
                break

        # check variable availability and dimensions
        if data.nc_desc.var_lon not in nc_vars:
            raise Exception(f'variable {data.nc_desc.var_lon} not in file {filename}')
        if data.nc_desc.var_lat not in nc_vars:
            raise Exception(f'variable {data.nc_desc.var_lat} not in file {filename}')
        if data.nc_desc.var_time not in nc_vars:
            raise Exception(f'variable {data.nc_desc.var_time} not in file {filename}')
        if switch_forecast and data.nc_desc.var_forecast not in nc_vars:
            raise Exception(f'variable {data.nc_desc.var_forecast} not in file {filename}')
        if switch_eps and data.nc_desc.var_data in nc_vars and data.nc_desc.var_eps not in nc_vars:
            raise Exception(f'variable {data.nc_desc.var_eps} not in file {filename}')

        nc_dims = ncfile.dimensions
        for nc_dim in nc_dims:
            if nc_dim not in [data.nc_desc.dim_time, data.nc_desc.dim_forecast, data.nc_desc.dim_lon,
                              data.nc_desc.dim_lat, data.nc_desc.dim_coord, data.nc_desc.dim_eps]:
                raise Exception(f'in file {filename} the dimension {nc_dim.name} cannot be interpreted')

        variable_names = []
        if switch_eps:
            data.eps_member = []
            if data.nc_desc.var_data in nc_vars:
                data.eps_member = nc_vars[data.nc_desc.var_eps][:]
                variable_names.append(data.nc_desc.var_data)
            else:
                for key in nc_vars.keys():
                    if key.startswith(data.nc_desc.var_data):
                        data.eps_member.append(int(key.split('_')[-1]))
                        variable_names.append(key)
        else:
            variable_names.append(data.nc_desc.var_data)

        switch_coords_array = False
        if data.nc_desc.dim_coord in nc_dims:
            # coordinates dimension available -> vector of coordinates
            switch_coords_array = True

        # fill_value handling
        if fill_value_result is None:
            fill_value_result = nc_vars[variable_names[0]]._FillValue

        # data type handling
        if short_result is None:
            short_result = 'f8'
        elif short_result == 'int32':
            short_result = 'i4'
            data.short = 'int32'
        elif short_result == 'int16':
            short_result = 'i2'
            data.short = 'int16'

        # fill data variables
        time_description = vd.TimeDescription(calendar=nc_vars[data.nc_desc.var_time].calendar,
                                              units=nc_vars[data.nc_desc.var_time].units)
        data.time_value = LonLatTime(data=nc_vars[data.nc_desc.var_time][:].data, data_description=time_description)
        if switch_forecast:
            forecast_description = vd.TimeDescription(units=nc_vars[data.nc_desc.var_forecast].units)
            data.forecast_value = LonLatTime(data=nc_vars[data.nc_desc.var_forecast][:].data,
                                             data_description=forecast_description)

        lon_description = vd.DataDescription(units=nc_vars[data.nc_desc.var_lon].units,
                                             long_name=nc_vars[data.nc_desc.var_lon].long_name,
                                             coordinate_system=nc_vars[data.nc_desc.var_lon].coordinate_system)
        lat_description = vd.DataDescription(units=nc_vars[data.nc_desc.var_lat].units,
                                             long_name=nc_vars[data.nc_desc.var_lat].long_name,
                                             coordinate_system=nc_vars[data.nc_desc.var_lat].coordinate_system)
        precipitation_units = nc_vars[variable_names[0]].units.split(' ')

        scale_factor_nc = nc_vars[variable_names[0]].scale_factor
        scale_factor_internal = nc_vars[variable_names[0]].scale_factor_internal
        if scale_factor_result is not None and scale_factor_nc * scale_factor_internal != scale_factor_result:
            # multiplication of internal and netcdf scale_factor to preserve precision
            scale_factor_result_multiplicator = scale_factor_nc / scale_factor_result
        else:
            scale_factor_result = scale_factor_nc * scale_factor_internal
            scale_factor_result_multiplicator = 1

        precipitation_units = f'{float(precipitation_units[0]) / scale_factor_nc} {" ".join(precipitation_units[1:])}'
        data_description = vd.DataDescription(fill_value=fill_value_result,
                                              scale_factor=scale_factor_result,
                                              long_name=nc_vars[variable_names[0]].long_name,
                                              time_note=f'start at {nc_vars[data.nc_desc.var_time].units.split(" ")[-1]}',
                                              units=precipitation_units)

        if not hasattr(data, 'eps_member') or data.eps_member is None or \
                (hasattr(data, 'eps_member') and len(data.eps_member) == 1):
            # no eps member or only one eps_member
            data.gr_data = GeoReferencedData()
            data.gr_data.lon = LonLatTime(data=nc_vars[data.nc_desc.var_lon][:].data, data_description=lon_description)
            data.gr_data.lat = LonLatTime(data=nc_vars[data.nc_desc.var_lat][:].data, data_description=lat_description)
            data.gr_data.data_description = data_description
            if switch_coords_array:
                dims = nc_vars[variable_names[0]].dimensions
                if switch_forecast:
                    dim_transpose_map = [dims.index(data.nc_desc.dim_coord), dims.index(data.nc_desc.dim_forecast),
                                         dims.index(data.nc_desc.dim_time)]
                else:
                    dim_transpose_map = [dims.index(data.nc_desc.dim_coord), dims.index(data.nc_desc.dim_time)]
            else:
                dims = nc_vars[variable_names[0]].dimensions
                if switch_forecast:
                    dim_transpose_map = [dims.index(data.nc_desc.dim_lat), dims.index(data.nc_desc.dim_lon),
                                         dims.index(data.nc_desc.dim_forecast), dims.index(data.nc_desc.dim_time)]
                else:
                    dim_transpose_map = [dims.index(data.nc_desc.dim_lat), dims.index(data.nc_desc.dim_lon),
                                         dims.index(data.nc_desc.dim_time)]
            data.gr_data.data = (
                    np.squeeze(np.transpose(nc_vars[variable_names[0]][:].data, axes=dim_transpose_map)) /
                    scale_factor_nc * scale_factor_result_multiplicator).astype(short_result)
            data.gr_data.data[data.gr_data.data == nc_vars[variable_names[0]]._FillValue /
                              scale_factor_nc * scale_factor_result_multiplicator] = fill_value_result

            # remove forecast time steps that only contain fill values (only for datasets with one single time step)
            if switch_forecast and nc_vars[variable_names[0]].shape[dims.index(data.nc_desc.dim_time)] == 1:
                num_forecasts = data.forecast_value.data.size
                num_coord_cells = 1
                for size_act in data.gr_data.lon.data.shape:
                    num_coord_cells = num_coord_cells * size_act
                idx_empty_forecast = np.all(
                    data.gr_data.data.reshape(num_coord_cells, num_forecasts) == fill_value_result, axis=0)
                if any(idx_empty_forecast):
                    data.forecast_value.data = data.forecast_value.data[~idx_empty_forecast]
                    if switch_coords_array:
                        data.gr_data.data = data.gr_data.data[:, ~idx_empty_forecast]
                    else:
                        data.gr_data.data = data.gr_data.data[:, :, ~idx_empty_forecast]

        else:
            # multiple eps members
            data.gr_data = []
            idx_empty_forecast = False
            for ct in range(len(data.eps_member)):
                gr_data = GeoReferencedData()
                gr_data.lon = LonLatTime(data=nc_vars[data.nc_desc.var_lon][:].data, data_description=lon_description)
                gr_data.lat = LonLatTime(data=nc_vars[data.nc_desc.var_lat][:].data, data_description=lat_description)
                data_description.fill_value = fill_value_result
                data_description.eps_note = f'eps member {data.eps_member[ct]}'
                gr_data.data_description = data_description

                if switch_coords_array:
                    dims = nc_vars[variable_names[ct]].dimensions
                    dim_transpose_map = [dims.index(data.nc_desc.dim_coord), dims.index(data.nc_desc.dim_forecast),
                                         dims.index(data.nc_desc.dim_time)]
                else:
                    dims = nc_vars[variable_names[ct]].dimensions
                    dim_transpose_map = [dims.index(data.nc_desc.dim_lat), dims.index(data.nc_desc.dim_lon),
                                         dims.index(data.nc_desc.dim_forecast), dims.index(data.nc_desc.dim_time)]
                gr_data.data = (
                        np.squeeze(np.transpose(nc_vars[variable_names[ct]][:].data, axes=dim_transpose_map)) /
                        scale_factor_nc * scale_factor_result_multiplicator).astype(short_result)
                gr_data.data[gr_data.data == nc_vars[variable_names[ct]]._FillValue /
                             scale_factor_nc * scale_factor_result_multiplicator] = fill_value_result

                # remove forecast time steps that only contain fill values (only for datasets with one single time step)
                if ct == 0 and nc_vars[variable_names[0]].shape[dims.index(data.nc_desc.dim_time)] == 1:
                    num_forecasts = data.forecast_value.data.size
                    num_coord_cells = 1
                    for size_act in gr_data.lon.data.shape:
                        num_coord_cells = num_coord_cells * size_act
                    idx_empty_forecast = np.all(
                        gr_data.data.reshape(num_coord_cells, num_forecasts) == fill_value_result, axis=0)
                if any(idx_empty_forecast):
                    data.forecast_value.data = data.forecast_value.data[~idx_empty_forecast]
                    if switch_coords_array:
                        gr_data.data = gr_data.data[:, ~idx_empty_forecast]
                    else:
                        gr_data.data = gr_data.data[:, :, ~idx_empty_forecast]

                data.gr_data.append(gr_data)

    return data
