#   @author: Michael Wagner
#   @organization: TU Dresden
#   @contact: michael.wagner@tu-dresden.de
import os

import datetime as dt
import numpy as np
import netCDF4 as nc

from met_entities.LonLatTime import LonLatTime
from met_entities.MetEntities import MetEntities
from met_entities.GeoReferencedData import GeoReferencedData
from read_icon import readIcon as readIcon
from met_entities import VariableDescription as vd
from met_entities.Exceptions import *
from aux_tools import grib2_tools as gt
from aux_tools import scaling_tools as st


class IconD2(MetEntities):
    """
    IconD2 class provides all relevant data of IconD2 data from DWD.
    """

    def __init__(self,
                 time_value: LonLatTime = None,
                 forecast_value: LonLatTime = None,
                 gr_data: GeoReferencedData = None,
                 short: str = None):
        """
        Initialize IconD2 class.

        :param time_value: time
        :type time_value: LonLatTime.LonLatTime, optional
        :param forecast_value: forecast time
        :type forecast_value: LonLatTime.LonLatTime, optional
        :param gr_data: spatial data
        :type gr_data: GeoReferencedData.GeoReferencedData, optional
        :param short: if int16 or int32, the large data variables are cast to numpy.int16/int32 to minimize memory
            usage; the user must pay attention to the scale_factor to include the necessary precision
        :type short: str, optional
        """
        super().__init__(time_value=time_value, forecast_value=forecast_value, gr_data=gr_data, short=short)

        # added netcdf dimension and variable description, it can be changed from outside but not with the constructor
        self.nc_desc = vd.NcDimVarDescription(dim_time='time', dim_forecast='forecast_time', dim_lon='lon',
                                              dim_lat='lat', dim_coord='coord', var_time='time',
                                              var_forecast='forecast_time', var_lon='longitude', var_lat='latitude',
                                              var_data='icond2')

    def read_file(self, start_datetime, directory='./', dir_time_descriptor=None, forecast_hours=48, scale_factor=1,
                  fill_value=np.nan, short=None):
        """
        Reading in all available files of an IconD2 dataset. If a file for a specific forecast time is not available
        yet, the function automatically stops reading in and delivers a proper datastructure.

        :param start_datetime: time to start reading in IconD2 files, it is used to build the filenames following the
            convention of the german weather service (DWD)
        :type start_datetime: datetime.datetime
        :param directory: directory with all IconD2 files from start_datetime
        :type directory: str, optional
        :param dir_time_descriptor: list of datetime.strftime time descriptors for an arbitrary number of additional
            time dependent folders with the data, e.g. ['%Y', '%Y%m%d'] for */yyyy/yyyymmdd/*; the final path is built
            to directory/dir_time_directory[0]/../dir_time_directory[n]/filename
        :type dir_time_descriptor: list, optional
        :param forecast_hours: the number of forecast hours to be read in, default is 48 (max for IconD2)
        :type forecast_hours: int, optional
        :param scale_factor: the final data has to be multiplied with this value
        :type scale_factor: float, optional
        :param fill_value: missing data is filled with that value
        :type fill_value: float, optional
        :param short: if int16 or int32, the large data variables are cast to numpy.int16/int32 to minimize memory
            usage; the user must pay attention to the scale_factor to include the necessary precision
        :type short: str, optional
        """
        if short:
            self.short = short

        # get coordinates
        clon, clat = readIcon.get_lonlat(variant='d2')

        # build filenames and read in data
        filenames = []
        icon_data = []
        forecast_value = []
        for i in range(forecast_hours + 1):
            filename = 'icon-d2_germany_icosahedral_single-level_' + start_datetime.strftime('%Y%m%d%H') + '_' + \
                       '{:03d}'.format(i) + '_2d_tot_prec.grib2.bz2'
            print(f'load {filename}')
            if dir_time_descriptor is not None:
                dir_time = ''
                for dir_act in dir_time_descriptor:
                    dir_time = os.path.join(dir_time, start_datetime.strftime(dir_act))
            else:
                dir_time = ''
            filenames.append(os.path.join(directory, dir_time, filename))
            try:
                icon_data.append(readIcon.read_icon_d2(filename=filenames[i], scale_factor=scale_factor,
                                                       fill_value=fill_value, variant='d2', short=self.short))
            except IOError:
                if i > 0:
                    print(f'Warning: file {filename} not available - stopping at {filenames[i - 1]}')
                else:
                    raise ForecastFileNotAvailable(f'first forecast file {filename} not available - import not successful')
                break

        icon_data_flattened = gt.accum_to_instantaneous_flattened(icon_data, fill_value=fill_value)
        for icon_data_act in icon_data_flattened:
            forecast_value.append(icon_data_act.metadata.prediction_time.total_seconds() / 60)
        forecast_description = vd.TimeDescription(units='minutes')

        self.forecast_value = LonLatTime(data=np.array(forecast_value), data_description=forecast_description)
        self.gr_data = GeoReferencedData()
        self.gr_data.lon = clon
        self.gr_data.lat = clat
        if self.short == 'int16':
            self.gr_data.data = np.empty((icon_data_flattened[0].data.size, len(icon_data_flattened)), dtype=np.short)
        elif self.short == 'int32':
            self.gr_data.data = np.empty((icon_data_flattened[0].data.size, len(icon_data_flattened)), dtype=np.int32)
        else:
            self.gr_data.data = np.empty((icon_data_flattened[0].data.size, len(icon_data_flattened)))
        for i in range(len(icon_data_flattened)):
            # convert list entries to ndarray
            self.gr_data.data[:, i] = icon_data_flattened[i].data
        self.gr_data.data_description = vd.DataDescription(fill_value=fill_value, scale_factor=scale_factor,
                                                           long_name='IconD2 precipitation data',
                                                           units=f'{1 / scale_factor} * '
                                                                 f'{icon_data_flattened[0].metadata.units}',
                                                           time_note=f'start at {icon_data_flattened[0].metadata.datum_iso}')
        time_value = 0
        time_unit = f"hours since {icon_data_flattened[0].metadata.datum_iso}"
        time_description = vd.TimeDescription(calendar="standard", units=time_unit)
        self.time_value = LonLatTime(data=time_value, data_description=time_description)

    def regrid(self, regrid_description=None, lon_target=None, lat_target=None, neighbors=3, file_nearest=None):
        """
        Interpolate the gridded data onto a new raster in the same coordinate system. It uses an inverse distance
        weighted (IDW) method with an arbitrary number of neighbors. If regrid_description is given it is prioritized.

        :param regrid_description: regrid description with some of lon_target, lat_target, neighbors, file_nearest
            variables
        :type regrid_description: VariableDescription.RegridDescription, optional
        :param lon_target: longitudes of the target raster, given as raster with shape (num_lat, num_lon)
        :type lon_target: LonLatTime.LonLatTime, optional
        :param lat_target: latitudes of the target raster, given as raster with shape (num_lat, num_lon)
        :type lat_target: LonLatTime.LonLatTime, optional
        :param neighbors: number of neighbors for IDW
        :type neighbors: int, optional
        :param file_nearest: npz file with indexes and lengths for the current setup
        :type file_nearest: str, optional
        """
        print('regridding')
        if regrid_description is not None:
            self.gr_data.regrid_idw(lon_target=regrid_description.lon_target, lat_target=regrid_description.lat_target,
                                    neighbors=regrid_description.neighbors,
                                    file_nearest=regrid_description.file_nearest, short=self.short)
        else:
            self.gr_data.regrid_idw(lon_target=lon_target, lat_target=lat_target, neighbors=neighbors,
                                    file_nearest=file_nearest, short=self.short)

    def crop(self, crop_description=None, lon_west=None, lon_east=None, lat_south=None, lat_north=None, idx_west=None,
             idx_east=None, idx_south=None, idx_north=None, idx_array=None):
        """
        Cropping of data of this IconD2 class. Usage of indexes directly if given. Otherwise, use lon/lat with the
        guarantee that the whole requested area is within the cropped region. If crop_description is given it is
        prioritized.

        :param crop_description: crop description with some of lon_west, lon_east, lat_south, lat_north, idx_west,
            idx_east, idx_south, idx_north, idx_array variables
        :type crop_description: VariableDescription.CropDescription, optional
        :param lon_west: western longitude limit
        :type lon_west: float, optional
        :param lon_east: eastern longitude limit
        :type lon_east: float, optional
        :param lat_south: southern latitude limit
        :type lat_south: float, optional
        :param lat_north: northern latitude limit
        :type lat_north: float, optional
        :param idx_west: western index limit
        :type idx_west: int, optional
        :param idx_east: eastern index limit
        :type idx_east: int, optional
        :param idx_south: southern index limit
        :type idx_south: int, optional
        :param idx_north: northern index limit
        :type idx_north: int, optional
        :param idx_array: index for 1D array in lon and lat (e.g. original icond2 data)
        :type idx_array: np.ndarray, optional
        """
        print('cropping')
        if crop_description is not None:
            self.gr_data.crop(lon_west=crop_description.lon_west, lon_east=crop_description.lon_east,
                              lat_south=crop_description.lat_south, lat_north=crop_description.lat_north,
                              idx_west=crop_description.idx_west, idx_east=crop_description.idx_east,
                              idx_south=crop_description.idx_south, idx_north=crop_description.idx_north,
                              idx_array=crop_description.idx_array)
        else:
            self.gr_data.crop(lon_west=lon_west, lon_east=lon_east, lat_south=lat_south, lat_north=lat_north,
                              idx_west=idx_west, idx_east=idx_east, idx_south=idx_south, idx_north=idx_north,
                              idx_array=idx_array)

    def export_netcdf(self, filename, data_format='f8', institution=None, scale_factor_nc=1, scale_undo=False,
                      data_kwargs=None):
        """
        Export the relevant content of this IconD2 class to a new netcdf file.

        :param filename: filename of netcdf file
        :type filename: str
        :param data_format: format description of resulting data field in netcdf, the following specifiers are allowed:
            'f8', 'f4', 'i8', 'i4', 'i2', 'i1', 'u8', 'u4', 'u2', 'u1', 'S1' (f - float, i - integer, u - unsigned
            integer, S1 - single character string; the number specifies the number of bytes)
        :type data_format: str, optional
        :param institution: description of the institution generating the netcdf file; will be used for global netcdf
            attribute 'institution'
        :type institution: str
        :param scale_factor_nc: is used as scale_factor for netcdf file and can be used for storage saving purposes in
            conjunction with data_format (e.g. 'i2'); do not mix it up with scale_factor_internal
        :type scale_factor_nc: float, optional
        :param scale_undo: if True, the original internal scale_factor is taken back in order to get the original
            values.
        :type scale_undo: bool, optional
        :param data_kwargs: keyword arguments that are passed to netCDF4.createVariable for data variables, for
            supported arguments refer to https://unidata.github.io/netcdf4-python/#Dataset.createVariable;
            e.g. {'compression': 'zlib'} compresses the data with zlib alorithm and default complevel=4
        :type data_kwargs: dict, optional
        """
        print('exporting to netcdf')
        if data_kwargs is None:
            data_kwargs = {}
        if self.gr_data.lon.data.ndim == 2:
            # coordinates given as matrix
            coord_matrix = True
        else:
            # coordinates given as array
            coord_matrix = False

        if scale_undo and self.gr_data.data_description.scale_factor == 1:
            # roll back scaling called, but not necessary
            scale_undo = False

        # create file and dimensions
        ncfile = nc.Dataset(filename, mode='w', format='NETCDF4')
        if coord_matrix:
            ncfile.createDimension(self.nc_desc.dim_lon, self.gr_data.lon.data.shape[1])  # longitude axis
            ncfile.createDimension(self.nc_desc.dim_lat, self.gr_data.lat.data.shape[0])  # latitude axis
        else:
            ncfile.createDimension(self.nc_desc.dim_coord, self.gr_data.lon.data.size)  # coordinate axis
        ncfile.createDimension(self.nc_desc.dim_forecast, self.forecast_value.data.size)  # forecast axis
        ncfile.createDimension(self.nc_desc.dim_time, None)  # unlimited time axis
        if not self.gr_data.regridded and not self.gr_data.cropped:
            ncfile.title = 'IconD2 data'
        elif self.gr_data.regridded and not self.gr_data.cropped:
            ncfile.title = 'IconD2 data (regridded)'
        elif not self.gr_data.regridded and self.gr_data.cropped:
            ncfile.title = 'IconD2 data (cropped)'
        else:
            ncfile.title = 'IconD2 data (regridded and cropped)'
        ncfile.history = 'v0'
        if institution is not None:
            ncfile.institution = institution
        ncfile.generator = 'generated with weatherDataHarmonizer'
        ncfile.source = 'IconD2 data from DWD'
        ncfile.created = dt.datetime.now(tz=dt.timezone.utc).isoformat()

        # create variables in file
        if coord_matrix:
            lon = ncfile.createVariable(self.nc_desc.var_lon, np.float32, (self.nc_desc.dim_lat, self.nc_desc.dim_lon),
                                        fill_value=self.gr_data.lon.data_description.fill_value)
            lat = ncfile.createVariable(self.nc_desc.var_lat, np.float32, (self.nc_desc.dim_lat, self.nc_desc.dim_lon),
                                        fill_value=self.gr_data.lat.data_description.fill_value)
            icond2 = ncfile.createVariable(self.nc_desc.var_data, data_format,
                                           (self.nc_desc.dim_time, self.nc_desc.dim_lat, self.nc_desc.dim_lon,
                                            self.nc_desc.dim_forecast),
                                           fill_value=self.gr_data.data_description.fill_value, **data_kwargs)
        else:
            lon = ncfile.createVariable(self.nc_desc.var_lon, np.float32, (self.nc_desc.dim_coord,),
                                        fill_value=self.gr_data.lon.data_description.fill_value)
            lat = ncfile.createVariable(self.nc_desc.var_lat, np.float32, (self.nc_desc.dim_coord,),
                                        fill_value=self.gr_data.lon.data_description.fill_value)
            icond2 = ncfile.createVariable(self.nc_desc.var_data, data_format,
                                           (self.nc_desc.dim_time, self.nc_desc.dim_coord, self.nc_desc.dim_forecast),
                                           fill_value=self.gr_data.data_description.fill_value, **data_kwargs)
        lon.units = self.gr_data.lon.data_description.units
        lon.long_name = self.gr_data.lon.data_description.long_name
        lon.coordinate_system = self.gr_data.lon.data_description.coordinate_system

        lat.units = self.gr_data.lat.data_description.units
        lat.long_name = self.gr_data.lat.data_description.long_name
        lat.coordinate_system = self.gr_data.lat.data_description.coordinate_system

        time = ncfile.createVariable(self.nc_desc.var_time, np.float64, (self.nc_desc.dim_time,))
        time.units = self.time_value.data_description.units
        time.calendar = self.time_value.data_description.calendar

        forecast_time = ncfile.createVariable(self.nc_desc.var_forecast, 'i2', (self.nc_desc.dim_forecast,))
        forecast_time.units = self.forecast_value.data_description.units

        precipitation_units = self.gr_data.data_description.units.split(' ')
        if scale_undo:
            precip_multiplier = float(precipitation_units[0]) * self.gr_data.data_description.scale_factor
            scale_factor_internal = 1
        else:
            precip_multiplier = precipitation_units[0]
            scale_factor_internal = self.gr_data.data_description.scale_factor
        precipitation_units = f'{precip_multiplier} {" ".join(precipitation_units[1:])}'

        icond2.units = precipitation_units
        icond2.long_name = self.gr_data.data_description.long_name
        icond2.scale_factor_internal = scale_factor_internal
        icond2.scale_factor = scale_factor_nc

        # write data to variables
        if coord_matrix:
            lon[:, :] = self.gr_data.lon.data
            lat[:, :] = self.gr_data.lat.data
            icond2[0, :, :, :] = st.gr_data_scaling(self.gr_data, scale_undo,
                                                    self.gr_data.data_description.scale_factor,
                                                    scale_factor_nc)
        else:
            lon[:] = self.gr_data.lon.data
            lat[:] = self.gr_data.lat.data
            icond2[0, :, :] = st.gr_data_scaling(self.gr_data, scale_undo, self.gr_data.data_description.scale_factor,
                                                 scale_factor_nc)
        time[:] = self.time_value.data
        forecast_time[:] = self.forecast_value.data

        ncfile.close()

    def export_netcdf_append(self, filename):
        """
        Append the relevant content of this IconD2 class to an existing netcdf file.

        :param filename: filename of netcdf file
        :type filename: str
        """
        print('append to existing netcdf')
        # append data to existing file
        if not os.path.exists(filename):
            raise Exception(f'file {filename} does not exist yet')
        ncfile = nc.Dataset(filename, mode='a')

        # check dimensions, shape
        if ((self.nc_desc.dim_lon not in ncfile.dimensions or self.nc_desc.dim_lat not in ncfile.dimensions) and
                self.nc_desc.dim_coord not in ncfile.dimensions) or \
                self.nc_desc.dim_time not in ncfile.dimensions or \
                self.nc_desc.dim_forecast not in ncfile.dimensions or \
                self.nc_desc.var_time not in ncfile.variables or \
                self.nc_desc.var_data not in ncfile.variables or \
                self.nc_desc.var_forecast not in ncfile.variables:
            raise Exception(f'either dimensions or variables are not appropriate in {filename}')

        if self.gr_data.lon.data.ndim == 2:
            # coordinates given as matrix
            coord_matrix = True
        else:
            # coordinates given as array
            coord_matrix = False

        if coord_matrix:
            if ncfile[self.nc_desc.var_lon].shape != self.gr_data.lon.data.shape:
                raise Exception(f'the data shape in {filename} differs from object shape')
        else:
            if ncfile[self.nc_desc.var_lon].size != self.gr_data.lon.data.size:
                raise Exception(f'the data shape in {filename} differs from object shape')
        if ncfile[self.nc_desc.var_forecast].size != self.forecast_value.data.size:
            raise Exception(f'the forecast shape in {filename} differs from object shape')

        # check scale_factor and fill_value
        scale_undo = False
        if self.gr_data.data_description.scale_factor != 1 and ncfile[self.nc_desc.var_data].scale_factor_internal == 1:
            scale_undo = True
        if ncfile[self.nc_desc.var_data].scale_factor_internal != 1 and \
                self.gr_data.data_description.scale_factor != ncfile[self.nc_desc.var_data].scale_factor_internal:
            raise Exception('internal scale factor in netcdf file differs from actual internal scale factor')
        nc_fill_value = ncfile[self.nc_desc.var_data]._FillValue
        if np.isnan(self.gr_data.data_description.fill_value) and not np.isnan(nc_fill_value):
            raise Exception('fill value in netcdf file differs from actual fill value')
        elif not np.isnan(self.gr_data.data_description.fill_value) and \
                self.gr_data.data_description.fill_value != nc_fill_value:
            raise Exception('fill value in netcdf file differs from actual fill value')

        # check time and evaluate correct time value for netcdf content
        time_units_nc = ncfile[self.nc_desc.var_time].units.split(' ')
        time_units_internal = self.time_value.data_description.units.split(' ')
        if time_units_nc[0] != 'hours':
            # currently only hours supported as time step description
            raise Exception(f'time step {time_units_nc[0]} in {filename} not supported')
        if time_units_internal[0] != 'hours':
            raise Exception(f'internal time step {time_units_internal[0]} not supported')

        time_start_nc = dt.datetime.fromisoformat(time_units_nc[-1])
        time_start_internal = dt.datetime.fromisoformat(time_units_internal[-1])
        time_act = time_start_internal + dt.timedelta(hours=self.time_value.data) - time_start_nc

        # write data
        ncfile.last_modified = dt.datetime.now(tz=dt.timezone.utc).isoformat()
        num_time_nc = ncfile[self.nc_desc.var_time].size  # actual number of time steps
        ncfile[self.nc_desc.var_time][num_time_nc] = time_act.total_seconds() / 3600
        scale_factor_nc = ncfile[self.nc_desc.var_data].scale_factor
        if coord_matrix:
            ncfile[self.nc_desc.var_data][num_time_nc, :, :, :] = \
                st.gr_data_scaling(self.gr_data, scale_undo, self.gr_data.data_description.scale_factor,
                                   scale_factor_nc)
        else:
            ncfile[self.nc_desc.var_data][num_time_nc, :, :] = \
                st.gr_data_scaling(self.gr_data, scale_undo, self.gr_data.data_description.scale_factor,
                                   scale_factor_nc)

        ncfile.close()
