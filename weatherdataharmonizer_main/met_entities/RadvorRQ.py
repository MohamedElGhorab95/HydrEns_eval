#   @author: Michael Wagner
#   @organization: TU Dresden
#   @contact: michael.wagner@tu-dresden.de
import datetime as dt
import os

import numpy as np
import netCDF4 as nc

from met_entities.GeoReferencedData import GeoReferencedData
from met_entities.MetEntities import MetEntities
from read_radolan import readRadolan as readRadolan
from met_entities import VariableDescription as vd
from met_entities.LonLatTime import LonLatTime
from met_entities.Exceptions import *
from aux_tools import scaling_tools as st


class RadvorRQ(MetEntities):
    """
    RadvorRQ class provides all relevant data of RadvorRQ data from DWD.
    """

    def __init__(self,
                 time_value: LonLatTime = None,
                 forecast_value: LonLatTime = None,
                 gr_data: GeoReferencedData = None,
                 short: str = None):
        """
        Initialize RadvorRQ class.

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
                                              dim_lat='lat', var_time='time', var_forecast='forecast_time',
                                              var_lon='longitude', var_lat='latitude', var_data='radvorrq')

    def read_file(self, start_datetime, directory='./', dir_time_descriptor=None, scale_factor=1, fill_value=np.nan,
                  short=None):
        """
        Reading in three files of RadvorRQ data: at 0, +60, +120 minutes.

        :param start_datetime: time to start reading in RadvorRQ files, it is used to build the filenames following the
            convention of the german weather service (DWD)
        :type start_datetime: datetime.datetime
        :param directory: directory with all RadvorRQ files from start_datetime
        :type directory: str, optional
        :param dir_time_descriptor: list of datetime.strftime time descriptors for an arbitrary number of additional
            time dependent folders with the data, e.g. ['%Y', '%Y%m%d'] for */yyyy/yyyymmdd/*; the final path is built
            to directory/dir_time_directory[0]/../dir_time_directory[n]/filename
        :type dir_time_descriptor: list, optional
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

        self.gr_data = GeoReferencedData()

        # get metadata exemplary for the first file to obtain forecast times
        if dir_time_descriptor is not None:
            dir_time = ''
            for dir_act in dir_time_descriptor:
                dir_time = os.path.join(dir_time, start_datetime.strftime(dir_act))
        else:
            dir_time = ''
        file_prefix = os.path.join(directory, dir_time, 'RQ' + start_datetime.strftime('%y%m%d%H%M'))
        print(f'load {file_prefix}_+++.gz')
        try:
            metadata_data_000 = readRadolan.read_radolan(f"{file_prefix}_000.gz", metadata_only=True)
            metadata_data_060 = readRadolan.read_radolan(f"{file_prefix}_060.gz", metadata_only=True)
            metadata_data_120 = readRadolan.read_radolan(f"{file_prefix}_120.gz", metadata_only=True)
        except FileNotFoundError:
            raise ForecastFileNotAvailable(f'forecast for {start_datetime} not available - import not successful')

        forecast_value = np.array(
            [metadata_data_000.prediction_time, metadata_data_060.prediction_time, metadata_data_120.prediction_time])
        forecast_description = vd.TimeDescription(units='minutes')
        self.forecast_value = LonLatTime(data=forecast_value, data_description=forecast_description)

        self.gr_data.lon, self.gr_data.lat = readRadolan.get_lonlat(metadata_data_000.format_version,
                                                                    grid_variant='radolanrx')

        # create a 3D matrix with lat x lon x num_forecast RadvorRQ values
        self.gr_data.data = np.empty((self.gr_data.lon.data.shape[0], self.gr_data.lon.data.shape[1], 3))
        self.gr_data.data[:, :, 0] = readRadolan.read_radolan(f"{file_prefix}_000.gz", scale_factor=scale_factor,
                                                              fill_value=fill_value, short=self.short).data
        self.gr_data.data[:, :, 1] = readRadolan.read_radolan(f"{file_prefix}_060.gz", scale_factor=scale_factor,
                                                              fill_value=fill_value, short=self.short).data
        self.gr_data.data[:, :, 2] = readRadolan.read_radolan(f"{file_prefix}_120.gz", scale_factor=scale_factor,
                                                              fill_value=fill_value, short=self.short).data

        self.gr_data.data_description = vd.DataDescription(fill_value=fill_value,
                                                           scale_factor=scale_factor,
                                                           long_name='RadvorRQ precipitation data',
                                                           units=f'{1/scale_factor} * mm/h',
                                                           time_note=f'start at {metadata_data_000.datum_iso}')

        time_value = 0
        time_unit = f"hours since {metadata_data_000.datum_iso}"
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
             idx_east=None, idx_south=None, idx_north=None):
        """
        Cropping of data of this RadvorRQ class. Usage of indexes directly if given. Otherwise, use lon/lat with the
        guarantee that the whole requested area is within the cropped region. If crop_description is given it is
        prioritized.

        :param crop_description: crop description with some of lon_west, lon_east, lat_south, lat_north, idx_west,
            idx_east, idx_south, idx_north variables
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
        """
        print('cropping')
        if crop_description is not None:
            self.gr_data.crop(lon_west=crop_description.lon_west, lon_east=crop_description.lon_east,
                              lat_south=crop_description.lat_south, lat_north=crop_description.lat_north,
                              idx_west=crop_description.idx_west, idx_east=crop_description.idx_east,
                              idx_south=crop_description.idx_south, idx_north=crop_description.idx_north)
        else:
            self.gr_data.crop(lon_west=lon_west, lon_east=lon_east, lat_south=lat_south, lat_north=lat_north,
                              idx_west=idx_west, idx_east=idx_east, idx_south=idx_south, idx_north=idx_north)

    def export_netcdf(self, filename, data_format='f8', institution=None, scale_factor_nc=1, scale_undo=False,
                      data_kwargs=None):
        """
        Export the relevant content of this RadvorRQ class to a new netcdf file.

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

        if scale_undo and self.gr_data.data_description.scale_factor == 1:
            # roll back scaling called, but not necessary
            scale_undo = False

        # create file and dimensions
        ncfile = nc.Dataset(filename, mode='w', format='NETCDF4')
        ncfile.createDimension(self.nc_desc.dim_lon, self.gr_data.lon.data.shape[1])  # longitude axis
        ncfile.createDimension(self.nc_desc.dim_lat, self.gr_data.lat.data.shape[0])  # latitude axis
        ncfile.createDimension(self.nc_desc.dim_forecast, self.forecast_value.data.size)  # forecast axis
        ncfile.createDimension(self.nc_desc.dim_time, None)  # unlimited time axis
        if not self.gr_data.regridded and not self.gr_data.cropped:
            ncfile.title = 'RadvorRQ data'
        elif self.gr_data.regridded and not self.gr_data.cropped:
            ncfile.title = 'RadvorRQ data (regridded)'
        elif not self.gr_data.regridded and self.gr_data.cropped:
            ncfile.title = 'RadvorRQ data (cropped)'
        else:
            ncfile.title = 'RadvorRQ data (regridded and cropped)'
        ncfile.history = 'v0'
        if institution is not None:
            ncfile.institution = institution
        ncfile.generator = 'generated with weatherDataHarmonizer'
        ncfile.source = 'RadvorRQ data from DWD'
        ncfile.created = dt.datetime.now(tz=dt.timezone.utc).isoformat()

        # create variables in file
        lon = ncfile.createVariable(self.nc_desc.var_lon, np.float32, (self.nc_desc.dim_lat, self.nc_desc.dim_lon),
                                    fill_value=self.gr_data.lon.data_description.fill_value)
        lon.units = self.gr_data.lon.data_description.units
        lon.long_name = self.gr_data.lon.data_description.long_name
        lon.coordinate_system = self.gr_data.lon.data_description.coordinate_system

        lat = ncfile.createVariable(self.nc_desc.var_lat, np.float32, (self.nc_desc.dim_lat, self.nc_desc.dim_lon),
                                    fill_value=self.gr_data.lat.data_description.fill_value)
        lat.units = self.gr_data.lat.data_description.units
        lat.long_name = self.gr_data.lat.data_description.long_name
        lat.coordinate_system = self.gr_data.lat.data_description.coordinate_system

        time = ncfile.createVariable(self.nc_desc.var_time, np.float64, (self.nc_desc.dim_time,))
        time.units = self.time_value.data_description.units
        time.calendar = self.time_value.data_description.calendar

        forecast_time = ncfile.createVariable(self.nc_desc.var_forecast, 'i2', (self.nc_desc.dim_forecast,))
        forecast_time.units = self.forecast_value.data_description.units

        radvorrq = ncfile.createVariable(self.nc_desc.var_data, data_format,
                                         (self.nc_desc.dim_time, self.nc_desc.dim_lat, self.nc_desc.dim_lon,
                                          self.nc_desc.dim_forecast),
                                         fill_value=self.gr_data.data_description.fill_value, **data_kwargs)
        precipitation_units = self.gr_data.data_description.units.split(' ')
        if scale_undo:
            precip_multiplier = float(precipitation_units[0]) * self.gr_data.data_description.scale_factor
            scale_factor_internal = 1
        else:
            precip_multiplier = precipitation_units[0]
            scale_factor_internal = self.gr_data.data_description.scale_factor
        precipitation_units = f'{precip_multiplier} {" ".join(precipitation_units[1:])}'
        radvorrq.units = precipitation_units
        radvorrq.long_name = self.gr_data.data_description.long_name
        radvorrq.scale_factor_internal = scale_factor_internal
        radvorrq.scale_factor = scale_factor_nc

        # write data to variables
        lon[:, :] = self.gr_data.lon.data
        lat[:, :] = self.gr_data.lat.data
        time[:] = self.time_value.data
        forecast_time[:] = self.forecast_value.data
        radvorrq[0, :, :, :] = st.gr_data_scaling(self.gr_data, scale_undo, self.gr_data.data_description.scale_factor,
                                                  scale_factor_nc)

        ncfile.close()

    def export_netcdf_append(self, filename):
        """
        Append the relevant content of this RadvorRQ class to an existing netcdf file.

        :param filename: filename of netcdf file
        :type filename: str
        """
        print('append to existing netcdf')

        # append data to existing file
        if not os.path.exists(filename):
            raise Exception(f'file {filename} does not exist yet')
        ncfile = nc.Dataset(filename, mode='a')

        # check dimensions, shape
        if self.nc_desc.dim_lon not in ncfile.dimensions or \
                self.nc_desc.dim_lat not in ncfile.dimensions or \
                self.nc_desc.dim_time not in ncfile.dimensions or \
                self.nc_desc.dim_forecast not in ncfile.dimensions or \
                self.nc_desc.var_time not in ncfile.variables or \
                self.nc_desc.var_data not in ncfile.variables or \
                self.nc_desc.var_forecast not in ncfile.variables:
            raise Exception(f'either dimensions or variables are not appropriate in {filename}')

        if ncfile[self.nc_desc.var_lon].shape != self.gr_data.lon.data.shape:
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
        ncfile[self.nc_desc.var_data][num_time_nc, :, :, :] = \
            st.gr_data_scaling(self.gr_data, scale_undo, self.gr_data.data_description.scale_factor, scale_factor_nc)

        ncfile.close()
