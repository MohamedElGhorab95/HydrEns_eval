#   @author: Michael Wagner
#   @organization: TU Dresden
#   @contact: michael.wagner@tu-dresden.de
import os
import numpy as np
import netCDF4 as nc
import datetime as dt

from met_entities.GeoReferencedData import GeoReferencedData
from met_entities.MetEntities import MetEntities
from read_radolan import readRadolan as readRadolan
from met_entities import VariableDescription as vd
from met_entities import LonLatTime as LonLatTime
from met_entities.Exceptions import *
from aux_tools import scaling_tools as st


class RadolanRW(MetEntities):
    """
    RadolanRW class provides all relevant data of RadolanRW data from DWD.
    """

    def __init__(self,
                 time_value: LonLatTime = None,
                 gr_data: GeoReferencedData = None,
                 short: str = None):
        """
        Initialize RadolanRW class.

        :param time_value: time
        :type time_value: LonLatTime.LonLatTime, optional
        :param gr_data: spatial data
        :type gr_data: GeoReferencedData.GeoReferencedData, optional
        :param short: if int16 or int32, the large data variables are cast to numpy.int16/int32 to minimize memory
            usage; the user must pay attention to the scale_factor to include the necessary precision
        :type short: str, optional
        """
        super().__init__(time_value=time_value, gr_data=gr_data, short=short)

        # added netcdf dimension and variable description, it can be changed from outside but not with the constructor
        self.nc_desc = vd.NcDimVarDescription(dim_time='time', dim_lon='lon', dim_lat='lat',
                                              var_time='time', var_lon='longitude', var_lat='latitude',
                                              var_data='radolanrw')

    def read_file(self, filename=None, start_datetime=None, directory='./', dir_time_descriptor=None, scale_factor=1,
                  fill_value=np.nan, short=None):
        """
        Reading in a RadolanRW file.

        :param filename: RW file, typically in the form raa01-rw_10000-yymmddHHMM-dwd---bin.bz2; if not given
            start_datetime and directory elements must be provided
        :type filename: str, optional
        :param start_datetime: time to start reading in RadolanRW files, it is used to build the filenames following the
            convention of the german weather service (DWD); can be given in xx:50 or xx:00 (full hours), in the second
            case a -10 minutes shift is done internally to meet the RadolanRW delivering convention
        :type start_datetime: datetime.datetime, optional
        :param directory: directory with the RadolanRW file at start_datetime
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

        # shift start_datetime by 10 min as radolanrw data delivered at xx:50
        if start_datetime is not None:
            if start_datetime.minute == 50:
                # correct observation time
                start_datetime_radolan_shift = start_datetime
            elif start_datetime.minute == 0:
                # full hours
                start_datetime_radolan_shift = start_datetime - dt.timedelta(minutes=10)
            else:
                raise Exception(f'only 50 or 0 minutes can be interpreted')
        else:
            start_datetime_radolan_shift = None

        # get data and metadata
        radolan_data = None
        if filename is not None:
            try:
                print(f'load {filename}')
                radolan_data = readRadolan.read_radolan(filename=filename, scale_factor=scale_factor,
                                                        fill_value=fill_value, short=self.short)
            except FileNotFoundError:
                raise RadarFileNotAvailable(f'file {filename} not available - import not successful')
        else:
            # rely on start_datetime and directory and three typical file names and extensions
            filename_extensions = ['', '.gz', '.bz2']
            if dir_time_descriptor is not None:
                dir_time = ''
                for dir_act in dir_time_descriptor:
                    dir_time = os.path.join(dir_time, start_datetime.strftime(dir_act))
            else:
                dir_time = ''
            filename_base = os.path.join(directory, dir_time, 'raa01-rw_10000-' + start_datetime_radolan_shift.strftime(
                '%y%m%d%H%M') + '-dwd---bin')
            data_loaded = False
            for filename_extension in filename_extensions:
                filename_act = filename_base + filename_extension
                if os.path.exists(filename_act):
                    print(f'load {filename_act}')
                    radolan_data = readRadolan.read_radolan(filename=filename_act, scale_factor=scale_factor,
                                                            fill_value=fill_value, short=self.short)
                    data_loaded = True
                    break
            if not data_loaded:
                raise RadarFileNotAvailable(f'file {filename} not available - import not successful')

        self.gr_data.lon, self.gr_data.lat = readRadolan.get_lonlat(radolan_data.metadata.format_version,
                                                                    grid_variant='radolanrx')
        self.gr_data.data = radolan_data.data
        self.gr_data.data_description = vd.DataDescription(fill_value=fill_value,
                                                           scale_factor=scale_factor,
                                                           long_name='RadolanRW precipitation data',
                                                           units=f'{1/scale_factor} * mm/h',
                                                           time_note=f'start at {radolan_data.metadata.datum_iso}')

        time_value = 0
        time_unit = f"hours since {radolan_data.metadata.datum_iso}"
        time_description = vd.TimeDescription(calendar="standard", units=time_unit)
        self.time_value = LonLatTime.LonLatTime(data=time_value, data_description=time_description)

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
        Export the relevant content of this RadolanRW class to a new netcdf file.

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
        ncfile.createDimension(self.nc_desc.dim_time, None)  # unlimited time axis
        if not self.gr_data.regridded and not self.gr_data.cropped:
            ncfile.title = 'RadolanRW data'
        elif self.gr_data.regridded and not self.gr_data.cropped:
            ncfile.title = 'RadolanRW data (regridded)'
        elif not self.gr_data.regridded and self.gr_data.cropped:
            ncfile.title = 'RadolanRW data (cropped)'
        else:
            ncfile.title = 'RadolanRW data (regridded and cropped)'
        ncfile.history = 'v0'
        if institution is not None:
            ncfile.institution = institution
        ncfile.generator = 'generated with weatherDataHarmonizer'
        ncfile.source = 'RadolanRW data from DWD'
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

        radolanrw = ncfile.createVariable(self.nc_desc.var_data, data_format,
                                          (self.nc_desc.dim_time, self.nc_desc.dim_lat, self.nc_desc.dim_lon),
                                          fill_value=self.gr_data.data_description.fill_value, **data_kwargs)
        precipitation_units = self.gr_data.data_description.units.split(' ')
        if scale_undo:
            precip_multiplier = float(precipitation_units[0]) * self.gr_data.data_description.scale_factor
            scale_factor_internal = 1
        else:
            precip_multiplier = precipitation_units[0]
            scale_factor_internal = self.gr_data.data_description.scale_factor
        precipitation_units = f'{precip_multiplier} {" ".join(precipitation_units[1:])}'
        radolanrw.units = precipitation_units
        radolanrw.long_name = self.gr_data.data_description.long_name
        radolanrw.scale_factor_internal = scale_factor_internal
        radolanrw.scale_factor = scale_factor_nc

        # write data to variables
        lon[:, :] = self.gr_data.lon.data
        lat[:, :] = self.gr_data.lat.data
        time[:] = self.time_value.data
        radolanrw[0, :, :] = st.gr_data_scaling(self.gr_data, scale_undo, self.gr_data.data_description.scale_factor,
                                                   scale_factor_nc)

        ncfile.close()

    def export_netcdf_append(self, filename):
        """
        Append the relevant content of this RadolanRW class to an existing netcdf file.

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
                self.nc_desc.var_time not in ncfile.variables or \
                self.nc_desc.var_data not in ncfile.variables:
            raise Exception(f'either dimensions or variables are not appropriate in {filename}')

        if ncfile[self.nc_desc.var_lon].shape != self.gr_data.lon.data.shape:
            raise Exception(f'the data shape in {filename} differs from object shape')

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
        ncfile[self.nc_desc.var_data][num_time_nc, :, :] = \
            st.gr_data_scaling(self.gr_data, scale_undo, self.gr_data.data_description.scale_factor, scale_factor_nc)

        ncfile.close()

    def import_netcdf(self, filename, nc_desc=None):
        """
        Import a netcdf file to instantiate a RadolanRW object.

        :param filename: name of the netcdf file.
        :type filename: str
        :param nc_desc: detailed description of dimension and variable names in the netcdf file; if omitted the standard
            is taken (see constructor __init__)
        :type nc_desc: vd.NcDimVarDescription, optional
        """
        print('import netcdf')

        # read netcdf file and check dimension and variable availability
        if self.gr_data is not None or self.time_value is not None:
            print('this RadolanRW instance is not empty, available content will be overwritten')
        if nc_desc is not None:
            self.nc_desc = nc_desc
        ncfile = nc.Dataset(filename, mode='r')

        nc_dims = ncfile.dimensions
        for nc_dim in nc_dims:
            if nc_dim not in [self.nc_desc.dim_lon, self.nc_desc.dim_lat, self.nc_desc.dim_time]:
                raise Exception(f'in file {filename} the dimension {nc_dim.name} cannot be interpreted')

        nc_vars = ncfile.variables
        if self.nc_desc.var_lon not in nc_vars:
            raise Exception(f'variable {self.nc_desc.var_lon} not in file {filename}')
        if self.nc_desc.var_lat not in nc_vars:
            raise Exception(f'variable {self.nc_desc.var_lat} not in file {filename}')
        if self.nc_desc.var_time not in nc_vars:
            raise Exception(f'variable {self.nc_desc.var_time} not in file {filename}')
        if self.nc_desc.var_data not in nc_vars:
            raise Exception(f'variable {self.nc_desc.var_data} not in file {filename}')

        # fill self variables
        time_description = vd.TimeDescription(calendar=nc_vars[self.nc_desc.var_time].calendar,
                                              units=nc_vars[self.nc_desc.var_time].units)
        self.time_value = LonLatTime.LonLatTime(data=nc_vars[self.nc_desc.var_time][:],
                                                data_description=time_description)

        self.gr_data = GeoReferencedData()
        lon_description = vd.DataDescription(units=nc_vars[self.nc_desc.var_lon].units,
                                             long_name=nc_vars[self.nc_desc.var_lon].long_name,
                                             coordinate_system=nc_vars[self.nc_desc.var_lon].coordinate_system)
        self.gr_data.lon = LonLatTime.LonLatTime(data=nc_vars[self.nc_desc.var_lon][:],
                                                 data_description=lon_description)
        lat_description = vd.DataDescription(units=nc_vars[self.nc_desc.var_lat].units,
                                             long_name=nc_vars[self.nc_desc.var_lat].long_name,
                                             coordinate_system=nc_vars[self.nc_desc.var_lat].coordinate_system)
        self.gr_data.lat = LonLatTime.LonLatTime(data=nc_vars[self.nc_desc.var_lat][:],
                                                 data_description=lat_description)
        self.gr_data.data = np.squeeze(nc_vars[self.nc_desc.var_data])
        self.gr_data.data_description = vd.DataDescription(fill_value=nc_vars[self.nc_desc.var_data]._FillValue,
                                                           scale_factor=nc_vars[self.nc_desc.var_data].scale_factor,
                                                           long_name=nc_vars[self.nc_desc.var_data].long_name,
                                                           units=nc_vars[self.nc_desc.var_data].units)

        ncfile.close()
