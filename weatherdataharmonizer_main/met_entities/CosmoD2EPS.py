#   @author: Michael Wagner
#   @organization: TU Dresden
#   @contact: michael.wagner@tu-dresden.de
import os
import copy

import datetime as dt
import numpy as np
import netCDF4 as nc

from met_entities.LonLatTime import LonLatTime
from met_entities.MetEntities import MetEntities
from met_entities.GeoReferencedData import GeoReferencedData
from read_cosmo import readCosmo as readCosmo
from met_entities import VariableDescription as vd
from met_entities.Exceptions import *
from aux_tools import grib2_tools as gt
from aux_tools import scaling_tools as st


class CosmoD2EPS(MetEntities):
    """
    CosmoD2EPS class provides all relevant data of CosmoD2EPS data from DWD.
    """

    def __init__(self,
                 time_value: LonLatTime = None,
                 forecast_value: LonLatTime = None,
                 gr_data=None,
                 eps_member: list = None,
                 short: str = None):
        """
        Initialize CosmoD2EPS class.

        :param time_value: time
        :type time_value: LonLatTime.LonLatTime, optional
        :param forecast_value: forecast time
        :type forecast_value: LonLatTime.LonLatTime, optional
        :param gr_data: list or a single instance of spatial data with list members as
            GeoreferencedData.GeoreferencedData class
        :type gr_data: list or GeoreferencedData.GeoreferencedData, optional
        :param eps_member: list of ensemble members
        :type eps_member: list, optional
        :param short: if int16 or int32, the large data variables are cast to numpy.int16/int32 to minimize memory
            usage; the user must pay attention to the scale_factor to include the necessary precision
        :type short: str, optional
        """
        super().__init__(time_value=time_value, forecast_value=forecast_value, gr_data=gr_data, eps_member=eps_member,
                         short=short)

        # added netcdf dimension and variable description, it can be changed from outside but not with the constructor
        self.nc_desc = vd.NcDimVarDescription(dim_time='time', dim_forecast='forecast_time', dim_lon='lon',
                                              dim_lat='lat', dim_coord='coord', dim_eps='eps', var_time='time',
                                              var_forecast='forecast_time', var_lon='longitude', var_lat='latitude',
                                              var_eps='eps', var_data='cosmod2eps')

    def read_file(self, start_datetime, directory='./', dir_time_descriptor=None, forecast_hours=27, eps_member=None,
                  scale_factor=1, fill_value=np.nan, short=None):
        """
        Reading in all available files of a CosmoD2EPS dataset. If a file for a specific forecast time is not available
        yet, the function automatically stops reading in and delivers a proper datastructure.

        :param start_datetime: time to start reading in CosmoD2EPS files, it is used to build the filenames following
            the convention of the german weather service (DWD)
        :type start_datetime: datetime.datetime
        :param directory: directory with all CosmoD2EPS files from start_datetime
        :type directory: str, optional
        :param dir_time_descriptor: list of datetime.strftime time descriptors for an arbitrary number of additional
            time dependent folders with the data, e.g. ['%Y', '%Y%m%d'] for */yyyy/yyyymmdd/*; the final path is built
            to directory/dir_time_directory[0]/../dir_time_directory[n]/filename
        :type dir_time_descriptor: list, optional
        :param forecast_hours: the number of forecast hours to be read in, default is 27 (max for CosmoD2EPS)
        :type forecast_hours: int, optional
        :param eps_member: list with numbers of eps members to read (0 to 19); if not given, take all 20 realisations
        :type eps_member: list, optional
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

        # define default eps members
        if eps_member is None:
            eps_member = range(20)
        self.eps_member = eps_member

        # get coordinates
        clon, clat = readCosmo.get_lonlat()

        # build filenames and read in data
        filenames = []
        cosmo_data_realis = []
        forecast_value = []
        for i in range(forecast_hours + 1):
            filename = 'cosmo-d2-eps_germany_rotated-lat-lon_single-level_' + start_datetime.strftime('%Y%m%d%H') + \
                       '_' + '{:03d}'.format(i) + '_tot_prec.grib2.bz2'
            print(f'load {filename}')
            if dir_time_descriptor is not None:
                dir_time = ''
                for dir_act in dir_time_descriptor:
                    dir_time = os.path.join(dir_time, start_datetime.strftime(dir_act))
            else:
                dir_time = ''
            filenames.append(os.path.join(directory, dir_time, filename))
            try:
                cosmo_data = readCosmo.read_cosmo_d2(filename=filenames[i], scale_factor=scale_factor,
                                                     fill_value=fill_value, variant='d2eps', short=self.short)

                for realis in range(len(cosmo_data)):
                    if i == 0:
                        cosmo_data_realis.append([cosmo_data[realis]])
                    else:
                        cosmo_data_realis[realis].append(cosmo_data[realis])

            except IOError:
                if i > 0:
                    print(f'Warning: file {filename} not available - stopping at {filenames[i - 1]}')
                else:
                    raise ForecastFileNotAvailable(
                        f'first forecast file {filename} not available - import not successful')
                break

        if len(eps_member) > 1:
            self.gr_data = []
        for member in eps_member:
            # flatten each realisation separately
            cosmo_data_flattened = gt.accum_to_instantaneous_flattened(cosmo_data_realis[member], fill_value=fill_value)
            if len(forecast_value) == 0:
                for cosmo_data_act in cosmo_data_flattened:
                    forecast_value.append(cosmo_data_act.metadata.prediction_time.total_seconds() / 60)
                forecast_description = vd.TimeDescription(units='minutes')

                self.forecast_value = LonLatTime(data=np.array(forecast_value), data_description=forecast_description)

                time_value = 0
                time_unit = f"hours since {cosmo_data_flattened[0].metadata.datum_iso}"
                time_description = vd.TimeDescription(calendar="standard", units=time_unit)
                self.time_value = LonLatTime(data=time_value, data_description=time_description)

            gr_data = GeoReferencedData()
            gr_data.lon = copy.copy(clon)  # copy with another storage address hex(id(gr_data.lon)) != hex(id(clon))
            gr_data.lat = copy.copy(clat)
            if self.short == 'int16':
                gr_data.data = np.empty((cosmo_data_flattened[0].data.size, len(cosmo_data_flattened)), dtype=np.short)
            elif self.short == 'int32':
                gr_data.data = np.empty((cosmo_data_flattened[0].data.size, len(cosmo_data_flattened)), dtype=np.int32)
            else:
                gr_data.data = np.empty((cosmo_data_flattened[0].data.shape + (len(cosmo_data_flattened),)))
            for i in range(len(cosmo_data_flattened)):
                # convert list entries to ndarray
                gr_data.data[:, :, i] = cosmo_data_flattened[i].data
            gr_data.data_description = vd.DataDescription(fill_value=fill_value, scale_factor=scale_factor,
                                                          long_name='CosmoD2EPS precipitation data',
                                                          units=f'{1 / scale_factor} * '
                                                                f'{cosmo_data_flattened[0].metadata.units}',
                                                          time_note=f'start at {cosmo_data_flattened[0].metadata.datum_iso}',
                                                          eps_note=f'eps member {member}')
            if len(eps_member) > 1:
                self.gr_data.append(gr_data)
            else:
                # if only one eps member is available
                self.gr_data = gr_data

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
            if isinstance(self.gr_data, GeoReferencedData):
                # if only one eps member is available
                self.gr_data.regrid_idw(lon_target=regrid_description.lon_target,
                                        lat_target=regrid_description.lat_target,
                                        neighbors=regrid_description.neighbors,
                                        file_nearest=regrid_description.file_nearest,
                                        short=self.short)
            else:
                for i in range(len(self.gr_data)):
                    print(f'eps member = {self.eps_member[i]}')
                    self.gr_data[i].regrid_idw(lon_target=regrid_description.lon_target,
                                               lat_target=regrid_description.lat_target,
                                               neighbors=regrid_description.neighbors,
                                               file_nearest=regrid_description.file_nearest,
                                               short=self.short)
        else:
            if isinstance(self.gr_data, GeoReferencedData):
                # if only one eps member is available
                self.gr_data.regrid_idw(lon_target=lon_target, lat_target=lat_target, neighbors=neighbors,
                                        file_nearest=file_nearest, short=self.short)
            else:
                for i in range(len(self.gr_data)):
                    print(f'eps member = {self.eps_member[i]}')
                    self.gr_data[i].regrid_idw(lon_target=lon_target, lat_target=lat_target, neighbors=neighbors,
                                               file_nearest=file_nearest, short=self.short)

    def crop(self, crop_description=None, lon_west=None, lon_east=None, lat_south=None, lat_north=None, idx_west=None,
             idx_east=None, idx_south=None, idx_north=None):
        """
        Cropping of data of this CosmoD2 class. Usage of indexes directly if given. Otherwise, use lon/lat with the
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
        """
        print('cropping')
        if crop_description is not None:
            if isinstance(self.gr_data, GeoReferencedData):
                # if only one eps member is available
                self.gr_data.crop(lon_west=crop_description.lon_west, lon_east=crop_description.lon_east,
                                  lat_south=crop_description.lat_south, lat_north=crop_description.lat_north,
                                  idx_west=crop_description.idx_west, idx_east=crop_description.idx_east,
                                  idx_south=crop_description.idx_south, idx_north=crop_description.idx_north,
                                  idx_array=crop_description.idx_array)
            else:
                for i in range(len(self.gr_data)):
                    self.gr_data[i].crop(lon_west=crop_description.lon_west, lon_east=crop_description.lon_east,
                                         lat_south=crop_description.lat_south, lat_north=crop_description.lat_north,
                                         idx_west=crop_description.idx_west, idx_east=crop_description.idx_east,
                                         idx_south=crop_description.idx_south, idx_north=crop_description.idx_north,
                                         idx_array=crop_description.idx_array)

        else:
            if isinstance(self.gr_data, GeoReferencedData):
                # if only one eps member is available
                self.gr_data.crop(lon_west=lon_west, lon_east=lon_east, lat_south=lat_south, lat_north=lat_north,
                                  idx_west=idx_west, idx_east=idx_east, idx_south=idx_south, idx_north=idx_north)
            else:
                for i in range(len(self.gr_data)):
                    self.gr_data[i].crop(lon_west=lon_west, lon_east=lon_east, lat_south=lat_south, lat_north=lat_north,
                                         idx_west=idx_west, idx_east=idx_east, idx_south=idx_south, idx_north=idx_north)

    def export_netcdf(self, filename, data_format='f8', version='separated', institution=None, scale_factor_nc=1,
                      scale_undo=False, data_kwargs=None):
        """
        Export the relevant content of this CosmoD2EPS class to a new netcdf file. The function expects the same
        coordinates for all realisations (i.e. the same regridding/cropping for each realisation).

        :param filename: filename of netcdf file
        :type filename: str
        :param data_format: format description of resulting data field in netcdf, the following specifiers are allowed:
            'f8', 'f4', 'i8', 'i4', 'i2', 'i1', 'u8', 'u4', 'u2', 'u1', 'S1' (f - float, i - integer, u - unsigned
            integer, S1 - single character string; the number specifies the number of bytes)
        :type data_format: str, optional
        :param version: the netcdf output is enabled in two different versions: 1. 'separated' (default) with an own
            variable for each realisation, 2. 'integrated' with one large data matrix with eps as the last dimension
        :type version: str, optional
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
        num_eps = len(self.eps_member)

        if num_eps > 1:
            gr_data_first = self.gr_data[0]
        else:
            # if only one eps member is available
            gr_data_first = self.gr_data

        variable_names = []  # only used for separated variables
        if version.lower() == 'separated':
            separation = True
            for realis in self.eps_member:
                variable_names.append(f'{self.nc_desc.var_data}_{realis}')
        else:
            separation = False

        if scale_undo and gr_data_first.data_description.scale_factor == 1:
            # roll back scaling called, but not necessary
            scale_undo = False

        # create file and dimensions
        ncfile = nc.Dataset(filename, mode='w', format='NETCDF4')
        ncfile.createDimension(self.nc_desc.dim_lon, gr_data_first.lon.data.shape[1])  # longitude axis
        ncfile.createDimension(self.nc_desc.dim_lat, gr_data_first.lat.data.shape[0])  # latitude axis
        ncfile.createDimension(self.nc_desc.dim_forecast, self.forecast_value.data.size)  # forecast axis
        ncfile.createDimension(self.nc_desc.dim_time, None)  # unlimited time axis
        if not separation:
            # the eps dimension is only created in 'integrated' version (one large data matrix) and if more than one
            # member is stored
            if num_eps > 1:
                ncfile.createDimension(self.nc_desc.dim_eps, num_eps)

        if not gr_data_first.regridded and not gr_data_first.cropped:
            ncfile.title = 'Cosmo-D2-EPS data'
        elif gr_data_first.regridded and not gr_data_first.cropped:
            ncfile.title = 'Cosmo-D2-EPS data (regridded)'
        elif not gr_data_first.regridded and gr_data_first.cropped:
            ncfile.title = 'Cosmo-D2-EPS data (cropped)'
        else:
            ncfile.title = 'Cosmo-D2-EPS data (regridded and cropped)'
        ncfile.history = 'v0'
        if institution is not None:
            ncfile.institution = institution
        ncfile.generator = 'generated with weatherDataHarmonizer'
        ncfile.source = 'Cosmo-D2-EPS data from DWD'
        ncfile.created = dt.datetime.now(tz=dt.timezone.utc).isoformat()

        # create variables in file
        cosmod2eps = None
        lon = ncfile.createVariable(self.nc_desc.var_lon, np.float32, (self.nc_desc.dim_lat, self.nc_desc.dim_lon),
                                    fill_value=gr_data_first.lon.data_description.fill_value)
        lat = ncfile.createVariable(self.nc_desc.var_lat, np.float32, (self.nc_desc.dim_lat, self.nc_desc.dim_lon),
                                    fill_value=gr_data_first.lat.data_description.fill_value)
        if separation:
            for realis in range(num_eps):
                vars()[variable_names[realis]] = ncfile.createVariable(variable_names[realis], data_format,
                                                                       (self.nc_desc.dim_time,
                                                                        self.nc_desc.dim_lat,
                                                                        self.nc_desc.dim_lon,
                                                                        self.nc_desc.dim_forecast),
                                                                       fill_value=gr_data_first
                                                                       .data_description.fill_value,
                                                                       **data_kwargs)
        else:
            if num_eps > 1:
                cosmod2eps = ncfile.createVariable(self.nc_desc.var_data, data_format,
                                                   (self.nc_desc.dim_time, self.nc_desc.dim_lat,
                                                    self.nc_desc.dim_lon, self.nc_desc.dim_forecast,
                                                    self.nc_desc.dim_eps),
                                                   fill_value=gr_data_first.data_description.fill_value,
                                                   **data_kwargs)
            else:
                # if only one eps member is available
                cosmod2eps = ncfile.createVariable(self.nc_desc.var_data, data_format,
                                                   (self.nc_desc.dim_time, self.nc_desc.dim_lat,
                                                    self.nc_desc.dim_lon, self.nc_desc.dim_forecast),
                                                   fill_value=gr_data_first.data_description.fill_value,
                                                   **data_kwargs)

        lon.units = gr_data_first.lon.data_description.units
        lon.long_name = gr_data_first.lon.data_description.long_name
        lon.coordinate_system = gr_data_first.lon.data_description.coordinate_system

        lat.units = gr_data_first.lat.data_description.units
        lat.long_name = gr_data_first.lat.data_description.long_name
        lat.coordinate_system = gr_data_first.lat.data_description.coordinate_system

        time = ncfile.createVariable(self.nc_desc.var_time, np.float64, (self.nc_desc.dim_time,))
        time.units = self.time_value.data_description.units
        time.calendar = self.time_value.data_description.calendar

        forecast_time = ncfile.createVariable(self.nc_desc.var_forecast, 'i2', (self.nc_desc.dim_forecast,))
        forecast_time.units = self.forecast_value.data_description.units

        precipitation_units = gr_data_first.data_description.units.split(' ')
        if scale_undo:
            precip_multiplier = float(precipitation_units[0]) * gr_data_first.data_description.scale_factor
            scale_factor_internal = 1
        else:
            precip_multiplier = precipitation_units[0]
            scale_factor_internal = gr_data_first.data_description.scale_factor
        precipitation_units = f'{precip_multiplier} {" ".join(precipitation_units[1:])}'

        if separation:
            for realis in range(num_eps):
                vars()[variable_names[realis]].units = precipitation_units
                vars()[variable_names[realis]].long_name = gr_data_first.data_description.long_name
                vars()[variable_names[realis]].scale_factor_internal = scale_factor_internal
                vars()[variable_names[realis]].scale_factor = scale_factor_nc
        else:
            cosmod2eps.units = precipitation_units
            cosmod2eps.long_name = gr_data_first.data_description.long_name
            cosmod2eps.scale_factor_internal = scale_factor_internal
            cosmod2eps.scale_factor = scale_factor_nc

        # write data to variables
        lon[:, :] = gr_data_first.lon.data
        lat[:, :] = gr_data_first.lat.data
        if separation:
            if isinstance(self.gr_data, GeoReferencedData):
                # if only one eps member is available
                vars()[variable_names[0]][0, :, :, :] = \
                    st.gr_data_scaling(self.gr_data, scale_undo, self.gr_data.data_description.scale_factor,
                                       scale_factor_nc)
            else:
                for realis in range(num_eps):
                    vars()[variable_names[realis]][0, :, :, :] = \
                        st.gr_data_scaling(self.gr_data[realis], scale_undo,
                                           self.gr_data[realis].data_description.scale_factor, scale_factor_nc)
        else:
            if isinstance(self.gr_data, GeoReferencedData):
                # if only one eps member is available
                cosmod2eps[0, :, :, :] = st.gr_data_scaling(self.gr_data, scale_undo,
                                                            self.gr_data.data_description.scale_factor,
                                                            scale_factor_nc)
            else:
                ct = -1
                for gr_data in self.gr_data:
                    ct = ct + 1
                    cosmod2eps[0, :, :, :, ct] = st.gr_data_scaling(gr_data, scale_undo,
                                                                    gr_data.data_description.scale_factor,
                                                                    scale_factor_nc)
        time[:] = self.time_value.data
        forecast_time[:] = self.forecast_value.data
        if not separation and num_eps > 1:
            eps = ncfile.createVariable(self.nc_desc.var_eps, 'i2', (self.nc_desc.dim_eps,))
            eps.description = 'ensemble member number'
            eps[:] = self.eps_member

        ncfile.close()

    def export_netcdf_append(self, filename):
        """
        Append the relevant content of this CosmoD2EPS class to an existing netcdf file.

        :param filename: filename of netcdf file
        :type filename: str
        """
        print('append to existing netcdf')

        num_eps = len(self.eps_member)
        if num_eps > 1:
            gr_data_first = self.gr_data[0]
        else:
            # if only one eps member is available
            gr_data_first = self.gr_data

        # append data to existing file
        if not os.path.exists(filename):
            raise Exception(f'file {filename} does not exist yet')
        ncfile = nc.Dataset(filename, mode='a')

        # get variable names and separation flag
        variable_names = [key for key in ncfile.variables.keys() if key.startswith(self.nc_desc.var_data)]
        if variable_names[0] == self.nc_desc.var_data:
            separation = False
        else:
            separation = True

        # check dimensions, shape
        if ((self.nc_desc.dim_lon not in ncfile.dimensions or self.nc_desc.dim_lat not in ncfile.dimensions) and
            self.nc_desc.dim_coord not in ncfile.dimensions) or \
                self.nc_desc.dim_time not in ncfile.dimensions or \
                self.nc_desc.dim_forecast not in ncfile.dimensions or \
                self.nc_desc.var_time not in ncfile.variables or \
                self.nc_desc.var_forecast not in ncfile.variables:
            raise Exception(f'either dimensions or variables are not appropriate in {filename}')
        if separation:
            for variable_name in variable_names:
                if variable_name not in ncfile.variables:
                    raise Exception(f'data variables are not appropriate in {filename}')
        else:
            if self.nc_desc.var_data not in ncfile.variables:
                raise Exception(f'data variables are not appropriate in {filename}')
            if num_eps > 1 and (
                    self.nc_desc.dim_eps not in ncfile.dimensions or self.nc_desc.var_eps not in ncfile.variables):
                raise Exception(f'either dimensions or variables are not appropriate in {filename}')

        if ncfile[self.nc_desc.var_lon].shape != gr_data_first.lon.data.shape:
            raise Exception(f'the data shape in {filename} differs from object shape')
        if ncfile[self.nc_desc.var_forecast].size != self.forecast_value.data.size:
            raise Exception(f'the forecast shape in {filename} differs from object shape')

        # check scale_factor and fill_value
        scale_undo = False
        if separation:
            if ncfile.variables[variable_names[0]].scale_factor_internal == 1 and gr_data_first.data_description.scale_factor != 1:
                scale_undo = True
            if ncfile.variables[variable_names[0]].scale_factor_internal != 1 and \
                    gr_data_first.data_description.scale_factor != ncfile[variable_names[0]].scale_factor_internal:
                raise Exception('internal scale factor in netcdf file differs from actual internal scale factor')
            nc_fill_value = ncfile[variable_names[0]]._FillValue
            if np.isnan(gr_data_first.data_description.fill_value) and not np.isnan(nc_fill_value):
                raise Exception('fill value in netcdf file differs from actual fill value')
            elif not np.isnan(gr_data_first.data_description.fill_value) and \
                    gr_data_first.data_description.fill_value != nc_fill_value:
                raise Exception('fill value in netcdf file differs from actual fill value')
        else:
            if ncfile.variables[self.nc_desc.var_data].scale_factor_internal == 1 and gr_data_first.data_description.scale_factor != 1:
                scale_undo = True
            if ncfile.variables[self.nc_desc.var_data].scale_factor_internal != 1 and \
                    gr_data_first.data_description.scale_factor != ncfile[self.nc_desc.var_data].scale_factor_internal:
                raise Exception('internal scale factor in netcdf file differs from actual internal scale factor')
            nc_fill_value = ncfile[self.nc_desc.var_data]._FillValue
            if np.isnan(gr_data_first.data_description.fill_value) and not np.isnan(nc_fill_value):
                raise Exception('fill value in netcdf file differs from actual fill value')
            elif not np.isnan(gr_data_first.data_description.fill_value) and \
                    gr_data_first.data_description.fill_value != nc_fill_value:
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
        if separation:
            scale_factor_nc = ncfile[variable_names[0]].scale_factor
        else:
            scale_factor_nc = ncfile[self.nc_desc.var_data].scale_factor
        if separation:
            if isinstance(self.gr_data, GeoReferencedData):
                # if only one eps member is available
                ncfile[variable_names[0]][num_time_nc, :, :, :] = \
                    st.gr_data_scaling(self.gr_data, scale_undo, self.gr_data.data_description.scale_factor,
                                       scale_factor_nc)
            else:
                for realis in range(num_eps):
                    ncfile[variable_names[realis]][num_time_nc, :, :, :] = \
                        st.gr_data_scaling(self.gr_data[realis], scale_undo,
                                           self.gr_data[realis].data_description.scale_factor, scale_factor_nc)
        else:
            if isinstance(self.gr_data, GeoReferencedData):
                # if only one eps member is available
                ncfile[self.nc_desc.var_data][num_time_nc, :, :, :] = \
                    st.gr_data_scaling(self.gr_data, scale_undo, self.gr_data.data_description.scale_factor,
                                       scale_factor_nc)
            else:
                for realis in range(num_eps):
                    ncfile[self.nc_desc.var_data][num_time_nc, :, :, :, realis] = \
                        st.gr_data_scaling(self.gr_data[realis], scale_undo,
                                           self.gr_data[realis].data_description.scale_factor, scale_factor_nc)

        ncfile.close()
