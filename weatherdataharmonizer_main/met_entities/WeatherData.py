#   @author: Michael Wagner
#   @organization: TU Dresden
#   @contact: michael.wagner@tu-dresden.de

import os
import datetime as dt
import numpy as np
import netCDF4 as nc
import copy

from aux_tools import scaling_tools as st
from met_entities.GeoReferencedData import GeoReferencedData
from met_entities.IconD2 import IconD2
from met_entities.IconD2EPS import IconD2EPS
from met_entities.IconEU import IconEU
from met_entities.IconEUEPS import IconEUEPS
from met_entities.RadolanRV import RadolanRV
from met_entities.RadolanRW import RadolanRW
from met_entities.RadvorRQ import RadvorRQ
import met_entities.VariableDescription as vd
from met_entities.Exceptions import *
from aux_tools.products import *
from met_entities.read_netcdf import read_netcdf


class WeatherData:
    """
    WeatherData class is a container for mixed weather data from radolan, radvor, icond2, icond2eps, iconeu, iconeueps
    data sources. It comprises methods to sort the date in a chronological order from past to future.
    """
    def __init__(self,
                 time_now: dt.datetime,
                 delta_t: dt.timedelta,
                 data=None,
                 scale_factor=1,
                 fill_value=-999,
                 short=None,
                 regrid_description: dict = None,
                 crop_description: vd.CropDescription = None):
        """
        Initialize WeatherData class.

        :param time_now: time until when observed data shall be used if available, otherwise it is changed to an earlier
            time
        :type time_now: datetime.datetime
        :param delta_t: time step length; valid values comprise 5 min, 15 min, 60 min, and 120 min
        :type delta_t: datetime.timedelta
        :param data: collected data (GeoreferencedData instances) with content filled from all collect_xxx methods; all
            list entries have the same size with central regridding or cropping descriptions; regridding is way better
            suited to ensure the same grid in the result as cropping does not change the original grid base and further
            depends on the overall grid because of bent row and line coordinates
        :type data: list
        :param scale_factor: the final data has to be multiplied with this value
        :type scale_factor: float, optional
        :param fill_value: missing data is filled with that value
        :type fill_value: float, optional
        :param short: if int16 or int32, the large data variables are cast to numpy.int16/int32 to minimize memory
            usage; the user must pay attention to the scale_factor to include the necessary precision
        :type short: str, optional
        :param regrid_description: the complete description for regridding the data; it shall consist of a dictionary
            with the product's name as key (str) and the regrid description as value
            (VariableDescription.RegridDescription)
        :type regrid_description: dict, optional
        :param crop_description: the complete description for cropping the data
        :type crop_description: VariableDescription.CropDescription, optional
        """
        if data is None:
            data = []
        self.time_now = time_now
        self.delta_t = delta_t
        self.data = data
        self.scale_factor = scale_factor
        self.fill_value = fill_value
        self.short = short
        self.regrid_description = regrid_description
        self.crop_description = crop_description

        allowed_dt = [dt.timedelta(minutes=t) for t in [5, 15, 60, 120]]
        if self.delta_t not in allowed_dt:
            raise Exception(f'currently, only delta t {allowed_dt} are supported')

        if self.time_now.minute != 0:
            raise Exception('only full hour datetimes are allowed as time_now')

        self.time = []
        self.time_protected = []
        self.radolanrw_loaded = False
        self.radvorrq_loaded = False
        self.radolanrv_loaded = False
        self.icond2_loaded = False
        self.icond2eps_loaded = False
        self.iconeu_loaded = False
        self.iconeueps_loaded = False
        self.eps_member = []
        self.eps_member_source = None

        # added netcdf dimension and variable description, it can be changed from outside but not with the constructor
        self.nc_desc = vd.NcDimVarDescription(dim_time='time', dim_lon='lon', dim_lat='lat',
                                              var_time='time', var_lon='longitude', var_lat='latitude',
                                              var_data='precipitation')
        self.data_source_dict = {'radolanrw': 0, 'radvorrq': 1, 'radolanrv': 2, 'icond2': 3, 'icond2eps': 4,
                                 'iconeu': 5, 'iconeueps': 6}

    def collect_radolanrw(self,
                          time_start,
                          directory=None,
                          dir_time_descriptor=None):
        """
        Collect all RadolanRW data into the variable data and construct the according time vector.

        :param time_start: time to start with RadolanRW data; it must avoid xx:50, instead give full hours; the time
            shift is done function-internal
        :type time_start: datetime.datetime
        :param directory: directory with RadolanRW data
        :type directory: str, optional
        :param dir_time_descriptor: list of datetime.strftime time descriptors for an arbitrary number of additional
            time dependent folders with the data, e.g. ['%Y', '%Y%m%d'] for */yyyy/yyyymmdd/*; the final path is built
            to directory/dir_time_directory[0]/../dir_time_directory[n]/filename
        :type dir_time_descriptor: list, optional
        """
        if self.radolanrw_loaded:
            print('RadolanRW already loaded; nothing will be done')
            return
        if time_start > self.time_now:
            raise Exception(f'time_start ({time_start}) must not be later than time_now ({self.time_now})')
        if any([self.radvorrq_loaded, self.radolanrv_loaded, self.icond2_loaded, self.icond2eps_loaded]):
            raise Exception('prediction data must not be loaded before radolan data')
        if time_start.minute != 0 and time_start.minute != 50:
            raise Exception(f'found starting time minutes {time_start.minute} - for radolanrw data it must be 0 or 50')
        if ((self.time_now - time_start) % self.delta_t).total_seconds() != 0:
            raise Exception(
                f'time_now {self.time_now} is not reachable from time_start {time_start} with dt {self.delta_t}')
        if time_start.minute == 0:
            time_start_tmp = time_start - dt.timedelta(minutes=10)
        else:
            time_start_tmp = time_start
        filename_base = os.path.join(directory, 'raa01-rw_10000-' + time_start_tmp.strftime('%y%m%d%H%M') + '-dwd---bin')
        if not os.path.exists(filename_base) and not os.path.exists(filename_base + '.gz') and \
                not os.path.exists(filename_base + '.bz2'):
            raise RadarFileNotAvailable(
                f'even the first RadolanRW file at {time_start_tmp} is not available - no data loading possible')

        print('collect RadolanRW data')

        time_start_act = time_start
        if self.delta_t == dt.timedelta(minutes=5):
            # add a time step after time_now to get 5-minute values until full hour
            latest_event = self.time_now + dt.timedelta(minutes=50)
            if time_start.minute == 0:
                time_start_act = time_start - dt.timedelta(minutes=10)
        else:
            # by default replace xx:50 with full hour xx:00
            latest_event = self.time_now
            if time_start.minute == 50:
                time_start_act = time_start + dt.timedelta(minutes=10)

        data, time, time_now_changed = collect_fc(product='radolanrw', directory=directory,
                                                  dir_time_descriptor=dir_time_descriptor, time_now=time_start_act,
                                                  latest_event=latest_event, update_time=RadolanRWProp().update_time,
                                                  delta_t=self.delta_t, regrid_description=self.regrid_description,
                                                  crop_description=self.crop_description,
                                                  scale_factor=self.scale_factor,
                                                  fill_value=self.fill_value, short=self.short)

        # if time_now is not covered by radolanrw change time_now to earlier time step
        if time_now_changed is not None:
            self.time_now = time_now_changed
        # include only time steps from time_start to self.time_now
        if self.time_now in time:
            idx = time.index(self.time_now)
            time = time[0:idx + 1]
            data = data[0:idx + 1]
        # else:
            # if time_now is not covered by radolanrw change time_now to earlier time step
            # self.time_now = time[-1]
            # print(f'time_now changed to {self.time_now} due to missing radolanrw data')
        if time_start in time:
            idx = time.index(time_start)
            time = time[idx:]
            data = data[idx:]
        self.time = time
        self.time_protected = time.copy()  # protected time that must not be replaced (observed values)
        self.data = data

        self.radolanrw_loaded = True

    def collect_radvorrq(self,
                         latest_event,
                         directory=None,
                         dir_time_descriptor=None):
        """
        Collect all RadvorRQ data into the variable data and construct the according time vector. collect_radolanrw must
        be invoked before collecting prediction data.

        :param latest_event: datetime of the latest forecast start to be tried to load
        :type latest_event: datetime.datetime
        :param directory: directory with RadvorRQ data
        :type directory: str, optional
        :param dir_time_descriptor: list of datetime.strftime time descriptors for an arbitrary number of additional
            time dependent folders with the data, e.g. ['%Y', '%Y%m%d'] for */yyyy/yyyymmdd/*; the final path is built
            to directory/dir_time_directory[0]/../dir_time_directory[n]/filename
        :type dir_time_descriptor: list, optional
        """
        if not self.radolanrw_loaded:
            raise Exception('you must first load RadolanRW data before loading forecast data')
        if self.radvorrq_loaded:
            print('RadvorRQ already loaded; nothing will be done')
            return

        update_time = RadvorRQProp().update_time
        update_time_in_minutes = update_time.total_seconds() / 60
        if latest_event.minute % update_time_in_minutes != 0:
            # latest event does not fit to update cycle replace the latest event to the next future one
            latest_event = latest_event + dt.timedelta(minutes=(update_time_in_minutes -
                                                                (latest_event.minute % update_time_in_minutes)))

        print('collect RadvorRQ data')

        fc_data, fc_time, _ = collect_fc(product='radvorrq', directory=directory,
                                         dir_time_descriptor=dir_time_descriptor, time_now=self.time_now,
                                         delta_t=self.delta_t, latest_event=latest_event,
                                         fc_max_time=dt.timedelta(hours=2), update_time=update_time,
                                         regrid_description=self.regrid_description,
                                         crop_description=self.crop_description,
                                         scale_factor=self.scale_factor, fill_value=self.fill_value, short=self.short)
        if fc_time is None:
            # no RadvorRQ data available
            return
        self.amend_time_data(time=fc_time, data=fc_data)
        self.radvorrq_loaded = True

    def collect_radolanrv(self,
                          latest_event,
                          directory=None,
                          dir_time_descriptor=None):
        """
        Collect all RadolanRV data into the variable data and construct the according time vector.
        collect_radolanrw must be invoked before collecting prediction data.

        :param latest_event: datetime of the latest forecast start to be tried to load
        :type latest_event: datetime.datetime
        :param directory: directory with RadolanRV data
        :type directory: str, optional
        :param dir_time_descriptor: list of datetime.strftime time descriptors for an arbitrary number of additional
            time dependent folders with the data, e.g. ['%Y', '%Y%m%d'] for */yyyy/yyyymmdd/*; the final path is built
            to directory/dir_time_directory[0]/../dir_time_directory[n]/filename
        :type dir_time_descriptor: list, optional
        """
        if not self.radolanrw_loaded:
            raise Exception('you must first load RadolanRW data before loading forecast data')
        if self.radolanrv_loaded:
            print('RadolanRV already loaded; nothing will be done')
            return

        update_time = RadolanRVProp().update_time
        update_time_in_minutes = update_time.total_seconds() / 60
        if latest_event.minute % update_time_in_minutes != 0:
            # latest event does not fit to update cycle replace the latest event to the next future one
            latest_event = latest_event + dt.timedelta(minutes=(update_time_in_minutes -
                                                                (latest_event.minute % update_time_in_minutes)))

        print('collect RadolanRV data')

        fc_data, fc_time, _ = collect_fc(product='radolanrv', directory=directory,
                                         dir_time_descriptor=dir_time_descriptor, time_now=self.time_now,
                                         delta_t=self.delta_t, latest_event=latest_event,
                                         fc_max_time=dt.timedelta(hours=2), update_time=update_time,
                                         regrid_description=self.regrid_description,
                                         crop_description=self.crop_description,
                                         scale_factor=self.scale_factor, fill_value=self.fill_value, short=self.short)
        if fc_time is None:
            # no RadolanRV data available
            return
        self.amend_time_data(time=fc_time, data=fc_data)
        self.radolanrv_loaded = True


    def collect_icond2(self,
                       latest_event,
                       directory=None,
                       dir_time_descriptor=None,
                       forecast_hours=48):
        """
        Collect all IconD2 data into the variable data and construct the according time vector. collect_radolanrw must
        be invoked before collecting prediction data.

        :param latest_event: datetime of the latest forecast start to be tried to load
        :type latest_event: datetime.datetime
        :param directory: directory with IconD2 data
        :type directory: str, optional
        :param dir_time_descriptor: list of datetime.strftime time descriptors for an arbitrary number of additional
            time dependent folders with the data, e.g. ['%Y', '%Y%m%d'] for */yyyy/yyyymmdd/*; the final path is built
            to directory/dir_time_directory[0]/../dir_time_directory[n]/filename
        :type dir_time_descriptor: list, optional
        :param forecast_hours: try to load this number of hours at maximum
        :type forecast_hours: int, optional
        """
        if not self.radolanrw_loaded:
            raise Exception('you must first load RadolanRW data before loading forecast data')
        if self.icond2_loaded:
            print('IconD2 already loaded; nothing will be done')
            return

        update_time = IconD2Prop().update_time
        update_time_in_hours = update_time.total_seconds()/3600
        if latest_event.hour % update_time_in_hours != 0:
            # latest event does not fit to update cycle replace the latest event to the next future one
            latest_event = latest_event + dt.timedelta(hours=(update_time_in_hours -
                                                              (latest_event.hour % update_time_in_hours)))

        print('collect IconD2 data')

        fc_data, fc_time, _ = collect_fc(product='icond2', directory=directory, dir_time_descriptor=dir_time_descriptor,
                                         time_now=self.time_now, delta_t=self.delta_t, latest_event=latest_event,
                                         fc_max_time=dt.timedelta(hours=forecast_hours), update_time=update_time,
                                         regrid_description=self.regrid_description,
                                         crop_description=self.crop_description, scale_factor=self.scale_factor,
                                         fill_value=self.fill_value, short=self.short)
        if fc_time is None:
            # no IconD2 data available
            return
        self.amend_time_data(time=fc_time, data=fc_data)
        self.icond2_loaded = True

    def collect_icond2eps(self,
                          latest_event,
                          directory=None,
                          dir_time_descriptor=None,
                          forecast_hours=48,
                          eps_member=range(20)):
        """
        Collect all IconD2EPS data for one ensemble member into the variable data and construct the according time
        vector. collect_radolanrw must be invoked before collecting prediction data.

        :param latest_event: datetime of the latest forecast start to be tried to load
        :type latest_event: datetime.datetime
        :param directory: directory with IconD2EPS data
        :type directory: str, optional
        :param dir_time_descriptor: list of datetime.strftime time descriptors for an arbitrary number of additional
            time dependent folders with the data, e.g. ['%Y', '%Y%m%d'] for */yyyy/yyyymmdd/*; the final path is built
            to directory/dir_time_directory[0]/../dir_time_directory[n]/filename
        :type dir_time_descriptor: list, optional
        :param forecast_hours: try to load this number of hours at maximum
        :type forecast_hours: int, optional
        :param eps_member: list of eps members to load (0 to 19)
        :type eps_member: list, optional
        """
        if not self.radolanrw_loaded:
            raise Exception('you must first load RadolanRW data before loading forecast data')
        if self.icond2eps_loaded:
            print('IconD2EPS already loaded; nothing will be done')
            return

        update_time = IconD2EPSProp().update_time
        update_time_in_hours = update_time.total_seconds() / 3600
        if latest_event.hour % update_time_in_hours != 0:
            # latest event does not fit to update cycle replace the latest event to the next future one
            latest_event = latest_event + dt.timedelta(hours=(update_time_in_hours -
                                                              (latest_event.hour % update_time_in_hours)))

        print('collect IconD2EPS data')

        fc_data, fc_time, _ = collect_fc(product='icond2eps', directory=directory,
                                         dir_time_descriptor=dir_time_descriptor, time_now=self.time_now,
                                         delta_t=self.delta_t, latest_event=latest_event,
                                         fc_max_time=dt.timedelta(hours=forecast_hours), update_time=update_time,
                                         regrid_description=self.regrid_description,
                                         crop_description=self.crop_description, scale_factor=self.scale_factor,
                                         fill_value=self.fill_value, eps_member=eps_member, short=self.short)
        if fc_time is None:
            # no IconD2EPS data available
            return
        self.amend_time_data(time=fc_time, data=fc_data)
        self.eps_member = eps_member
        self.eps_member_source = 'IconD2EPS'
        self.icond2eps_loaded = True

    def collect_iconeu(self,
                       latest_event,
                       directory=None,
                       dir_time_descriptor=None,
                       forecast_hours=120):
        """
        Collect all IconEU data into the variable data and construct the according time vector. collect_radolanrw must
        be invoked before collecting prediction data.

        :param latest_event: datetime of the latest forecast start to be tried to load
        :type latest_event: datetime.datetime
        :param directory: directory with IconEU data
        :type directory: str, optional
        :param dir_time_descriptor: list of datetime.strftime time descriptors for an arbitrary number of additional
            time dependent folders with the data, e.g. ['%Y', '%Y%m%d'] for */yyyy/yyyymmdd/*; the final path is built
            to directory/dir_time_directory[0]/../dir_time_directory[n]/filename
        :type dir_time_descriptor: list, optional
        :param forecast_hours: try to load this number of hours at maximum
        :type forecast_hours: int, optional
        """
        if not self.radolanrw_loaded:
            raise Exception('you must first load RadolanRW data before loading forecast data')
        if self.iconeu_loaded:
            print('IconEU already loaded; nothing will be done')
            return

        update_time = IconEUProp().update_time
        update_time_in_hours = update_time.total_seconds()/3600
        if latest_event.hour % update_time_in_hours != 0:
            # latest event does not fit to update cycle replace the latest event to the next future one
            latest_event = latest_event + dt.timedelta(hours=(update_time_in_hours -
                                                              (latest_event.hour % update_time_in_hours)))

        print('collect IconEU data')

        fc_data, fc_time, _ = collect_fc(product='iconeu', directory=directory, dir_time_descriptor=dir_time_descriptor,
                                         time_now=self.time_now, delta_t=self.delta_t, latest_event=latest_event,
                                         fc_max_time=dt.timedelta(hours=forecast_hours), update_time=update_time,
                                         regrid_description=self.regrid_description,
                                         crop_description=self.crop_description, scale_factor=self.scale_factor,
                                         fill_value=self.fill_value, short=self.short)
        if fc_time is None:
            # no IconEU data available
            return
        self.amend_time_data(time=fc_time, data=fc_data)
        self.iconeu_loaded = True

    def collect_iconeueps(self,
                          latest_event,
                          directory=None,
                          dir_time_descriptor=None,
                          forecast_hours=120,
                          eps_member=range(40)):
        """
        Collect all IconEUEPS data for one ensemble member into the variable data and construct the according time
        vector. collect_radolanrw must be invoked before collecting prediction data.

        :param latest_event: datetime of the latest forecast start to be tried to load
        :type latest_event: datetime.datetime
        :param directory: directory with IconEUEPS data
        :type directory: str, optional
        :param dir_time_descriptor: list of datetime.strftime time descriptors for an arbitrary number of additional
            time dependent folders with the data, e.g. ['%Y', '%Y%m%d'] for */yyyy/yyyymmdd/*; the final path is built
            to directory/dir_time_directory[0]/../dir_time_directory[n]/filename
        :type dir_time_descriptor: list, optional
        :param forecast_hours: try to load this number of hours at maximum
        :type forecast_hours: int, optional
        :param eps_member: list of eps members to load (0 to 39)
        :type eps_member: list, optional
        """
        if not self.radolanrw_loaded:
            raise Exception('you must first load RadolanRW data before loading forecast data')
        if self.iconeueps_loaded:
            print('IconEUEPS already loaded; nothing will be done')
            return

        update_time = IconEUEPSProp().update_time
        update_time_in_hours = update_time.total_seconds() / 3600
        if latest_event.hour % update_time_in_hours != 0:
            # latest event does not fit to update cycle replace the latest event to the next future one
            latest_event = latest_event + dt.timedelta(hours=(update_time_in_hours -
                                                              (latest_event.hour % update_time_in_hours)))

        print('collect IconEUEPS data')

        fc_data, fc_time, _ = collect_fc(product='iconeueps', directory=directory,
                                         dir_time_descriptor=dir_time_descriptor, time_now=self.time_now,
                                         delta_t=self.delta_t, latest_event=latest_event,
                                         fc_max_time=dt.timedelta(hours=forecast_hours), update_time=update_time,
                                         regrid_description=self.regrid_description,
                                         crop_description=self.crop_description, scale_factor=self.scale_factor,
                                         fill_value=self.fill_value, eps_member=eps_member, short=self.short)
        if fc_time is None:
            # no IconEUEPS data available
            return
        self.amend_time_data(time=fc_time, data=fc_data)
        self.eps_member = eps_member
        self.eps_member_source = 'IconEUEPS'
        self.iconeueps_loaded = True

    def amend_time_data(self, time, data):
        """
        Amend the time vector in self.time from existing datetimes until latest_time. self.data is extended and filled
        accordingly. If (single or all) datetime already available, the existing forecast data is overwritten. This can
        be used to load a longer forecasting product and overwrite some time steps with a (potentially better)
        short-time forecast product. Time steps which are in self.time_protected (typically observed data) are never
        changed.

        :param time: the latest datetime to be included in self.time
        :type time: list
        :param data: list of data as GeoreferencedData instances
        :type data: list
        """
        if time[0] < self.time[0]:
            # no amending before available time
            raise Exception('try amending forecast time start before available time')
        if time[-1] > self.time[-1]:
            # amend time vector
            while self.time[-1] < time[-1]:
                self.time.append(self.time[-1] + self.delta_t)
                self.data.append(None)

        num_members = 1
        if isinstance(self.data[0], GeoReferencedData) and isinstance(data[0], GeoReferencedData):
            # one ensemble member in self.data and data -> no multiplication necessary
            multiplication_case = 0
        elif isinstance(self.data[0], GeoReferencedData) and isinstance(data[0], list):
            # one ensemble member in self.data and multiple ensemble members in data -> multiplicate self.data
            multiplication_case = 1
            num_members = len(data)
            self.data = [self.data]
            for member in range(1, num_members):
                self.data.append(copy.deepcopy(self.data[0]))
        elif isinstance(self.data[0], list) and isinstance(data[0], GeoReferencedData):
            # multiple ensemble members in self.data and one ensemble member in data -> multiply data
            multiplication_case = 2
            num_members = len(self.data)
        elif isinstance(self.data[0], list) and isinstance(data[0], list):
            # multiple ensemble members in self.data and data -> exception
            raise Exception('Collected data has already multiple ensemble members. Adding further ensemble members not possible')
        else:
            raise Exception('Something went wrong with the dimensionality of amending of collected data')

        for i in range(len(time)):
            if time[i] in self.time_protected:
                # no overwriting of protected time steps
                continue
            idx = self.time.index(time[i])
            if multiplication_case == 0:
                self.data[idx] = data[i]
            elif multiplication_case == 1:
                for member in range(num_members):
                    self.data[member][idx] = data[member][i]
            elif multiplication_case == 2:
                for member in range(num_members):
                    self.data[member][idx] = data[i]


    def export_netcdf(self, filename, data_format='f8', institution=None, scale_factor_nc=1, scale_undo=False,
                      data_kwargs=None):
        """
        Export the relevant content of this WeatherData class to a new netcdf file.

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

        if scale_undo and self.scale_factor == 1:
            # roll back scaling called, but not necessary
            scale_undo = False

        precipitation_names = []  # only used for separated variables
        if len(self.eps_member) > 1:
            # multiple ensemble members available
            separation = True
            for member in self.eps_member:
                precipitation_names.append(f'{self.nc_desc.var_data}_{member}')
            example_data = self.data[0]
        elif len(self.eps_member) == 1:
            # one ensemble member available
            separation = True
            precipitation_names.append(f'{self.nc_desc.var_data}_{self.eps_member[0]}')
            example_data = self.data
        else:
            # no ensemble product used
            separation = False
            example_data = self.data

        # create file and dimensions
        ncfile = nc.Dataset(filename, mode='w', format='NETCDF4')
        ncfile.createDimension(self.nc_desc.dim_lon, example_data[0].lon.data.shape[1])  # longitude axis
        ncfile.createDimension(self.nc_desc.dim_lat, example_data[0].lon.data.shape[0])  # latitude axis
        ncfile.createDimension(self.nc_desc.dim_time, None)  # unlimited time axis
        ncfile.title = 'weather data with'
        if self.radvorrq_loaded:
            ncfile.title = f'{ncfile.title} RadvorRQ data and'
        if self.radolanrv_loaded:
            ncfile.title = f'{ncfile.title} RadolanRV data and'
        if self.icond2_loaded:
            ncfile.title = f'{ncfile.title} IconD2 data and'
        if self.icond2eps_loaded:
            ncfile.title = f'{ncfile.title} IconD2EPS data and'
        ncfile.title = f'{ncfile.title} RadolanRW data'
        if institution is not None:
            ncfile.institution = institution
        ncfile.generator = 'generated with weatherDataHarmonizer'
        ncfile.source = 'data from DWD available on https://opendata.dwd.de/'
        ncfile.created = dt.datetime.now(tz=dt.timezone.utc).isoformat()

        # create variables in file
        lon = ncfile.createVariable(self.nc_desc.var_lon, 'f4', (self.nc_desc.dim_lat, self.nc_desc.dim_lon),
                                    fill_value=self.fill_value)
        lon.units = example_data[0].lon.data_description.units
        lon.long_name = example_data[0].lon.data_description.long_name
        lon.coordinate_system = example_data[0].lon.data_description.coordinate_system

        lat = ncfile.createVariable(self.nc_desc.var_lat, 'f4', (self.nc_desc.dim_lat, self.nc_desc.dim_lon),
                                    fill_value=self.fill_value)
        lat.units = example_data[0].lat.data_description.units
        lat.long_name = example_data[0].lat.data_description.long_name
        lat.coordinate_system = example_data[0].lat.data_description.coordinate_system

        time = ncfile.createVariable(self.nc_desc.var_time, 'f8', (self.nc_desc.dim_time,))
        time.units = f'hours since {self.time[0].isoformat()}'
        time.calendar = 'standard'

        precipitation_units = example_data[0].data_description.units.split(' ')
        if scale_undo:
            precip_multiplier = float(precipitation_units[0]) * self.scale_factor
            scale_factor_internal = 1
        else:
            precip_multiplier = precipitation_units[0]
            scale_factor_internal = self.scale_factor
        precipitation_units = f'{precip_multiplier} {" ".join(precipitation_units[1:])}'
        precipitation = None
        if separation:
            for member in range(len(self.eps_member)):
                vars()[precipitation_names[member]] = ncfile.createVariable(precipitation_names[member], data_format,
                                               (self.nc_desc.dim_time, self.nc_desc.dim_lat, self.nc_desc.dim_lon),
                                                fill_value=self.fill_value, **data_kwargs)
                vars()[precipitation_names[member]].units = precipitation_units
                vars()[precipitation_names[member]].long_name = 'harmonized precipitation data from different sources'
                vars()[precipitation_names[member]].scale_factor_internal = scale_factor_internal
                vars()[precipitation_names[member]].scale_factor = scale_factor_nc
                vars()[precipitation_names[member]].ensemble_member = f'{self.eps_member[member]} for source {self.eps_member_source}'
        else:
            precipitation = ncfile.createVariable(self.nc_desc.var_data, data_format,
                                               (self.nc_desc.dim_time, self.nc_desc.dim_lat, self.nc_desc.dim_lon),
                                                fill_value=self.fill_value, **data_kwargs)
            precipitation.units = precipitation_units
            precipitation.long_name = 'harmonized precipitation data from different sources'
            precipitation.scale_factor_internal = scale_factor_internal
            precipitation.scale_factor = scale_factor_nc

        data_source = ncfile.createVariable('data_source', 'i2', (self.nc_desc.dim_time,), fill_value=self.fill_value)
        data_source.description = 'data source indication per time step'
        data_source.ind_0 = 'RadolanRW'
        data_source.ind_1 = 'RadvorRQ'
        data_source.ind_2 = 'RadolanRV'
        data_source.ind_3 = 'IconD2'
        data_source.ind_4 = 'IconD2EPS'
        data_source.ind_5 = 'IconEU'
        data_source.ind_6 = 'IconEUEPS'

        time_data_source = ncfile.createVariable('time_data_source', 'f8', (self.nc_desc.dim_time,), fill_value=self.fill_value)
        time_data_source.description = 'reference time of data source per time step; for RadolanRW product only the first time step is given'
        time_data_source.units = f'hours since {self.time[0].isoformat()}'

        # write data to variables
        lon[:, :] = example_data[0].lon.data
        lat[:, :] = example_data[0].lat.data

        # data sources defined with indicator
        num_time = len(self.time)
        for t in range(num_time):
            time[t] = (self.time[t] - self.time[0]).total_seconds() / 3600
            if separation and len(self.eps_member) > 1:
                for member in range(len(self.eps_member)):
                    vars()[precipitation_names[member]][t, :, :] = \
                        st.gr_data_scaling(self.data[member][t], scale_undo, self.scale_factor, scale_factor_nc)
            elif separation and len(self.eps_member) == 1:
                vars()[precipitation_names[0]][t, :, :] = st.gr_data_scaling(self.data[t], scale_undo,
                                                                             self.scale_factor, scale_factor_nc)
            else:
                precipitation[t, :, :] = st.gr_data_scaling(self.data[t], scale_undo, self.scale_factor,
                                                            scale_factor_nc)

            product = example_data[t].data_description.long_name.split(' ')[0]
            data_source[t] = self.data_source_dict[product.lower()]
            time_data_source[t] = (dt.datetime.fromisoformat(example_data[t].data_description.time_note.split(' ')[-1])
                                   - self.time[0]).total_seconds() / 3600

        ncfile.close()

    def export_netcdf_append(self, filename):
        """
        Append or replace the relevant content of this WeatherData class to an existing netcdf file. It is possible to
        broadcast one (here) to many (netcdf), but it is not allowed to store many (here) to one (netcdf). In the latter
        case the standard variable name in the netcdf would have to be removed. That is not supported in netcdf4
        standard.

        Newer times (here) can be amended to netcdf, extending older time is not supported.

        New precipitation data (here) will overwrite existing data in netcdf with the same datetime. If multiple
        ensemble members exist in netcdf, newer data without ensemble members will overwrite all ensemble members at the
        appropriate times.

        Several checks are conducted to ensure the compatibility of the netcdf with the data to be amended before the
        netcdf file is changed.

        :param filename: filename of an existing netcdf file with the same general properties
        :type filename: str
        """
        print('append to existing netcdf')

        if not os.path.exists(filename):
            raise Exception(f'file {filename} does not exist')

        precipitation_names = []  # only used for separated variables
        if len(self.eps_member) > 1:
            # multiple ensemble members available
            separation_internal = True
            for member in self.eps_member:
                precipitation_names.append(f'{self.nc_desc.var_data}_{member}')
            example_data = self.data[0]
        elif len(self.eps_member) == 1:
            # one ensemble member available
            separation_internal = True
            precipitation_names.append(f'{self.nc_desc.var_data}_{self.eps_member[0]}')
            example_data = self.data
        else:
            # no ensemble product used
            separation_internal = False
            precipitation_names.append(self.nc_desc.var_data)
            example_data = self.data

        ncfile = nc.Dataset(filename, mode='a')

        # check variable names
        required_variables = [self.nc_desc.var_lon, self.nc_desc.var_lat, self.nc_desc.var_time, 'data_source', 'time_data_source']
        for required_variable in required_variables:
            if required_variable not in ncfile.variables:
                raise Exception(f'variable {required_variable} not in {filename}')

        # check shape
        if ncfile.variables[self.nc_desc.var_lon].shape != example_data[0].data.shape:
            raise Exception(f'data shape differs from {filename} content')

        # check corner coordinates
        if not np.isclose(example_data[0].lon.data[0, 0], float(ncfile.variables[self.nc_desc.var_lon][0, 0])) or \
            not np.isclose(example_data[0].lon.data[0, -1], float(ncfile.variables[self.nc_desc.var_lon][0, -1])) or \
            not np.isclose(example_data[0].lon.data[-1, 0], float(ncfile.variables[self.nc_desc.var_lon][-1, 0])) or \
            not np.isclose(example_data[0].lon.data[-1, -1], float(ncfile.variables[self.nc_desc.var_lon][-1, -1])) or \
            not np.isclose(example_data[0].lat.data[0, 0], float(ncfile.variables[self.nc_desc.var_lat][0, 0])) or \
            not np.isclose(example_data[0].lat.data[0, -1], float(ncfile.variables[self.nc_desc.var_lat][0, -1])) or \
            not np.isclose(example_data[0].lat.data[-1, 0], float(ncfile.variables[self.nc_desc.var_lat][-1, 0])) or \
            not np.isclose(example_data[0].lat.data[-1, -1], float(ncfile.variables[self.nc_desc.var_lat][-1, -1])):
            raise Exception(f'corner coordinates different from corners in {filename}')

        # check eps members in netcdf
        precipitation_names_nc = [key for key in ncfile.variables if key.startswith(self.nc_desc.var_data)]
        if len(precipitation_names_nc) == 0:
            # no matching variable name in netcdf
            raise Exception(f'no {self.nc_desc.var_data} like variable in {filename}')
        elif self.nc_desc.var_data in precipitation_names_nc:
            # no ensemble product used
            separation_nc = False
        else:
            # ensemble members available
            separation_nc = True
        if separation_nc and separation_internal:
            if len(precipitation_names_nc) != len(precipitation_names):
                raise Exception(f'internal number of ensemble members does not match the number in {filename}')
            for precipitation_name in precipitation_names:
                if precipitation_name not in precipitation_names_nc:
                    raise Exception(f'internal variable {precipitation_name} not in {filename}')
        elif not separation_nc and separation_internal:
            # no ensemble member in netcdf but ensemble members in internal version is not allowed, because the original
            # variable would have to be deleted -> that is not included in netcdf functionality
            raise Exception(f'{filename} has no ensemble members, but internally ensemble members are requested - not allowed')

        # check some attributes - roll back short mode if scale_factor_internal different from self.scale_factor
        scale_undo = False
        if ncfile.variables[precipitation_names_nc[0]].scale_factor_internal == 1 and self.scale_factor != 1:
            scale_undo = True
        if ncfile.variables[precipitation_names_nc[0]].scale_factor_internal != 1 and \
                not np.isclose(ncfile.variables[precipitation_names_nc[0]].scale_factor_internal, self.scale_factor):
            raise Exception(f'internal scale factor {self.scale_factor} does not match internal scale factor in {filename}')
        if self.fill_value is not np.nan:
            if ncfile.variables[precipitation_names_nc[0]]._FillValue != self.fill_value:
                raise Exception(f'internal fill value {self.fill_value} does not match fill value in {filename}')
        else:
            if ncfile.variables[precipitation_names_nc[0]]._FillValue is not np.nan:
                raise Exception(f'internal fill value {self.fill_value} does not match fill value in {filename}')

        # check time and time resolution
        time_nc = ncfile.variables['time']
        time_nc_units = time_nc.units.split(' ')
        base_time = dt.datetime.fromisoformat(time_nc_units[-1])
        time_nc_values = []
        if time_nc_units[0] == 'hours':
            for time_act in time_nc:
                time_nc_values.append(base_time + dt.timedelta(hours=float(time_act.data)))
        elif time_nc_units[0] == 'minutes':
            for time_act in time_nc:
                time_nc_values.append(base_time + dt.timedelta(minutes=float(time_act.data)))
        else:
            raise Exception(f'base time description {time_nc_units[0]} not supported')
        if time_nc_values[0] > self.time[0]:
            raise Exception(f'start time in {filename} is after start time to be appended - only appending of later time steps supported')
        delta_time_nc_values = np.array(time_nc_values[1:]) - np.array(time_nc_values[0:-1])
        if not np.all(delta_time_nc_values == delta_time_nc_values[0]):
            raise Exception(f'time resolution in {filename} not constant')
        if delta_time_nc_values[0] != self.delta_t:
            raise Exception(f'time resolution in {filename} does not match internal resolution')

        # amend internal time to time in netcdf, extend it necessary
        ct = len(time_nc) - 1
        while self.time[-1] > time_nc_values[-1]:
            ct = ct + 1
            time_nc_values.append(time_nc_values[ct - 1] + self.delta_t)
            if time_nc_units[0] == 'hours':
                time_nc[ct] = time_nc[ct - 1] + self.delta_t.total_seconds() / 3600
            else:
                time_nc[ct] = time_nc[ct - 1] + self.delta_t.total_seconds() / 60

        # add data to netcdf file
        scale_factor_nc = ncfile.variables[precipitation_names_nc[0]].scale_factor
        for time_ct in range(len(self.time)):
            idx = time_nc_values.index(self.time[time_ct])

            # add precipitation
            if not separation_internal and not separation_nc:
                ncfile.variables[self.nc_desc.var_data][idx, :, :] = \
                    st.gr_data_scaling(self.data[time_ct], scale_undo, self.scale_factor, scale_factor_nc)
            elif not separation_internal and separation_nc:
                for precipitation_name_nc in precipitation_names_nc:
                    ncfile.variables[precipitation_name_nc][idx, :, :] = \
                        st.gr_data_scaling(self.data[time_ct], scale_undo, self.scale_factor, scale_factor_nc)
            else:
                for precipitation_name_nc in precipitation_names_nc:
                    idx_eps_member = self.eps_member.index(int(precipitation_name_nc.split('_')[-1]))
                    ncfile.variables[precipitation_name_nc][idx, :, :] = \
                        st.gr_data_scaling(self.data[idx_eps_member][time_ct], scale_undo, self.scale_factor,
                                           scale_factor_nc)

            # add data source
            product = example_data[time_ct].data_description.long_name.split(' ')[0]
            ncfile.variables['data_source'][idx] = self.data_source_dict[product.lower()]
            reference_time_act = dt.datetime.fromisoformat(example_data[time_ct].data_description.time_note.split(' ')[-1])
            if time_nc_units[0] == 'hours':
                ncfile.variables['time_data_source'][idx] = (reference_time_act - base_time).total_seconds() / 3600
            else:
                ncfile.variables['time_data_source'][idx] = (reference_time_act - base_time).total_seconds() / 60

        ncfile.last_modified = dt.datetime.now(tz=dt.timezone.utc).isoformat()

        ncfile.close()


def collect_fc(product, directory, time_now, delta_t, latest_event, update_time, dir_time_descriptor=None,
               fc_max_time=dt.timedelta(minutes=0), regrid_description=None, crop_description=None, scale_factor=1,
               fill_value=np.nan, eps_member=range(20), short=None):
    """
    Collection of available forecast data including loading, aggregation/disaggregation as needed, and a consistent time
    vector.

    :param product: product name; currently supported: radvorrq, radolanrv, icond2, icond2eps, iconeu, iconeueps
    :type product: str
    :param directory: directory with the forecast data
    :type directory: str
    :param time_now: current time step (ergo start of forecast time)
    :type time_now: datetime.datetime
    :param delta_t: required temporal resolution
    :type: delta_t: datetime.timedelta
    :param latest_event: time of the latest forecast start event that shall potentially be included; must be a time
        that is potentially available (e.g. 15:00 for icond2, or 10:15 for radvorrq)
    :type latest_event: datetime.datetime
    :param fc_max_time: potential maximum forecast time for the specific forecast product
    :type fc_max_time: datetime.timedelta
    :param update_time: time between updates of the forecast product
    :type update_time: datetime.timedelta
    :param dir_time_descriptor: list of datetime.strftime time descriptors for an arbitrary number of additional
        time dependent folders with the data, e.g. ['%Y', '%Y%m%d'] for */yyyy/yyyymmdd/*; the final path is built
        to directory/dir_time_directory[0]/../dir_time_directory[n]/filename
    :type dir_time_descriptor: list, optional
    :param regrid_description: dictionary with product key and a description as VariableDescription.RegridDescription
        instance
    :type regrid_description: dict
    :param crop_description: crop description
    :type crop_description: VariableDescription.CropDescription
    :param scale_factor: the final data has to be multiplied with this value
    :type scale_factor: float, optional
    :param fill_value: missing data is filled with that value
    :type fill_value: float, optional
    :param eps_member: list of eps members to load (e.g. in IconD2EPS max possible: 0 to 19 - default,
        in IconEUEPS: 0 to 39)
    :type eps_member: list, optional
    :param short: if int16 or int32, the large data variables are cast to numpy.int16/int32 to minimize memory
        usage; the user must pay attention to the scale_factor to include the necessary precision
    :type short: str, optional
    :return: (i) list of GeoReferencedData with the most recent forecast value for each time step, (ii) list of time
        steps that belong to the data
    :rtype: list, list
    """
    # time vector for reading data
    start_potential = [latest_event]
    while (start_potential[-1] + fc_max_time - update_time) >= time_now - delta_t:
        start_potential.append(start_potential[-1] - update_time)

    # read data
    data_list = []
    time_now_changed = None
    ct = -1
    for start in start_potential:
        print(f'try loading of {product} data at {start}')
        ct = ct + 1
        try:
            if product == 'radolanrw':
                data = RadolanRW()
                data.read_file(start_datetime=start, directory=directory, dir_time_descriptor=dir_time_descriptor,
                               scale_factor=scale_factor, fill_value=fill_value, short=short)
            elif product == 'radvorrq':
                data = RadvorRQ()
                data.read_file(start_datetime=start, directory=directory, dir_time_descriptor=dir_time_descriptor,
                               scale_factor=scale_factor, fill_value=fill_value, short=short)
            elif product == 'radolanrv':
                data = RadolanRV()
                data.read_file(start_datetime=start, directory=directory, dir_time_descriptor=dir_time_descriptor,
                               scale_factor=scale_factor, fill_value=fill_value, short=short)
            elif product == 'icond2':
                # enable netcdf import, if available
                nc_filename = get_potential_netcdf_filename(product_class=IconD2Prop(), time_step=start,
                                                            directory=directory,
                                                            dir_time_descriptor=dir_time_descriptor)
                try:
                    data = read_netcdf(filename=nc_filename, scale_factor_result=scale_factor,
                                       fill_value_result=fill_value, short_result=short)
                except IOError:
                    data = IconD2()
                    data.read_file(start_datetime=start, directory=directory, dir_time_descriptor=dir_time_descriptor,
                                   forecast_hours=int(fc_max_time.total_seconds() / 3600), scale_factor=scale_factor,
                                   fill_value=fill_value, short=short)
            elif product == 'icond2eps':
                # enable netcdf import, if available
                nc_filename = get_potential_netcdf_filename(product_class=IconD2EPSProp(), time_step=start,
                                                            directory=directory,
                                                            dir_time_descriptor=dir_time_descriptor)
                try:
                    data = read_netcdf(filename=nc_filename, scale_factor_result=scale_factor,
                                       fill_value_result=fill_value, short_result=short)
                except IOError:
                    data = IconD2EPS()
                    data.read_file(start_datetime=start, directory=directory, dir_time_descriptor=dir_time_descriptor,
                                   forecast_hours=int(fc_max_time.total_seconds() / 3600), eps_member=eps_member,
                                   scale_factor=scale_factor, fill_value=fill_value, short=short)
            elif product == 'iconeu':
                # enable netcdf import, if available
                nc_filename = get_potential_netcdf_filename(product_class=IconEUProp(), time_step=start,
                                                            directory=directory,
                                                            dir_time_descriptor=dir_time_descriptor)
                try:
                    data = read_netcdf(filename=nc_filename, scale_factor_result=scale_factor,
                                       fill_value_result=fill_value, short_result=short)
                except IOError:
                    data = IconEU()
                    data.read_file(start_datetime=start, directory=directory, dir_time_descriptor=dir_time_descriptor,
                                   forecast_hours=int(fc_max_time.total_seconds() / 3600), scale_factor=scale_factor,
                                   fill_value=fill_value, short=short)
            elif product == 'iconeueps':
                # enable netcdf import, if available
                nc_filename = get_potential_netcdf_filename(product_class=IconEUEPSProp(), time_step=start,
                                                            directory=directory,
                                                            dir_time_descriptor=dir_time_descriptor)
                try:
                    data = read_netcdf(filename=nc_filename, scale_factor_result=scale_factor,
                                       fill_value_result=fill_value, short_result=short)
                except IOError:
                    data = IconEUEPS()
                    data.read_file(start_datetime=start, directory=directory, dir_time_descriptor=dir_time_descriptor,
                                   forecast_hours=int(fc_max_time.total_seconds() / 3600), eps_member=eps_member,
                                   scale_factor=scale_factor, fill_value=fill_value, short=short)
            else:
                raise Exception(f'product {product} not supported')

            if regrid_description is not None:
                data.regrid(regrid_description=regrid_description[product])
            if crop_description is not None:
                data.crop(crop_description)
            data_list.append(data)
            if not data_list_time_interpretation(data_list, time_now, fc_max_time, update_time):
                # no further loading of data necessary
                print(f'stop importing here because forecast times already fully cover the forecast period from '
                      f'{time_now} on')
                break
        except ForecastFileNotAvailable:
            print('could not be loaded')
            continue
        except ForecastFileFormatError:
            print('could not be loaded due to format issues')
            continue
        except RadarFileNotAvailable:
            print('could not be loaded')
            if ct < len(start_potential) - 1 and len(data_list) == 0:
                # change time_now due to missing first radar file
                time_now_changed = start_potential[ct + 1]
                print(f'time_now changed to {time_now_changed} due to missing radolanrw data')
            elif ct == len(start_potential) - 1:
                # first potential start for radolanrw data typically not necessary -> do not raise exception
                pass
            else:
                raise Exception('no radolanrw file available')
            continue

    if len(data_list) == 0:
        print(f'no data was load for product {product}')
        return None, None

    if product != 'radolanrw':
        # forecast relevant times
        if data_list[0].forecast_value.data_description.units == 'minutes':
            fc_delta_t = dt.timedelta(minutes=int(data_list[0].forecast_value.data[1] - data_list[0].forecast_value.data[0]))
        elif data_list[0].forecast_value.data_description.units == 'hours':
            fc_delta_t = dt.timedelta(hours=int(data_list[0].forecast_value.data[1] - data_list[0].forecast_value.data[0]))
        else:
            raise Exception(f'time unit {data_list[0].forecast_value.data_description.units} not known')
    else:
        # in case or radolanrw data use update_time as forecast_delta_t and concatenate all data from data_list in
        # chronological order (reverse from reading state) and shift datetime by 10 min for full hour matching in all
        # cases but 5 min temporal resolution
        if delta_t == dt.timedelta(minutes=5):
            radrw_shift = dt.timedelta(minutes=0)
        else:
            radrw_shift = dt.timedelta(minutes=10)
        fc_delta_t = update_time
        data_list[-1].time_value.data = [data_list[-1].time_value.data + radrw_shift.total_seconds() / 3600]
        if short:
            data_list_data = np.empty([data_list[-1].gr_data.data.shape[0], data_list[-1].gr_data.data.shape[1],
                                       len(data_list)], dtype=np.short)
        else:
            data_list_data = np.empty([data_list[-1].gr_data.data.shape[0], data_list[-1].gr_data.data.shape[1],
                                       len(data_list)])
        data_list_data[:, :, 0] = data_list[-1].gr_data.data
        tmp_time_start = data_list[-1].time_value.data_description.units.split(' ')
        ct = 0
        for i in range(len(data_list)-2, -1, -1):
            ct = ct + 1
            tmp_time_act = data_list[i].time_value.data_description.units.split(' ')
            if tmp_time_act[0] != 'hours':
                raise Exception(f'found {tmp_time_act[0]} as unit of time description - only hours are supported')
            time_act = ((dt.datetime.fromisoformat(tmp_time_act[2]) + radrw_shift + dt.timedelta(
                hours=data_list[i].time_value.data)) - (
                                  dt.datetime.fromisoformat(tmp_time_start[2]))).total_seconds() / 3600
            data_list[-1].time_value.data.append(time_act)
            data_list_data[:, :, ct] = data_list[i].gr_data.data
        data_list[-1].time_value.data = np.array(data_list[-1].time_value.data)
        data_list[-1].gr_data.data = data_list_data
        data_list = [data_list[-1]]

    # get time list for external usage
    num_time_external = int((latest_event + fc_max_time - time_now) / delta_t)
    if product != 'radolanrw':
        time_external = [time_now + delta_t + x * delta_t for x in range(num_time_external)]
        pivot_time = time_now
    else:
        # also add the very first time step for observed time steps
        time_external = [time_now + x * delta_t for x in range(num_time_external + 1)]
        pivot_time = time_external[-1]

    # assign internal data to external time
    if isinstance(data_list[0].gr_data, list):
        # in case of multiple ensemble members process each member separately (e.g. in IconD2EPS data)
        data_collected = []
        for member in range(len(data_list[0].gr_data)):
            data_list_act = copy.deepcopy(data_list)
            for read_files in range(len(data_list_act)):
                data_list_act[read_files].gr_data = data_list_act[read_files].gr_data[member]
            data_list_transient, time_list_transient = create_transient_data_matrix(data_list_act, fc_delta_t, delta_t,
                                                                                    pivot_time, short=short)

            data_collected_part = []
            for time_external_act in time_external:
                for fc_act in range(len(time_list_transient)):
                    if time_external_act in time_list_transient[fc_act]:
                        # take the most recent data available
                        idx = time_list_transient[fc_act].index(time_external_act)
                        data_collected_act = GeoReferencedData(lon=data_list_transient[fc_act].lon,
                                                               lat=data_list_transient[fc_act].lat,
                                                               data=data_list_transient[fc_act].data[:, :, idx],
                                                               data_description=data_list_transient[
                                                                   fc_act].data_description)
                        data_collected_part.append(data_collected_act)
                        break
            data_collected.append(data_collected_part)
        time_external = time_external[0:len(data_collected[0])]
    else:
        data_list_transient, time_list_transient = create_transient_data_matrix(data_list, fc_delta_t, delta_t,
                                                                                pivot_time, short=short)

        data_collected = []
        for time_external_act in time_external:
            for fc_act in range(len(time_list_transient)):
                if time_external_act in time_list_transient[fc_act]:
                    # take the most recent data available
                    idx = time_list_transient[fc_act].index(time_external_act)
                    data_collected_act = GeoReferencedData(lon=data_list_transient[fc_act].lon,
                                                           lat=data_list_transient[fc_act].lat,
                                                           data=data_list_transient[fc_act].data[:, :, idx],
                                                           data_description=data_list_transient[fc_act].data_description)
                    data_collected.append(data_collected_act)
                    break
        time_external = time_external[0:len(data_collected)]

    return data_collected, time_external, time_now_changed


def create_transient_data_matrix(data_list, fc_delta_t, delta_t, time_now, short=None):
    """
    Creation of a transient data matrix with target temporal resolution. Aggregation/disaggregation is done as needed.
    Currently only works for two-dimensional forecast data and a third dimension with the forecast time.

    :param data_list: list of forecast data
    :type data_list: list
    :param fc_delta_t: temporal resolution of forecast data
    :type fc_delta_t: datetime.timedelta
    :param delta_t: temporal resolution of target data
    :type delta_t: datetime.timedelta
    :param time_now: last time step of observations, afterwards forecast starts; only used for aggregation
    :type time_now: datetime.datetime
    :param short: if int16 or int32, the large data variables are cast to numpy.int16/int32 to minimize memory
        usage; the user must pay attention to the scale_factor to include the necessary precision
    :type short: str, optional
    :return: transient list of forecast data and time vectors
    :rtype: list, list
    """
    if isinstance(data_list[0], RadolanRW):
        # RadolanRW has no forecast dimension -> special handling necessary
        radolanrw = True
    else:
        radolanrw = False

    divider = fc_delta_t / delta_t
    data_list_transient = []
    time_list_transient = []
    if divider > 1:
        # disaggregation necessary -> create revamped data_list with target resolution
        for data in data_list:
            if not radolanrw:
                if data.forecast_value.data[0] == 0:
                    # delete forecast time = 0, this is not a forecast (and in the case of IconD2 always zero)
                    data.gr_data.data = data.gr_data.data[:, :, 1:]
                    data.forecast_value.data = data.forecast_value.data[1:]
                num_steps_orig = data.forecast_value.data.size

                # replace forecast time vector with disaggregated times
                delta_t_val = get_delta_t_val(data, fc_delta_t, divider)
                data.forecast_value.data = np.arange(data.forecast_value.data[0] - (divider - 1) * delta_t_val,
                                                     data.forecast_value.data[-1] + delta_t_val, delta_t_val)
                num_steps_disag = data.forecast_value.data.size
            else:
                # skip forecast specific parts in radolanrw
                num_steps_orig = data.time_value.data.size

                # replace time vector with disaggregated times
                delta_t_val = get_delta_t_val(data, fc_delta_t, divider)
                data.time_value.data = np.arange(data.time_value.data[0] - (divider - 1) * delta_t_val,
                                                 data.time_value.data[-1] + delta_t_val, delta_t_val)
                num_steps_disag = data.time_value.data.size


            # replace forecast data matrix with disaggregated matrices
            if short:
                tmp1 = np.empty([data.gr_data.data.shape[0], data.gr_data.data.shape[1], num_steps_disag], dtype=np.short)
            else:
                tmp1 = np.empty([data.gr_data.data.shape[0], data.gr_data.data.shape[1], num_steps_disag])
            ct = -1
            for time_step in range(num_steps_orig):
                if short:
                    tmp2 = np.short(np.around(data.gr_data.data[:, :, time_step] / divider))
                else:
                    tmp2 = data.gr_data.data[:, :, time_step] / divider
                if not np.isnan(data.gr_data.data_description.fill_value):
                    # due to divider cells with fill_value changed -> restore original fill_value
                    tmp2[data.gr_data.data[:, :, time_step] == data.gr_data.data_description.fill_value] = data.gr_data.data_description.fill_value
                for sub_step in range(int(divider)):
                    ct = ct + 1
                    tmp1[:, :, ct] = tmp2
            data.gr_data.data = tmp1

    elif divider < 1:
        # aggregation necessary -> create revamped data_list with target resolution
        data_remove = []
        for data in data_list:
            if not radolanrw:
                if data.forecast_value.data[0] == 0:
                    # delete forecast time = 0, this is not a forecast (and in the case of IconD2 always zero)
                    data.gr_data.data = data.gr_data.data[:, :, 1:]
                    data.forecast_value.data = data.forecast_value.data[1:]

                # replace forecast time vector with aggregated times
                time_value_data_datetime = get_timestamps(data)  # timestamps for comparison with time_now
                idx_time_value_data = get_time_value_data(time_value_data_datetime, time_now, delta_t)
                time_value_data = data.forecast_value.data[idx_time_value_data]
                num_steps_aggr = time_value_data.size
                if num_steps_aggr == 0:
                    # no aggregated steps usable for time_now + multiples of delta_t
                    data_remove.append(data)
                    continue
            else:
                # skip forecast specific parts in radolanrw
                # replace time vector with aggregated times
                time_value_data_datetime = get_timestamps(data)  # timestamps for comparison with time_now
                delta_t_val = get_delta_t_val(data, fc_delta_t, divider)
                # radolanrw data must be aggregated backward to ensure the availability of time_now
                idx_time_now = time_value_data_datetime.index(time_now)
                time_value_data = np.flip(np.arange(data.time_value.data[idx_time_now], data.time_value.data[0],
                                                    -delta_t_val))
                num_steps_aggr = time_value_data.size

            # replace forecast data matrix with aggregated matrices
            if short:
                tmp1 = np.empty([data.gr_data.data.shape[0], data.gr_data.data.shape[1], num_steps_aggr], dtype=np.short)
            else:
                tmp1 = np.empty([data.gr_data.data.shape[0], data.gr_data.data.shape[1], num_steps_aggr])
            for time_step in range(num_steps_aggr):
                tmp2 = np.zeros([data.gr_data.data.shape[0], data.gr_data.data.shape[1]])  # temporary aggregated matrix
                if not radolanrw:
                    idx_time_step = np.where(np.isclose(data.forecast_value.data, time_value_data[time_step]))[0][0]
                else:
                    idx_time_step = np.where(np.isclose(data.time_value.data, time_value_data[time_step]))[0][0]
                idx_missing = np.zeros([data.gr_data.data.shape[0], data.gr_data.data.shape[1]], dtype=bool)
                for idx in range(int(idx_time_step - 1 / divider + 1), int(idx_time_step + 1)):
                    tmp2 = tmp2 + data.gr_data.data[:, :, idx]
                    # if any missing value in aggregated time steps, the target time step also becomes missing
                    idx_missing[data.gr_data.data[:, :, idx] == data.gr_data.data_description.fill_value] = True
                tmp2[idx_missing] = data.gr_data.data_description.fill_value
                if short:
                    tmp1[:, :, time_step] = np.short(np.around(tmp2))
                else:
                    tmp1[:, :, time_step] = tmp2
            if not radolanrw:
                data.forecast_value.data = time_value_data
            else:
                data.time_value.data = time_value_data
            data.gr_data.data = tmp1

        if len(data_remove) > 0:
            # remove data that has no more elements due to aggregation
            for data in data_remove:
                data_list.remove(data)

    # dt matches (after aggregation or disaggregation if necessary)
    for data in data_list:
        start = dt.datetime.fromisoformat(data.time_value.data_description.units.split(' ')[-1])
        data.gr_data.data_description.units = f'{1 / data.gr_data.data_description.scale_factor} ' \
                                              f'mm/{delta_t.total_seconds() / 60} min'  # insert correct unit
        if not radolanrw:
            if data.forecast_value.data_description.units == 'minutes':
                event_time = [start + dt.timedelta(minutes=int(x)) for x in data.forecast_value.data]
            elif data.forecast_value.data_description.units == 'hours':
                event_time = [start + dt.timedelta(hours=int(x)) for x in data.forecast_value.data]
            else:
                raise Exception(f'time unit {data.forecast_value.data_description.units} not known')

            if data.forecast_value.data[0] == 0:
                # delete forecast time = 0, this is not a forecast (and in the case of IconD2 always zero)
                data.gr_data.data = data.gr_data.data[:, :, 1:]
                data_list_transient.append(data.gr_data)
                time_list_transient.append(event_time[1:])
            else:
                data_list_transient.append(data.gr_data)
                time_list_transient.append(event_time)
        else:
            # replace forecast specific parts for general time in radolanrw
            if data.time_value.data_description.units.split(' ')[0] == 'minutes':
                event_time = [start + dt.timedelta(minutes=x) for x in data.time_value.data]
            elif data.time_value.data_description.units.split(' ')[0] == 'hours':
                event_time = [start + dt.timedelta(hours=x) for x in data.time_value.data]
            else:
                raise Exception(f'time unit {data.time_value.data_description.units.split(" ")[0]} not known')
            data_list_transient.append(data.gr_data)
            time_list_transient.append(event_time)

    return data_list_transient, time_list_transient


def data_list_time_interpretation(data_list, time_now, fc_max_time, update_time):
    """
    Interpretation of time aspects of data in a list regarding full coverage until time_now and the potential maximum
    forecast time. If the beginning of the first forecast data is equal to or earlier than last observation, further
    data loading is advised. Further data loading is also advised if the newest forecast time plus its maximum forecast
    time is older than a potential forecast with its maximum forecast time one update step earlier than already loaded.
    The second case appears if, e.g. the actual (newest) forecast is not completely provided by DWD, yet.

    Loading observed data is always advised.

    :param data_list: list of forecast data
    :type data_list: list
    :param time_now: actual time when forecast starts
    :type time_now: datetime.datetime
    :param fc_max_time: potential maximum forecast time for the specific forecast product
    :type fc_max_time: datetime.timedelta
    :param update_time: time between updates of the forecast product
    :type update_time: datetime.timedelta
    :return: logical value whether further loading data is advised (True) or no not necessary (False)
    :rtype: bool
    """
    forecast_start = []
    forecast_end = []
    if hasattr(data_list[0], 'forecast_value') and data_list[0].forecast_value is not None:
        # check start of forecast time
        for data in data_list:
            time_start = dt.datetime.fromisoformat(data.time_value.data_description.units.split(' ')[-1])
            if data.forecast_value.data_description.units == 'minutes':
                forecast_start.append(time_start + dt.timedelta(minutes=float(data.forecast_value.data[0])))
                forecast_end.append(time_start + dt.timedelta(minutes=float(data.forecast_value.data[-1])))
            elif data.forecast_value.data_description.units == 'hours':
                forecast_start.append(time_start + dt.timedelta(hours=float(data.forecast_value.data[0])))
                forecast_end.append(time_start + dt.timedelta(hours=float(data.forecast_value.data[-1])))
            else:
                raise Exception(f'time unit {data.forecast_value.data_description.units} not known')

        # check potential end of forecast time
        forecast_end_latest = np.max(forecast_end)
        forecast_end_one_older = forecast_start[-1] - update_time + fc_max_time

        if any([forecast_start[x] <= time_now for x in range(len(forecast_start))]) and \
                forecast_end_one_older < forecast_end_latest:
            # no further loading seems necessary
            return False
        else:
            # further loading is advised
            return True
    else:
        # if no forecast_value is available (ergo observing data) further loading is advised
        return True


def get_timestamps(data):
    """
    Get datetime timestamps in time_data as a combination of time unit description and float value as distance.

    :param data: object with data; if attribute forecast_value exists this is taken, otherwise take attribute
        time_value.
    :type data: MetEntities.MetEntities
    :return: timestamps of the time object
    :rtype: list
    """
    tmp = data.time_value.data_description.units.split(' ')
    start_datetime = dt.datetime.fromisoformat(tmp[-1])

    if hasattr(data, 'forecast_value') and data.forecast_value is not None:
        if data.forecast_value.data_description.units == 'minutes':
            timestamps = [start_datetime + dt.timedelta(minutes=float(x)) for x in data.forecast_value.data]
        elif data.forecast_value.data_description.units == 'hours':
            timestamps = [start_datetime + dt.timedelta(hours=x) for x in data.forecast_value.data]
        else:
            raise Exception(f'time unit {data.forecast_value.data_description.units} not known')
    else:
        if tmp[0] == 'minutes':
            timestamps = [start_datetime + dt.timedelta(minutes=x) for x in data.time_value.data]
        elif tmp[0] == 'hours':
            timestamps = [start_datetime + dt.timedelta(hours=x) for x in data.time_value.data]
        else:
            raise Exception(f'time unit {data.time_value.data_description.units} not known')

    return timestamps


def get_time_value_data(time_value_data_datetime, time_now, delta_t):
    """
    Find indexes in timestamps which are met from time_now with a multiple delta_t distance, e.g. in a list 13:30,
    14:00, 14:30, 15:00, time_now = 13:00, delta_t = 1 h, find indexes 1 and 3. It can be used for aggregation.

    :param time_value_data_datetime: timestamps which shall be scanned through
    :type time_value_data_datetime: list
    :param time_now: time from which starts the multiple delta_t distances
    :type time_now: datetime.datetime
    :param delta_t: target delta t
    :type delta_t: datetime.timedelta
    :return: Indexes compliant to target resolution
    :rtype: list
    """
    idx_time_value = []
    delta_t_internal = time_value_data_datetime[1] - time_value_data_datetime[0]
    min_idx = int(delta_t / delta_t_internal) - 1  # minimum index to deliver right aggregations
    for idx in range(min_idx, len(time_value_data_datetime)):
        max_multiplier = int(np.ceil((time_value_data_datetime[-1] - time_now)/delta_t))
        time_list = [time_now + x * delta_t for x in range(max_multiplier + 1)]
        if time_value_data_datetime[idx] in time_list:
            idx_time_value.append(idx)

    return idx_time_value


def get_delta_t_val(data, fc_delta_t, divider):
    """
    Get delta_t_val, the delta_t of the values within data object.

    :param data: entity with meteorological data
    :type data: MetEntities.MetEntities
    :param fc_delta_t: delta_t of the data
    :type fc_delta_t: datetime.timedelta
    :param divider: division of the original data, > 1 means a disaggregation, < 1 an aggregation
    :type divider: float
    :return: the new delta_t for resulting values
    :rtype: float
    """
    if hasattr(data, 'forecast_value') and data.forecast_value is not None:
        if data.forecast_value.data_description.units == 'minutes':
            delta_t_val = fc_delta_t.total_seconds() / 60 / divider
        elif data.forecast_value.data_description.units == 'hours':
            delta_t_val = fc_delta_t.total_seconds() / 3600 / divider
        else:
            raise Exception(f'time unit {data.forecast_value.data_description.units} not known')
    else:
        if data.time_value.data_description.units.split(' ')[0] == 'minutes':
            delta_t_val = fc_delta_t.total_seconds() / 60 / divider
        elif data.time_value.data_description.units.split(' ')[0] == 'hours':
            delta_t_val = fc_delta_t.total_seconds() / 3600 / divider
        else:
            raise Exception(f'time unit {data.time_value.data_description.units.split(" ")[0]} not known')

    return delta_t_val


def get_potential_netcdf_filename(product_class, time_step, directory, dir_time_descriptor):
    """
    Build path string with a potentially existing netcdf containing the requested data.

    :param product_class: product object
    :type product_class: aux_tools.Product.Product
    :param time_step: time object of requested data
    :type time_step: datetime.datetime
    :param directory: directory with all IconD2 files from start_datetime
    :type directory: str, optional
    :param dir_time_descriptor: list of datetime.strftime time descriptors for an arbitrary number of additional
        time dependent folders with the data, e.g. ['%Y', '%Y%m%d'] for */yyyy/yyyymmdd/*; the final path is built
        to directory/dir_time_directory[0]/../dir_time_directory[n]/filename
    :type dir_time_descriptor: list, optional
    :return: path string
    :rtype: str
    """
    filename = f'{product_class.file_prefix}{time_step.strftime(product_class.file_mid)}_all{product_class.file_suffix}.nc'
    if dir_time_descriptor is not None:
        dir_time = ''
        for dir_act in dir_time_descriptor:
            dir_time = os.path.join(dir_time, time_step.strftime(dir_act))
    else:
        dir_time = ''
    return os.path.join(directory, dir_time, filename)
