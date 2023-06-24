#   @author: Michael Wagner
#   @organization: TU Dresden
#   @contact: michael.wagner@tu-dresden.de
import datetime as dt

import numpy as np
import matplotlib.pyplot as plt
import netCDF4 as nc

import read_radolan.readRadolan as readRadolan
import read_icon.readIcon as readIcon
import read_cosmo.readCosmo as readCosmo
from download_data.DownloadJob import DownloadJob
from met_entities.CosmoD2 import CosmoD2
from met_entities.CosmoD2EPS import CosmoD2EPS
from met_entities.RadolanRV import RadolanRV
from met_entities.RadvorRQ import RadvorRQ
from met_entities.RadolanRW import RadolanRW
from met_entities.IconD2 import IconD2
from met_entities.IconD2EPS import IconD2EPS
from met_entities.IconEU import IconEU
from met_entities.IconEUEPS import IconEUEPS
from met_entities.VariableDescription import CropDescription, RegridDescription
from met_entities.WeatherData import WeatherData
from met_entities import LonLatTime as LonLatTime
from met_entities.read_netcdf import read_netcdf
import time

lon = np.arange(4, 15, .3)
lat = np.arange(54, 48, -.3)
lon_target, lat_target = np.meshgrid(lon, lat)
lon_target = LonLatTime.LonLatTime(data=lon_target)
lat_target = LonLatTime.LonLatTime(data=lat_target)

def main_radolanrw():
    # RadolanRW example
    start_datetime = dt.datetime(2022, 12, 14, 6, 50)
    radrw = RadolanRW()
    # radrw.read_file(filename='data/radolanrw/raa01-rw_10000-2212140650-dwd---bin.bz2', scale_factor=0.1, fill_value=-999)
    radrw.read_file(start_datetime=start_datetime, directory='data/radolanrw', scale_factor=1, fill_value=-1)
    radrw.regrid(lon_target=lon_target, lat_target=lat_target, file_nearest='data/radrw_regridding.npz')
    # radrw.crop(lon_west=11.7, lon_east=15.2, lat_south=50.1, lat_north=51.8)
    radrw.export_netcdf('data/radrw_example.nc', data_format='i2', scale_factor_nc=0.1)
    # rd = readRadolan.read_radolan_binary('data/radolanrw/raa01-rw_10000-2210250650-dwd---bin')

    # continuous RadolanRW example
    start_datetime = dt.datetime(2022, 12, 14, 5, 50)
    end_datetime = dt.datetime(2022, 12, 14, 13, 50)
    period_step = dt.timedelta(hours=1)
    periods = [start_datetime]
    ct = 0
    while periods[-1] < end_datetime:
        ct = ct + 1
        periods.append(start_datetime + ct * period_step)
    filenames = []
    mode = 'create'
    for period in periods:
        print(period)
        filenames.append(f"data/radolanrw/raa01-rw_10000-{period.strftime('%y%m%d%H%M')}-dwd---bin.bz2")
        radrw = RadolanRW()
        radrw.read_file(start_datetime=period, directory='data/radolanrw', scale_factor=1, fill_value=-1)
        radrw.crop(lon_west=11.7, lon_east=15.2, lat_south=50.1, lat_north=51.8)
        if mode == 'create':
            radrw.export_netcdf('data/radrw_cont_example.nc', data_format='i2', scale_factor_nc=0.1)
        else:
            radrw.export_netcdf_append('data/radrw_cont_example.nc')
        mode = 'append'


def main_radvorrq():
    # RadvorRQ example
    start_datetime = dt.datetime(2022, 12, 14, 12, 30)
    radrq = RadvorRQ()
    radrq.read_file(start_datetime=start_datetime, directory='data/radvorrq', scale_factor=1, fill_value=-1)
    radrq.regrid(lon_target=lon_target, lat_target=lat_target, file_nearest='data/radrq_regridding.npz')
    # radrq.crop(lon_west=11.7, lon_east=15.2, lat_south=50.1, lat_north=51.8)
    # radrq.crop(idx_west=635, idx_east=899, idx_south=556, idx_north=324)
    radrq.export_netcdf("data/radrq_example.nc", data_format='i2', scale_factor_nc=0.1)

    # continuous RadvorRQ example
    start_datetime = dt.datetime(2022, 12, 14, 12, 30)
    end_datetime = dt.datetime(2022, 12, 14, 14, 0)
    period_step = dt.timedelta(minutes=15)
    periods = [start_datetime]
    ct = 0
    while periods[-1] < end_datetime:
        ct = ct + 1
        periods.append(start_datetime + ct * period_step)
    mode = 'create'
    for period in periods:
        print(period)
        radrq = RadvorRQ()
        radrq.read_file(start_datetime=period, directory='data/radvorrq', scale_factor=1, fill_value=-1)
        radrq.crop(lon_west=11.7, lon_east=15.2, lat_south=50.1, lat_north=51.8)
        if mode == 'create':
            radrq.export_netcdf('data/radrq_cont_example.nc', data_format='i2', scale_factor_nc=0.1)
        else:
            radrq.export_netcdf_append('data/radrq_cont_example.nc')
        mode = 'append'

def main_radolanrv():
    # RadolanRV example
    start_datetime = dt.datetime(2022, 12, 14, 12, 30)
    radrv = RadolanRV()
    radrv.read_file(start_datetime=start_datetime, directory='data/radolanrv', scale_factor=1, fill_value=-1)
    # radrv.read_file("data/radolanrv/DE1200_RV2212141230.tar.bz2", scale_factor=0.1, fill_value=-999)
    radrv.regrid(lon_target=lon_target, lat_target=lat_target, file_nearest='data/radrv_regridding.npz')
    # radrv.crop(lon_west=11.7, lon_east=15.2, lat_south=50.1, lat_north=51.8)
    radrv.export_netcdf("data/radrv_example.nc", data_format='i2', scale_factor_nc=0.1)

    # continuous RadolanRV example
    start_datetime = dt.datetime(2022, 12, 14, 12, 30, 0)
    end_datetime = dt.datetime(2022, 12, 14, 13, 10, 0)
    period_step = dt.timedelta(minutes=5)
    periods = [start_datetime]
    ct = 0
    while periods[-1] < end_datetime:
        ct = ct + 1
        periods.append(start_datetime + ct * period_step)
    filenames = []
    mode = 'create'
    for period in periods:
        print(period)
        filenames.append(f"data/radolanrv/DE1200_RV{period.strftime('%y%m%d%H%M')}.tar.bz2")
        radrv = RadolanRV()
        radrv.read_file(filenames[-1], scale_factor=1, fill_value=-1)
        radrv.crop(lon_west=11.7, lon_east=15.2, lat_south=50.1, lat_north=51.8)
        if mode == 'create':
            radrv.export_netcdf('data/radrv_cont_example.nc', data_format='i2', scale_factor_nc=0.1)
        else:
            radrv.export_netcdf_append('data/radrv_cont_example.nc')
        mode = 'append'

def main_icond2():
    # IconD2 example
    start_datetime = dt.datetime(2022, 12, 14, 9)
    icond2 = IconD2()
    icond2.read_file(start_datetime=start_datetime, directory='data/icon_d2', forecast_hours=3, scale_factor=1,
                     fill_value=-1)
    icond2.regrid(lon_target=lon_target, lat_target=lat_target, file_nearest='data/icond2_regridding.npz')
    # icond2.crop(lon_west=11.7, lon_east=15.2, lat_south=50.1, lat_north=51.8)
    icond2.export_netcdf('data/icond2_example.nc', data_format='i2', scale_factor_nc=0.1)

    # continuous IconD2 example
    start_datetime = dt.datetime(2022, 12, 14, 9)
    end_datetime = dt.datetime(2022, 12, 14, 15)
    period_step = dt.timedelta(hours=3)
    periods = [start_datetime]
    ct = 0
    while periods[-1] < end_datetime:
        ct = ct + 1
        periods.append(start_datetime + ct * period_step)
    mode = 'create'
    for period in periods:
        print(period)
        icond2 = IconD2()
        icond2.read_file(period, directory='data/icon_d2', forecast_hours=2, scale_factor=1, fill_value=-1)
        icond2.regrid(lon_target=lon_target, lat_target=lat_target, file_nearest='data/icond2_regridding.npz')
        if mode == 'create':
            icond2.export_netcdf('data/icond2_cont_example.nc', data_format='i2', scale_factor_nc=0.1)
        else:
            icond2.export_netcdf_append('data/icond2_cont_example.nc')
        mode = 'append'

def main_icond2eps():
    # IconD2EPS example
    start_datetime = dt.datetime(2022, 12, 14, 9)
    icond2eps = IconD2EPS()
    icond2eps.read_file(start_datetime=start_datetime, directory='data/icon_d2_eps', forecast_hours=1,
                        eps_member=[0, 1, 5], scale_factor=1, fill_value=-1)
    icond2eps.regrid(lon_target=lon_target, lat_target=lat_target, file_nearest='data/icond2_regridding.npz')
    # icond2eps.crop(lon_west=11.7, lon_east=15.2, lat_south=50.1, lat_north=51.8)
    icond2eps.export_netcdf('data/icond2eps_example.nc', data_format='i2', version='separated', scale_factor_nc=0.1)

    # continuous IconD2EPS example
    start_datetime = dt.datetime(2022, 12, 14, 9, 0, 0)
    end_datetime = dt.datetime(2022, 12, 14, 15, 0, 0)
    period_step = dt.timedelta(hours=3)
    periods = [start_datetime]
    ct = 0
    while periods[-1] < end_datetime:
        ct = ct + 1
        periods.append(start_datetime + ct * period_step)
    mode = 'create'
    for period in periods:
        print(period)
        icond2eps = IconD2EPS()
        icond2eps.read_file(period, directory='data/icon_d2_eps', forecast_hours=0, scale_factor=1, fill_value=-1)
        icond2eps.regrid(lon_target=lon_target, lat_target=lat_target, file_nearest='data/icond2eps_regridding.npz')
        if mode == 'create':
            icond2eps.export_netcdf('data/icond2eps_cont_example.nc', data_format='i2', version='separated',
                                    scale_factor_nc=0.1)
        else:
            icond2eps.export_netcdf_append('data/icond2eps_cont_example.nc')
        mode = 'append'

def main_icond2eps_second():
    # second IconD2EPS example with more arguments and a larger grid
    lon_target_act, lat_target_act = readRadolan.get_lonlat(4, 'radolanrx')
    start_datetime = dt.datetime(2022, 12, 14, 9)
    icond2eps = IconD2EPS()
    tic = time.time()
    icond2eps.read_file(start_datetime=start_datetime, directory='data/icon_d2_eps', forecast_hours=0,
                        scale_factor=0.1, fill_value=-1, short='int16', eps_member=[0])
    # icond2eps.regrid(lon_target=lon_target_act, lat_target=lat_target_act,
    #                  file_nearest='data/idw_rules/rad4_regridding_icond2.npz')
    print(f'{time.time() - tic} s')
    tic = time.time()
    data_kwargs = {'compression': 'zlib', 'complevel': 1}
    icond2eps.export_netcdf('data/icond2eps_example_i2.nc', data_format='i2', version='separated',
                            data_kwargs=data_kwargs, scale_factor_nc=0.1, scale_undo=True)
    print(f'{time.time() - tic} s')

def main_icond2eps_third():
    # third IconD2EPS example with reading in netcdf file
    tic = time.time()
    icond2eps = read_netcdf('data/icond2eps_example.nc')
    print(f'{time.time() - tic} s')
    icond2eps.export_netcdf('data/icond2eps_example2.nc', data_format='i4', version='separated', scale_factor_nc=0.01,
                            scale_undo=True, data_kwargs={'compression':'zlib', 'complevel':1})
    print(f'{time.time() - tic} s')

def main_iconeu():
    # IconEU example
    start_datetime = dt.datetime(2023, 3, 23, 0)
    iconeu = IconEU()
    iconeu.read_file(start_datetime=start_datetime, directory='data/icon_eu', forecast_hours=10, scale_factor=1,
                     fill_value=-1)
    # iconeu.regrid(lon_target=lon_target, lat_target=lat_target, file_nearest='data/iconeu_regridding.npz')
    iconeu.crop(lon_west=11.7, lon_east=15.2, lat_south=50.1, lat_north=51.8)
    iconeu.export_netcdf('data/iconeu_example.nc', data_format='i2', scale_factor_nc=0.1)

    # continuous IconEU example
    start_datetime = dt.datetime(2023, 3, 23, 0)
    end_datetime = dt.datetime(2023, 3, 23, 6)
    period_step = dt.timedelta(hours=3)
    periods = [start_datetime]
    ct = 0
    while periods[-1] < end_datetime:
        ct = ct + 1
        periods.append(start_datetime + ct * period_step)
    mode = 'create'
    for period in periods:
        print(period)
        iconeu = IconEU()
        iconeu.read_file(period, directory='data/icon_eu', forecast_hours=3, scale_factor=1, fill_value=-1)
        iconeu.regrid(lon_target=lon_target, lat_target=lat_target, file_nearest='data/iconeu_regridding.npz')
        if mode == 'create':
            iconeu.export_netcdf('data/iconeu_cont_example.nc', data_format='i2', scale_factor_nc=0.1)
        else:
            iconeu.export_netcdf_append('data/iconeu_cont_example.nc')
        mode = 'append'

def main_iconeueps():
    # IconEUEPS example
    start_datetime = dt.datetime(2023, 3, 23, 0)
    iconeueps = IconEUEPS()
    iconeueps.read_file(start_datetime=start_datetime, directory='data/icon_eu_eps', forecast_hours=5,
                        eps_member=[0, 1, 5], scale_factor=1, fill_value=-1)
    iconeueps.regrid(lon_target=lon_target, lat_target=lat_target, file_nearest='data/iconeueps_regridding.npz')
    # iconeueps.crop(lon_west=11.7, lon_east=15.2, lat_south=50.1, lat_north=51.8)
    iconeueps.export_netcdf('data/iconeueps_example.nc', data_format='i2', version='separated', scale_factor_nc=0.1)

    # continuous IconEUEPS example
    start_datetime = dt.datetime(2023, 3, 23, 0)
    end_datetime = dt.datetime(2023, 3, 23, 6)
    period_step = dt.timedelta(hours=6)
    periods = [start_datetime]
    ct = 0
    while periods[-1] < end_datetime:
        ct = ct + 1
        periods.append(start_datetime + ct * period_step)
    mode = 'create'
    for period in periods:
        print(period)
        iconeueps = IconEUEPS()
        iconeueps.read_file(period, directory='data/icon_eu_eps', forecast_hours=3, scale_factor=1, fill_value=-1)
        iconeueps.regrid(lon_target=lon_target, lat_target=lat_target, file_nearest='data/iconeueps_regridding.npz')
        if mode == 'create':
            iconeueps.export_netcdf('data/iconeueps_cont_example.nc', data_format='i2', version='separated',
                                    scale_factor_nc=0.1)
        else:
            iconeueps.export_netcdf_append('data/iconeueps_cont_example.nc')
        mode = 'append'

def main_cosmod2():
    # CosmoD2 example
    start_datetime = dt.datetime(2021, 1, 17, 0)
    cosmod2 = CosmoD2()
    cosmod2.read_file(start_datetime=start_datetime, directory='data/cosmo_d2', forecast_hours=3, scale_factor=1,
                      fill_value=-1)
    cosmod2.regrid(lon_target=lon_target, lat_target=lat_target, file_nearest='data/cosmod2_regridding.npz')
    # cosmod2.crop(lon_west=11.7, lon_east=15.2, lat_south=50.1, lat_north=51.8)
    cosmod2.export_netcdf('data/cosmod2_example.nc', data_format='i2', scale_factor_nc=0.1)

    # continuous CosmoD2 example
    start_datetime = dt.datetime(2021, 1, 17, 0)
    end_datetime = dt.datetime(2021, 1, 17, 6)
    period_step = dt.timedelta(hours=3)
    periods = [start_datetime]
    ct = 0
    while periods[-1] < end_datetime:
        ct = ct + 1
        periods.append(start_datetime + ct * period_step)
    mode = 'create'
    for period in periods:
        print(period)
        cosmod2 = CosmoD2()
        cosmod2.read_file(period, directory='data/cosmo_d2', forecast_hours=3, scale_factor=1, fill_value=-1)
        cosmod2.regrid(lon_target=lon_target, lat_target=lat_target, file_nearest='data/cosmod2_regridding.npz')
        if mode == 'create':
            cosmod2.export_netcdf('data/cosmod2_cont_example.nc', data_format='i2', scale_factor_nc=0.1)
        else:
            cosmod2.export_netcdf_append('data/cosmod2_cont_example.nc')
        mode = 'append'

def main_cosmod2eps():
    # CosmoD2EPS example
    start_datetime = dt.datetime(2021, 1, 17, 0)
    cosmod2eps = CosmoD2EPS()
    cosmod2eps.read_file(start_datetime=start_datetime, directory='data/cosmo_d2_eps', forecast_hours=1,
                         eps_member=[0, 1, 5], scale_factor=1, fill_value=-1)
    cosmod2eps.regrid(lon_target=lon_target, lat_target=lat_target, file_nearest='data/cosmod2_regridding.npz')
    # icond2eps.crop(lon_west=11.7, lon_east=15.2, lat_south=50.1, lat_north=51.8)
    cosmod2eps.export_netcdf('data/cosmod2eps_example.nc', data_format='i2', version='separated', scale_factor_nc=0.1)

    # continuous CosmoD2EPS example
    start_datetime = dt.datetime(2021, 1, 17, 0)
    end_datetime = dt.datetime(2021, 1, 17, 6)
    period_step = dt.timedelta(hours=3)
    periods = [start_datetime]
    ct = 0
    while periods[-1] < end_datetime:
        ct = ct + 1
        periods.append(start_datetime + ct * period_step)
    mode = 'create'
    for period in periods:
        print(period)
        cosmod2eps = CosmoD2EPS()
        cosmod2eps.read_file(period, directory='data/cosmo_d2_eps', forecast_hours=0, scale_factor=1, fill_value=-1)
        cosmod2eps.regrid(lon_target=lon_target, lat_target=lat_target, file_nearest='data/cosmo_regridding.npz')
        if mode == 'create':
            cosmod2eps.export_netcdf('data/cosmod2eps_cont_example.nc', data_format='i2', version='separated',
                                     scale_factor_nc=0.1)
        else:
            cosmod2eps.export_netcdf_append('data/cosmod2eps_cont_example.nc')
        mode = 'append'

def main_combined_data():
    # example with combined data (RadRW + RadRQ + RadRV + IconD2/EPS)
    time_now = dt.datetime(2022, 12, 14, 13, tzinfo=dt.timezone.utc)
    # crop_description = CropDescription(lon_west=11.7, lon_east=15.2, lat_south=50.1, lat_north=51.8)
    regrid_description = {'radolanrw': RegridDescription(lon_target=lon_target, lat_target=lat_target,
                                                         file_nearest='data/radrw_regridding.npz'),
                          'radvorrq': RegridDescription(lon_target=lon_target, lat_target=lat_target,
                                                        file_nearest='data/radrq_regridding.npz'),
                          'radolanrv': RegridDescription(lon_target=lon_target, lat_target=lat_target,
                                                        file_nearest='data/radrv_regridding.npz'),
                          'icond2': RegridDescription(lon_target=lon_target, lat_target=lat_target,
                                                      file_nearest='data/icond2_regridding.npz'),
                          'icond2eps': RegridDescription(lon_target=lon_target, lat_target=lat_target,
                                                      file_nearest='data/icond2eps_regridding.npz')}
    wd = WeatherData(time_now=time_now, delta_t=dt.timedelta(minutes=15), fill_value=-1, scale_factor=0.1,
                     regrid_description=regrid_description, short='int16')

    wd.collect_radolanrw(time_start=time_now - dt.timedelta(hours=2), directory='data/radolanrw')
    # wd.collect_icond2eps(latest_event=time_now + dt.timedelta(hours=1), directory='data/icon_d2_eps', forecast_hours=1,
    #                      eps_member=range(2))
    wd.collect_icond2(latest_event=time_now + dt.timedelta(hours=1), directory='data/icon_d2', forecast_hours=2)
    wd.collect_radvorrq(latest_event=time_now + dt.timedelta(hours=0.25), directory='data/radvorrq')
    # wd.collect_radolanrv(latest_event=time_now + dt.timedelta(hours=0.25), directory='data/radolanrv')

    filename_nc = f'data/precipitation_data_example_{time_now.strftime("%Y%m%d%H%M")}.nc'
    wd.export_netcdf(filename=filename_nc, institution='TU Dresden, Institute of Hydrology and Meteorology',
                     data_format='i2', scale_factor_nc=0.1, scale_undo=True)

def main_append_combined_data():
    time_now = dt.datetime(2022, 12, 14, 14, tzinfo=dt.timezone.utc)
    regrid_description = {'radolanrw': RegridDescription(lon_target=lon_target, lat_target=lat_target,
                                                         file_nearest='data/radrw_regridding.npz'),
                          'radvorrq': RegridDescription(lon_target=lon_target, lat_target=lat_target,
                                                        file_nearest='data/radrq_regridding.npz'),
                          'radolanrv': RegridDescription(lon_target=lon_target, lat_target=lat_target,
                                                         file_nearest='data/radrv_regridding.npz'),
                          'icond2': RegridDescription(lon_target=lon_target, lat_target=lat_target,
                                                      file_nearest='data/icond2_regridding.npz'),
                          'icond2eps': RegridDescription(lon_target=lon_target, lat_target=lat_target,
                                                         file_nearest='data/icond2eps_regridding.npz')}
    wd = WeatherData(time_now=time_now, delta_t=dt.timedelta(minutes=15), fill_value=-1, scale_factor=0.1,
                     regrid_description=regrid_description, short='int16')
    wd.collect_radolanrw(time_start=time_now - dt.timedelta(hours=1), directory='data/radolanrw')
    wd.collect_icond2(latest_event=time_now + dt.timedelta(hours=1), directory='data/icon_d2', forecast_hours=4)

    filename_nc = 'data/precipitation_data_example_202212141300.nc'
    wd.export_netcdf_append(filename=filename_nc)

def main_combined_data2():
    # example with combined data (RadRW + IconEU)
    time_now = dt.datetime(2023, 3, 23, 0, tzinfo=dt.timezone.utc)
    regrid_description = {'radolanrw': RegridDescription(lon_target=lon_target, lat_target=lat_target,
                                                         file_nearest='data/radrw_regridding.npz'),
                          'radvorrq': RegridDescription(lon_target=lon_target, lat_target=lat_target,
                                                        file_nearest='data/radrq_regridding.npz'),
                          'radolanrv': RegridDescription(lon_target=lon_target, lat_target=lat_target,
                                                        file_nearest='data/radrv_regridding.npz'),
                          'icond2': RegridDescription(lon_target=lon_target, lat_target=lat_target,
                                                      file_nearest='data/icond2_regridding.npz'),
                          'icond2eps': RegridDescription(lon_target=lon_target, lat_target=lat_target,
                                                      file_nearest='data/icond2eps_regridding.npz'),
                          'iconeu': RegridDescription(lon_target=lon_target, lat_target=lat_target,
                                                      file_nearest='data/iconeu_regridding.npz'),
                          'iconeueps': RegridDescription(lon_target=lon_target, lat_target=lat_target,
                                                         file_nearest='data/iconeueps_regridding.npz')}
    wd = WeatherData(time_now=time_now, delta_t=dt.timedelta(minutes=15), fill_value=-1, scale_factor=1,
                     regrid_description=regrid_description)

    wd.collect_radolanrw(time_start=time_now - dt.timedelta(hours=2), directory='data/radolanrw')
    wd.collect_iconeu(latest_event=time_now, directory='data/icon_eu')
    # wd.collect_iconeueps(latest_event=time_now, directory='data/icon_eu_eps', eps_member=range(10))

    filename_nc = f'data/precipitation_data_example_{time_now.strftime("%Y%m%d%H%M")}.nc'
    wd.export_netcdf(filename=filename_nc, institution='TU Dresden, Institute of Hydrology and Meteorology',
                     data_format='i2', scale_factor_nc=0.1, scale_undo=True)

def main_download_data():
    # examplary usage of download from opendata
    dj = DownloadJob(product='RadolanRW',
                     directory='data/radolanrw',
                     date_start=dt.datetime(2023, 3, 22, 20, 50, tzinfo=dt.timezone.utc))
    dj.download_files()
    # dj.delete_files()

def main_get_coordinates():
    coord_radrw = readRadolan.get_lonlat(4, 'radolanrx')
    coord_radrq = readRadolan.get_lonlat(5, 'radolanrx')
    coord_radrv = readRadolan.get_lonlat(5, 'de1200')
    coord_icond2 = readIcon.get_lonlat('d2')
    coord_iconeu = readIcon.get_lonlat('eu')
    coord_iconeueps = readIcon.get_lonlat('eueps')
    coord_cosmo = readCosmo.get_lonlat()

    with open('coord_radrw.txt', 'w') as fid:
        fid.write('lon\tlat\n')
        lon_rw = np.ravel(coord_radrw[0].data)
        lat_rw = np.ravel(coord_radrw[1].data)
        for ct in range(len(lon_rw)):
            fid.write(f'{lon_rw[ct]:.4f}\t{lat_rw[ct]:.4f}\n')

    with open('coord_radrq.txt', 'w') as fid:
        fid.write('lon\tlat\n')
        lon_rq = np.ravel(coord_radrq[0].data)
        lat_rq = np.ravel(coord_radrq[1].data)
        for ct in range(len(lon_rq)):
            fid.write(f'{lon_rq[ct]:.4f}\t{lat_rq[ct]:.4f}\n')

    with open('coord_radrv.txt', 'w') as fid:
        fid.write('lon\tlat\n')
        lon_rv = np.ravel(coord_radrv[0].data)
        lat_rv = np.ravel(coord_radrv[1].data)
        for ct in range(len(lon_rv)):
            fid.write(f'{lon_rv[ct]:.4f}\t{lat_rv[ct]:.4f}\n')

    with open('coord_icond2.txt', 'w') as fid:
        fid.write('lon\tlat\n')
        lon_icon = np.ravel(coord_icond2[0].data)
        lat_icon = np.ravel(coord_icond2[1].data)
        for ct in range(len(lon_icon)):
            fid.write(f'{lon_icon[ct]:.4f}\t{lat_icon[ct]:.4f}\n')

    with open('coord_iconeu.txt', 'w') as fid:
        fid.write('lon\tlat\n')
        lon_iconeu = np.ravel(coord_iconeu[0].data)
        lat_iconeu = np.ravel(coord_iconeu[1].data)
        for ct in range(len(lon_iconeu)):
            fid.write(f'{lon_iconeu[ct]:.4f}\t{lat_iconeu[ct]:.4f}\n')

    with open('coord_iconeueps.txt', 'w') as fid:
        fid.write('lon\tlat\n')
        lon_iconeueps = np.ravel(coord_iconeueps[0].data)
        lat_iconeueps = np.ravel(coord_iconeueps[1].data)
        for ct in range(len(lon_iconeueps)):
            fid.write(f'{lon_iconeueps[ct]:.4f}\t{lat_iconeueps[ct]:.4f}\n')

    with open('coord_cosmo.txt', 'w') as fid:
        fid.write('lon\tlat\n')
        lon_cosmo = np.ravel(coord_cosmo[0].data)
        lat_cosmo = np.ravel(coord_cosmo[1].data)
        for ct in range(len(lon_cosmo)):
            fid.write(f'{lon_cosmo[ct]:.4f}\t{lat_cosmo[ct]:.4f}\n')


if __name__ == '__main__':
    # main_radolanrw()
    # main_radvorrq()
    # main_radolanrv()
    # main_icond2()
    # main_icond2eps()
    # main_icond2eps_second()
    # main_icond2eps_third()
    # main_iconeu()
    # main_iconeueps()
    # main_cosmod2()
    # main_cosmod2eps()
    # main_combined_data()
    # main_append_combined_data()
    main_combined_data2()
    # main_download_data()
    # main_get_coordinates()
