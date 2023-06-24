#   @author: Michael Wagner
#   @organization: TU Dresden
#   @contact: michael.wagner@tu-dresden.de
import datetime as dt
import numpy as np

from aux_tools.Product import Product


class RadolanRWProp(Product):
    """
    Instantiate Product class with data for RadolanRW.
    """
    def __init__(self):
        super().__init__(update_time=dt.timedelta(hours=1),
                         base_time=dt.timedelta(minutes=50),
                         url_folder='https://opendata.dwd.de/climate_environment/CDC/grids_germany/hourly/radolan/recent/bin/',
                         file_prefix='raa01-rw_10000-',
                         file_mid='%y%m%d%H%M',
                         file_suffix='-dwd---bin',
                         file_extension='.gz')


class RadvorRQProp(Product):
    """
    Instantiate Product class with data for RadvorRQ.
    """
    def __init__(self):
        super().__init__(update_time=dt.timedelta(minutes=15),
                         url_folder='https://opendata.dwd.de/weather/radar/radvor/rq/',
                         file_prefix='RQ',
                         file_mid='%y%m%d%H%M',
                         file_mid_aux=['_000', '_060', '_120'],
                         file_suffix='',
                         file_extension='.gz')


class RadolanRVProp(Product):
    """
    Instantiate Product class with data for RadolanRV.
    """
    def __init__(self):
        super().__init__(update_time=dt.timedelta(minutes=5),
                         url_folder='https://opendata.dwd.de/weather/radar/composite/rv/',
                         file_prefix = 'DE1200_RV',
                         file_mid = '%y%m%d%H%M',
                         file_suffix = '',
                         file_extension = '.tar.bz2')


class IconD2Prop(Product):
    """
    Instantiate Product class with data for IconD2.
    """
    def __init__(self):
        super().__init__(update_time=dt.timedelta(hours=3),
                         url_folder='https://opendata.dwd.de/weather/nwp/icon-d2/grib/',
                         url_time_subfolder='%H',
                         url_subfolder='tot_prec',
                         file_prefix='icon-d2_germany_icosahedral_single-level_',
                         file_mid = '%Y%m%d%H',
                         file_mid_aux = [f'_{i:03d}' for i in range(49)],
                         file_suffix = '_2d_tot_prec',
                         file_extension = '.grib2.bz2')


class IconD2EPSProp(Product):
    """
    Instantiate Product class with data for IconD2EPS.
    """
    def __init__(self):
        super().__init__(update_time=dt.timedelta(hours=3),
                         url_folder = 'https://opendata.dwd.de/weather/nwp/icon-d2-eps/grib/',
                         url_time_subfolder='%H',
                         url_subfolder='tot_prec',
                         file_prefix = 'icon-d2-eps_germany_icosahedral_single-level_',
                         file_mid = '%Y%m%d%H',
                         file_mid_aux=[f'_{i:03d}' for i in range(49)],
                         file_suffix = '_2d_tot_prec',
                         file_extension = '.grib2.bz2')


class IconEUProp(Product):
    """
    Instantiate Product class with data for IconEU.
    """
    def __init__(self):
        super().__init__(update_time=dt.timedelta(hours=3),
                         url_folder='https://opendata.dwd.de/weather/nwp/icon-eu/grib/',
                         url_time_subfolder='%H',
                         url_subfolder='tot_prec',
                         file_prefix='icon-eu_europe_regular-lat-lon_single-level_',
                         file_mid = '%Y%m%d%H',
                         file_mid_aux = [f'_{i:03d}' for i in np.concatenate((np.arange(0, 79), np.arange(81, 121, 3)))],
                         file_suffix = '_TOT_PREC',
                         file_extension = '.grib2.bz2')


class IconEUEPSProp(Product):
    """
    Instantiate Product class with data for IconEUEPS.
    """
    def __init__(self):
        super().__init__(update_time=dt.timedelta(hours=6),
                         url_folder = 'https://opendata.dwd.de/weather/nwp/icon-eu-eps/grib/',
                         url_time_subfolder='%H',
                         url_subfolder='tot_prec',
                         file_prefix = 'icon-eu-eps_europe_icosahedral_single-level_',
                         file_mid = '%Y%m%d%H',
                         file_mid_aux=[f'_{i:03d}' for i in np.concatenate((np.arange(0, 49), np.arange(51, 73, 3), np.arange(78, 121, 6)))],
                         file_suffix = '_tot_prec',
                         file_extension = '.grib2.bz2')
