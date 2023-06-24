#   @author: Michael Wagner
#   @organization: TU Dresden
#   @contact: michael.wagner@tu-dresden.de

#   @author: Michael Wagner
#   @organization: TU Dresden
#   @contact: michael.wagner@tu-dresden.de

import datetime as dt


class Product:
    """
    Product class contains all relevant technical data for download of RadolanRW, RadvorRQ, RadolanRV, IconD2, IconD2EPS
    files (from the German Weather Service - DWD).
    """
    def __init__(self,
                 update_time: dt.timedelta = None,
                 base_time: dt.timedelta = dt.timedelta(minutes=0),
                 url_folder: str = None,
                 url_time_subfolder: str = None,
                 url_subfolder: str = None,
                 file_prefix: str = None,
                 file_mid: str = None,
                 file_mid_aux: list = None,
                 file_suffix: str = None,
                 file_extension: str = None):
        """
        Initialize Product class.

        :param update_time: update time at DWD
        :type update_time: datetime.timedelta, optional
        :param base_time: the product's base time, e.g. 50 min for RadolanRW or 0 min for IconD2
        :type base_time: datetime.timedelta, optional
        :param url_folder: url folder before the first time dependent part of the file url, e.g.
            https://opendata.dwd.de/weather/nwp/icon-d2/grib/ in case of IconD2
        :type url_folder: str, optional
        :param url_time_subfolder: datetime.strftime descriptor for the time dependent url folder part, e.g. '%H' for
            IconD2
        :type url_time_subfolder: str, optional
        :param url_subfolder: url folder after a potential time dependent url folder part (requires url_time_subfolder),
            e.g. 'tot_prec' for IconD2
        :type url_subfolder: str, optional
        :param file_prefix: file prefix, e.g. 'raa01-rw_10000-' for RadolanRW data
        :type file_prefix: str, optional
        :param file_mid: datetime.strftime descriptor for the time dependent file name part, e.g. '%y%m%d%H%M' for
            RadolanRW
        :type file_mid: str, optional
        :param file_mid_aux: auxiliary part after time dependent file name part as a list of strings, e.g.
            ['_000', '_060', '_120'] for RadvorRQ
        :type file_mid_aux: list, optional
        :param file_suffix: file suffix after middle part, e.g. '-dwd---bin' for recent part of RadolanRW
        :type file_suffix: str, optional
        :param file_extension: extension of file, e.g. `.gz` for RadolanRW or `.grib2.bz2` for IconD2
        :type file_extension: str, optional
        """
        self.update_time = update_time
        self.base_time = base_time
        self.url_folder = url_folder
        self.url_time_subfolder = url_time_subfolder
        self.url_subfolder = url_subfolder
        self.file_prefix = file_prefix
        self.file_mid = file_mid
        self.file_mid_aux = file_mid_aux
        self.file_suffix = file_suffix
        self.file_extension = file_extension
