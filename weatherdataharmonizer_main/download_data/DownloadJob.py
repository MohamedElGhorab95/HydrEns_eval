#   @author: Michael Wagner
#   @organization: TU Dresden
#   @contact: michael.wagner@tu-dresden.de

import requests
import os
import datetime as dt

from aux_tools.products import *


class DownloadJob:
    """
    DownloadJob class contains all data and methods to download RadolanRW, RadvorRQ, RadolanRV, IconD2, IconD2EPS,
    IconEU, IconEUEPS data from https://opendata.dwd.de/.
    """
    def __init__(self,
                 product: str,
                 directory: str,
                 date_start: dt.datetime = None,
                 date_end: dt.datetime = dt.datetime.now(tz=dt.timezone.utc)):
        """
        Initialize DownloadJob class.

        :param product: supported products are RadolanRW, RadvorRQ, RadolanRV, IconD2, IconD2EPS, IconEU, IconEUEPS
        :type product: str
        :param directory: target directory where the downloaded files are saved to
        :type directory: str
        :param date_start: first date for download, default is 3 product update time cycles before date_end
        :type date_start: datetime.datetime, optional
        :param date_end: last date for download, default is now
        :type date_end: datetime.datetime, optional
        """
        self.product = product
        self.date_start = date_start
        self.date_end = dt.datetime(date_end.year, date_end.month, date_end.day, date_end.hour, date_end.minute,
                                    tzinfo=dt.timezone.utc)
        self.directory = directory

        if self.product.lower() == 'radolanrw':
            self.product_data = RadolanRWProp()
        elif self.product.lower() == 'radvorrq':
            self.product_data = RadvorRQProp()
        elif self.product.lower() == 'radolanrv':
            self.product_data = RadolanRVProp()
        elif self.product.lower() == 'icond2':
            self.product_data = IconD2Prop()
        elif self.product.lower() == 'icond2eps':
            self.product_data = IconD2EPSProp()
        elif self.product.lower() == 'iconeu':
            self.product_data = IconEUProp()
        elif self.product.lower() == 'iconeueps':
            self.product_data = IconEUEPSProp()
        else:
            raise Exception(f'product {self.product} not supported')

        self.files_stored = None

        if self.date_start is None:
            # if no specific start is given, a start of self.date_start - 2 * update_time is defined
            self.date_start = self.date_end - 2 * self.product_data.update_time

    def download_files(self):
        """
        Downloads the files. If the start or the end date do not match dates of delivery at DWD the dates are floored to
        the next provision date.
        """
        # identify most recent starts/ends at or before self.date_start/end
        possible_times = [self.product_data.base_time + x * self.product_data.update_time for x in
                          range(int(dt.timedelta(hours=24) / self.product_data.update_time))]
        while dt.timedelta(hours=self.date_start.hour, minutes=self.date_start.minute) not in possible_times:
            self.date_start = self.date_start - dt.timedelta(minutes=1)
        while dt.timedelta(hours=self.date_end.hour, minutes=self.date_end.minute) not in possible_times:
            self.date_end = self.date_end - dt.timedelta(minutes=1)

        # potential times for file download
        potential_times = [self.date_start]
        while potential_times[-1] < self.date_end:
            potential_times.append(potential_times[-1] + self.product_data.update_time)

        # build filenames
        time_subfolder_names = None
        if self.product_data.file_mid_aux is not None:
            # auxiliary part in the middle of the file name
            file_names = [self.product_data.file_prefix + t.strftime(self.product_data.file_mid) + mid_aux +
                          self.product_data.file_suffix + self.product_data.file_extension for t in potential_times
                          for mid_aux in self.product_data.file_mid_aux]
            if self.product_data.url_time_subfolder is not None:
                # subfolder naming using time
                time_subfolder_names = [t.strftime(self.product_data.url_time_subfolder) for t in potential_times
                                   for x in self.product_data.file_mid_aux]

        else:
            file_names = [self.product_data.file_prefix + t.strftime(self.product_data.file_mid) +
                          self.product_data.file_suffix + self.product_data.file_extension for t in potential_times]
            if self.product_data.url_time_subfolder is not None:
                # subfolder naming using time
                time_subfolder_names = [t.strftime(self.product_data.url_time_subfolder) for t in potential_times]

        # build urls
        if time_subfolder_names is not None:
            if self.product_data.url_subfolder is not None:
                # constant subfolder after time dependent subfolder
                subfolder_names = [os.path.join(x, self.product_data.url_subfolder) for x in time_subfolder_names]
            else:
                subfolder_names = time_subfolder_names
        else:
            # get list
            subfolder_names = ['' for x in file_names]

        file_urls = [os.path.join(self.product_data.url_folder, subfolder_names[x], file_names[x]) for x in
                     range(len(file_names))]

        # build output filenames
        self.files_stored = [os.path.join(self.directory, file_name) for file_name in file_names]

        # download files
        for file_ct in range(len(file_names)):
            if os.path.exists(self.files_stored[file_ct]):
                print(f'file {self.files_stored[file_ct]} already exists; skip download')
                continue

            print(f'download {file_urls[file_ct]}')
            r = requests.get(file_urls[file_ct], stream=True)
            if r.status_code == requests.codes.ok:
                with open(self.files_stored[file_ct], 'wb') as file:
                    for chunk in r.iter_content(chunk_size=1024):
                        # writing one chunk at a time to target file to avoid large ram consumption
                        if chunk:
                            file.write(chunk)
            else:
                print(f'\twarning: url not available')

    def delete_files(self):
        """
        Delete all downloaded files from the current DownloadJob instance.
        """
        for file in self.files_stored:
            if os.path.exists(file):
                os.remove(file)
