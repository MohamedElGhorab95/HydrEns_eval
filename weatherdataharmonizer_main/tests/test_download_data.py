#   @author: Michael Wagner
#   @organization: TU Dresden
#   @contact: michael.wagner@tu-dresden.de
import datetime as dt
import unittest
import os

import read_radolan
from download_data.DownloadJob import DownloadJob
from read_radolan import readRadolan as readRadolan


class TestDownloadJob(unittest.TestCase):
    def test_download_files(self):
        """
        Test of download functionality. Test for correct number of downloaded files after download and correct deleting.
        """
        date_end = dt.datetime.now(tz=dt.timezone.utc) - dt.timedelta(hours=3)
        date_start = date_end - dt.timedelta(hours=1)
        download_path = 'test_download'

        # test download of RadolanRW product
        dj = DownloadJob(product='RadolanRW', directory=download_path, date_start=date_start, date_end=date_end)
        dj.download_files()
        ct = 0
        for files in os.listdir(download_path):
            # check if current path is a file
            if os.path.isfile(os.path.join(download_path, files)):
                ct = ct + 1
        self.assertEqual(ct, 2)

        # exemplary test of content of a downloaded file
        files = os.listdir(download_path)
        data = readRadolan.read_radolan(filename=os.path.join(download_path, files[0]))
        self.assertEqual(data.data.shape, (900, 900))

        dj.delete_files()
        ct = 0
        for files in os.listdir(download_path):
            # check if current path is a file
            if os.path.isfile(os.path.join(download_path, files)):
                ct = ct + 1
        self.assertEqual(ct, 0)


        # test download of RadvorRQ product (additional middle part in the file naming convention)
        date_start = date_end
        dj = DownloadJob(product='RadvorRQ', directory=download_path, date_start=date_start, date_end=date_end)
        dj.download_files()
        ct = 0
        for files in os.listdir(download_path):
            # check if current path is a file
            if os.path.isfile(os.path.join(download_path, files)):
                ct = ct + 1
        self.assertEqual(ct, 3)

        dj.delete_files()
        ct = 0
        for files in os.listdir(download_path):
            # check if current path is a file
            if os.path.isfile(os.path.join(download_path, files)):
                ct = ct + 1
        self.assertEqual(ct, 0)
