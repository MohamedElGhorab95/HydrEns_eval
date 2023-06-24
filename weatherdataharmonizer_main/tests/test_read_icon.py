#   @author: Michael Wagner
#   @organization: TU Dresden
#   @contact: michael.wagner@tu-dresden.de
import unittest

from read_icon import readIcon


class TestReadIcon(unittest.TestCase):
    def test_read_icon(self):
        """
        Test whether an Icon file is properly read in and provides the correct product name and the correct shape of
        data.
        """
        data = readIcon.read_icon_d2('data/icon-d2_germany_icosahedral_single-level_2022121412_000_2d_tot_prec.grib2.bz2')
        self.assertEqual(len(data), 4)
        self.assertEqual(data[0].metadata.name, 'Total Precipitation')
        self.assertEqual(data[0].data.shape, (542040,))

    def test_get_lonlat(self):
        """
        get_lonlat is tested for the right values and dimension.
        """
        lon, lat = readIcon.get_lonlat('d2')
        self.assertEqual(lon.data.shape, (542040,))
        self.assertEqual(lat.data.shape, (542040,))
        self.assertAlmostEqual(lon.data.min(), -4.16164, places=5)
        self.assertAlmostEqual(lon.data.max(), 20.54442, places=5)
        self.assertAlmostEqual(lat.data.min(), 43.04405, places=5)
        self.assertAlmostEqual(lat.data.max(), 58.16465, places=5)


