#   @author: Michael Wagner
#   @organization: TU Dresden
#   @contact: michael.wagner@tu-dresden.de
import unittest

from read_cosmo import readCosmo


class TestReadCosmo(unittest.TestCase):
    def test_read_cosmo(self):
        """
        Test whether a Cosmo file is properly read in and provides the correct product name and the correct shape of
        data.
        """
        data = readCosmo.read_cosmo_d2('data/cosmo-d2_germany_rotated-lat-lon_single-level_2021011700_000_TOT_PREC.grib2.bz2')
        self.assertEqual(len(data), 4)
        self.assertEqual(data[0].metadata.name, 'Total Precipitation')
        self.assertEqual(data[0].data.shape, (716, 651))

    def test_get_lonlat(self):
        """
        get_lonlat is tested for the right values and dimension.
        """
        sw = [-0.249426, 43.185953]
        ne = [20.209766, 57.624364]
        lon, lat = readCosmo.get_lonlat()
        self.assertEqual(lon.data.shape, (716, 651))
        self.assertEqual(lat.data.shape, (716, 651))
        self.assertAlmostEqual(lon.data[-1, 0], sw[0], places=5)
        self.assertAlmostEqual(lat.data[-1, 0], sw[1], places=5)
        self.assertAlmostEqual(lon.data[0, -1], ne[0], places=5)
        self.assertAlmostEqual(lat.data[0, -1], ne[1], places=5)
