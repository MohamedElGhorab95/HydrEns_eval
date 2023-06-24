#   @author: Michael Wagner
#   @organization: TU Dresden
#   @contact: michael.wagner@tu-dresden.de
import unittest

from read_radolan import readRadolan


class TestReadRadolan(unittest.TestCase):
    def test_read_radolan(self):
        """
        Test whether a radolan file is properly read in and provides the correct product name and the correct shape of
        data.
        """
        data = readRadolan.read_radolan('data/raa01-rw_10000-2212141250-dwd---bin.gz')
        self.assertEqual(data.metadata.product, 'RW')
        self.assertEqual(data.data.shape, (900, 900))

    def test_get_lonlat_sphere(self):
        """
        Test if get_lonlat_sphere delivers the right shape of coordinates. Further, get_lonlat shall return a ValueError
        if the number of rows is not 900 or 1100.
        """
        lon, lat = readRadolan.get_lonlat_sphere(1100)
        self.assertEqual(lon.data.shape, (1100, 900))

        with self.assertRaises(ValueError):
            readRadolan.get_lonlat_sphere(1000)

    def test_get_lonlat(self):
        """
        get_lonlat is tested for the right lower left and upper right coordinates for different radolan variants.
        """
        # corners in radolanrx grid, VS <= 4
        sw_900 = [3.588929951, 46.95258041]
        ne_900 = [15.72075586, 54.74054767]
        lon, lat = readRadolan.get_lonlat(4, grid_variant='RadolanRX', corners=True)
        self.assertAlmostEqual(lon.data[1, 0], sw_900[0], places=5)
        self.assertAlmostEqual(lat.data[1, 0], sw_900[1], places=5)
        self.assertAlmostEqual(lon.data[0, 1], ne_900[0], places=5)
        self.assertAlmostEqual(lat.data[0, 1], ne_900[1], places=5)

        # corners in radolanrx grid, VS >= 5
        sw_900 = [3.604382997, 46.95361533]
        ne_900 = [15.69697166, 54.73806893]
        lon, lat = readRadolan.get_lonlat(5, num_rows=900, corners=True)
        self.assertAlmostEqual(lon.data[1, 0], sw_900[0], places=5)
        self.assertAlmostEqual(lat.data[1, 0], sw_900[1], places=5)
        self.assertAlmostEqual(lon.data[0, 1], ne_900[0], places=5)
        self.assertAlmostEqual(lat.data[0, 1], ne_900[1], places=5)

        # corners in de1200 grid, VS <= 4
        sw_1200 = [3.551921296, 45.69587048]
        ne_1200 = [18.76728172, 55.84848692]
        lon, lat = readRadolan.get_lonlat(4, grid_variant='DE1200', corners=True)
        self.assertAlmostEqual(lon.data[1, 0], sw_1200[0], places=5)
        self.assertAlmostEqual(lat.data[1, 0], sw_1200[1], places=5)
        self.assertAlmostEqual(lon.data[0, 1], ne_1200[0], places=5)
        self.assertAlmostEqual(lat.data[0, 1], ne_1200[1], places=5)

        # corners in de1200 grid, VS >= 5
        sw_1200 = [3.566994635, 45.69642538]
        ne_1200 = [18.73161645, 55.84543856]
        lon, lat = readRadolan.get_lonlat(5, num_rows=1200, corners=True)
        self.assertAlmostEqual(lon.data[1, 0], sw_1200[0], places=5)
        self.assertAlmostEqual(lat.data[1, 0], sw_1200[1], places=5)
        self.assertAlmostEqual(lon.data[0, 1], ne_1200[0], places=5)
        self.assertAlmostEqual(lat.data[0, 1], ne_1200[1], places=5)

        # corners in de4800 grid, VS <= 4
        sw_4800 = [3.551921296, 45.69587048]
        ne_4800 = [18.76728172, 55.84848692]
        lon, lat = readRadolan.get_lonlat(4, grid_variant='DE4800', corners=True)
        self.assertAlmostEqual(lon.data[1, 0], sw_4800[0], places=5)
        self.assertAlmostEqual(lat.data[1, 0], sw_4800[1], places=5)
        self.assertAlmostEqual(lon.data[0, 1], ne_4800[0], places=5)
        self.assertAlmostEqual(lat.data[0, 1], ne_4800[1], places=5)

        # corners in de4800 grid, VS >= 5
        sw_4800 = [3.566994635, 45.69642538]
        ne_4800 = [18.73161645, 55.84543856]
        lon, lat = readRadolan.get_lonlat(5, num_rows=4800, corners=True)
        self.assertAlmostEqual(lon.data[1, 0], sw_4800[0], places=5)
        self.assertAlmostEqual(lat.data[1, 0], sw_4800[1], places=5)
        self.assertAlmostEqual(lon.data[0, 1], ne_4800[0], places=5)
        self.assertAlmostEqual(lat.data[0, 1], ne_4800[1], places=5)

        # corners in europecomposite grid, VS <= 4
        sw_2400 = [-6.563068217, 38.80143975]
        ne_2400 = [24.24544273, 60.26405993]
        lon, lat = readRadolan.get_lonlat(4, grid_variant='EuropeComposite', corners=True)
        self.assertAlmostEqual(lon.data[1, 0], sw_2400[0], places=5)
        self.assertAlmostEqual(lat.data[1, 0], sw_2400[1], places=5)
        self.assertAlmostEqual(lon.data[0, 1], ne_2400[0], places=5)
        self.assertAlmostEqual(lat.data[0, 1], ne_2400[1], places=5)

        # corners in europecomposite grid, VS >= 5
        sw_2400 = [-6.526855245, 38.79723387]
        ne_2400 = [24.18150939, 60.25549810]
        lon, lat = readRadolan.get_lonlat(5, num_rows=2400, corners=True)
        self.assertAlmostEqual(lon.data[1, 0], sw_2400[0], places=5)
        self.assertAlmostEqual(lat.data[1, 0], sw_2400[1], places=5)
        self.assertAlmostEqual(lon.data[0, 1], ne_2400[0], places=5)
        self.assertAlmostEqual(lat.data[0, 1], ne_2400[1], places=5)
