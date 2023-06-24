#   @author: Michael Wagner
#   @organization: TU Dresden
#   @contact: michael.wagner@tu-dresden.de
import unittest

import numpy as np
import datetime as dt
import netCDF4 as nc

from met_entities.CosmoD2 import CosmoD2
from met_entities.CosmoD2EPS import CosmoD2EPS
from met_entities.GeoReferencedData import GeoReferencedData
from met_entities.RadolanRW import RadolanRW
from met_entities.RadvorRQ import RadvorRQ
from met_entities.RadolanRV import RadolanRV
from met_entities.IconD2 import IconD2
from met_entities.IconD2EPS import IconD2EPS
from met_entities.IconEU import IconEU
from met_entities.IconEUEPS import IconEUEPS
from met_entities.VariableDescription import RegridDescription
from met_entities.LonLatTime import LonLatTime
from met_entities.WeatherData import WeatherData
from met_entities.read_netcdf import read_netcdf


class TestRadolanRW(unittest.TestCase):
    def test_read_file(self):
        """
        Test of proper reading of data and the correct usage of scaling and filling.
        """
        # read with filename
        rd = RadolanRW()
        rd.read_file('data/raa01-rw_10000-2212141250-dwd---bin.gz')
        self.assertEqual(rd.gr_data.data.shape, (900, 900))
        self.assertAlmostEqual(rd.gr_data.data[829, 364], 10.1)

        # read with start datetime
        rd_datum = RadolanRW()
        rd_datum.read_file(start_datetime=dt.datetime(2022, 12, 14, 13), directory='data')
        self.assertAlmostEqual(rd.gr_data.data[829, 364], rd_datum.gr_data.data[829, 364])

        # read with start datetime including minutes
        rd_datum_min = RadolanRW()
        rd_datum_min.read_file(start_datetime=dt.datetime(2022, 12, 14, 12, 50), directory='data')
        self.assertAlmostEqual(rd.gr_data.data[829, 364], rd_datum_min.gr_data.data[829, 364])

        # read with scale factor
        scale_factor = 0.1
        rd_scaled = RadolanRW()
        rd_scaled.read_file('data/raa01-rw_10000-2212141250-dwd---bin.gz', scale_factor=scale_factor)
        self.assertEqual(rd.gr_data.data[270, 704] / scale_factor, rd_scaled.gr_data.data[270, 704])

        # read with fill value
        fill_value = -9999
        rd_fill_value = RadolanRW()
        rd_fill_value.read_file('data/raa01-rw_10000-2212141250-dwd---bin.gz', fill_value=fill_value)
        self.assertEqual(rd_fill_value.gr_data.data[0, 0], fill_value)

    def test_export_netcdf(self):
        """
        Test of export to netcdf functionality.
        """
        # export to netcdf
        rd = RadolanRW()
        rd.read_file('data/raa01-rw_10000-2212141150-dwd---bin.gz')
        rd.crop(lon_west=11.7, lon_east=15.2, lat_south=50.1, lat_north=51.8)
        nc_filename = 'data/RadolanRW_example_export.nc'
        rd.export_netcdf(nc_filename, data_format='f8')

        with nc.Dataset(nc_filename, mode='r') as nc_file:
            rd_from_nc = np.squeeze(nc_file.variables[rd.nc_desc.var_data])
        self.assertEqual(rd.gr_data.data.shape, rd_from_nc.shape)

        # append to existing netcdf
        rd_false_cropped = RadolanRW()
        rd_false_cropped.read_file('data/raa01-rw_10000-2212141250-dwd---bin.gz')
        rd_false_cropped.crop(lon_west=11.0, lon_east=15.2, lat_south=50.1, lat_north=51.8)
        with self.assertRaises(Exception):
            rd_false_cropped.export_netcdf_append(nc_filename[1:-2])
        with self.assertRaises(Exception):
            rd_false_cropped.export_netcdf_append(nc_filename)
        rd = RadolanRW()
        rd.read_file('data/raa01-rw_10000-2212141250-dwd---bin.gz')
        rd.crop(lon_west=11.7, lon_east=15.2, lat_south=50.1, lat_north=51.8)
        rd.export_netcdf_append(nc_filename)
        with nc.Dataset(nc_filename, mode='r') as nc_file:
            rd_from_nc = np.squeeze(nc_file.variables[rd.nc_desc.var_data])
        self.assertEqual(rd_from_nc.shape[0], 2)
        self.assertEqual(rd.gr_data.data.shape, rd_from_nc.shape[1:])


class TestRadvorRQ(unittest.TestCase):
    def test_read_file(self):
        """
        Test of proper reading of data and the correct usage of scaling and filling.
        """
        # read with start datetime
        start_datetime = dt.datetime(2022, 12, 14, 13)
        rd = RadvorRQ()
        rd.read_file(start_datetime=start_datetime, directory='data')
        self.assertEqual(rd.gr_data.data.shape, (900, 900, 3))
        self.assertAlmostEqual(rd.gr_data.data[819, 368, 0], 8.3)

        # read with scale factor
        scale_factor = 0.1
        rd_scaled = RadvorRQ()
        rd_scaled.read_file(start_datetime=start_datetime, directory='data', scale_factor=scale_factor)
        self.assertEqual(rd.gr_data.data[468, 651, 0] / scale_factor, rd_scaled.gr_data.data[468, 651, 0])

        # read with fill value
        fill_value = -9999
        rd_fill_value = RadvorRQ()
        rd_fill_value.read_file(start_datetime=start_datetime, directory='data', fill_value=fill_value)
        self.assertEqual(rd_fill_value.gr_data.data[0, 0, 0], fill_value)

    def test_export_netcdf(self):
        """
        Test of export to netcdf functionality.
        """
        # export to netcdf
        start_datetime = dt.datetime(2022, 12, 14, 13)
        rd = RadvorRQ()
        rd.read_file(start_datetime=start_datetime, directory='data')
        rd.crop(lon_west=11.7, lon_east=15.2, lat_south=50.1, lat_north=51.8)
        nc_filename = 'data/RadvorRQ_example_export.nc'
        rd.export_netcdf(nc_filename, data_format='f8')

        with nc.Dataset(nc_filename, mode='r') as nc_file:
            rd_from_nc = np.squeeze(nc_file.variables[rd.nc_desc.var_data])
        self.assertEqual(rd.gr_data.data.shape, rd_from_nc.shape)

        # append to existing netcdf
        rd_false_cropped = RadvorRQ()
        rd_false_cropped.read_file(start_datetime=start_datetime, directory='data')
        rd_false_cropped.crop(lon_west=11.0, lon_east=15.2, lat_south=50.1, lat_north=51.8)
        with self.assertRaises(Exception):
            rd_false_cropped.export_netcdf_append(nc_filename[1:-2])
        with self.assertRaises(Exception):
            rd_false_cropped.export_netcdf_append(nc_filename)
        rd = RadvorRQ()
        rd.read_file(start_datetime=start_datetime, directory='data')
        rd.crop(lon_west=11.7, lon_east=15.2, lat_south=50.1, lat_north=51.8)
        rd.export_netcdf_append(nc_filename)
        with nc.Dataset(nc_filename, mode='r') as nc_file:
            rd_from_nc = np.squeeze(nc_file.variables[rd.nc_desc.var_data])
        self.assertEqual(rd_from_nc.shape[0], 2)
        self.assertEqual(rd.gr_data.data.shape, rd_from_nc.shape[1:])


class TestRadolanRV(unittest.TestCase):
    def test_read_file(self):
        """
        Test of proper reading of data and the correct usage of scaling and filling.
        """
        # read with start datetime
        start_datetime = dt.datetime(2022, 12, 14, 13)
        rd = RadolanRV()
        rd.read_file(start_datetime=start_datetime, directory='data')
        self.assertEqual(rd.gr_data.data.shape, (1200, 1100, 25))
        self.assertAlmostEqual(rd.gr_data.data[970, 405, 0], 1.28)

        # read with filename
        rd_datum = RadolanRV()
        rd_datum.read_file('data/DE1200_RV2212141300.tar.bz2')
        self.assertAlmostEqual(rd_datum.gr_data.data[970, 405, 0], rd.gr_data.data[970, 405, 0])

        # read with scale factor
        scale_factor = 0.1
        rd_scaled = RadolanRV()
        rd_scaled.read_file(start_datetime=start_datetime, directory='data', scale_factor=scale_factor)
        self.assertEqual(rd.gr_data.data[468, 651, 0] / scale_factor, rd_scaled.gr_data.data[468, 651, 0])

        # read with fill value
        fill_value = -9999
        rd_fill_value = RadolanRV()
        rd_fill_value.read_file(start_datetime=start_datetime, directory='data', fill_value=fill_value)
        self.assertEqual(rd_fill_value.gr_data.data[0, 0, 0], fill_value)

    def test_export_netcdf(self):
        """
        Test of export to netcdf functionality.
        """
        # export to netcdf
        start_datetime = dt.datetime(2022, 12, 14, 13)
        rd = RadolanRV()
        rd.read_file(start_datetime=start_datetime, directory='data')
        rd.crop(lon_west=11.7, lon_east=15.2, lat_south=50.1, lat_north=51.8)
        nc_filename = 'data/RadolanRV_example_export.nc'

        rd.export_netcdf(nc_filename, data_format='f8')
        with nc.Dataset(nc_filename, mode='r') as nc_file:
            rd_from_nc = np.squeeze(nc_file.variables[rd.nc_desc.var_data])
        self.assertEqual(rd.gr_data.data.shape, rd_from_nc.shape)

        # append to existing netcdf
        rd_false_cropped = RadolanRV()
        rd_false_cropped.read_file(start_datetime=start_datetime, directory='data')
        rd_false_cropped.crop(lon_west=11.0, lon_east=15.2, lat_south=50.1, lat_north=51.8)
        with self.assertRaises(Exception):
            rd_false_cropped.export_netcdf_append(nc_filename[1:-2])
        with self.assertRaises(Exception):
            rd_false_cropped.export_netcdf_append(nc_filename)
        rd = RadolanRV()
        rd.read_file(start_datetime=dt.datetime(2022, 12, 14, 13, 5), directory='data')
        rd.crop(lon_west=11.7, lon_east=15.2, lat_south=50.1, lat_north=51.8)
        rd.export_netcdf_append(nc_filename)
        with nc.Dataset(nc_filename, mode='r') as nc_file:
            rd_from_nc = np.squeeze(nc_file.variables[rd.nc_desc.var_data])
        self.assertEqual(rd_from_nc.shape[0], 2)
        self.assertEqual(rd.gr_data.data.shape, rd_from_nc.shape[1:])


class TestIconD2(unittest.TestCase):
    def test_read_file(self):
        """
        Test of proper reading of data and the correct usage of scaling and filling.
        """
        # read with start datetime
        start_datetime = dt.datetime(2022, 12, 14, 12)
        icond2 = IconD2()
        icond2.read_file(start_datetime=start_datetime, forecast_hours=2, directory='data')
        self.assertEqual(icond2.gr_data.data.shape, (542040, 12))
        self.assertAlmostEqual(icond2.gr_data.data[30000, 2], 0.1221313, places=7)

        # read with scale factor
        scale_factor = 0.1
        icond2_scaled = IconD2()
        icond2_scaled.read_file(start_datetime=start_datetime, directory='data', forecast_hours=0,
                                scale_factor=scale_factor)
        self.assertAlmostEqual(icond2_scaled.gr_data.data[30000, 2] * scale_factor, icond2.gr_data.data[30000, 2])

        # read with fill value
        fill_value = -9999
        icond2_fill_value = IconD2()
        icond2_fill_value.read_file(start_datetime=start_datetime, directory='data', forecast_hours=0,
                                    fill_value=fill_value)
        self.assertEqual(icond2_fill_value.gr_data.data[0, 0], fill_value)

    def test_export_netcdf(self):
        """
        Test of export to netcdf functionality.
        """
        # export to netcdf
        start_datetime = dt.datetime(2022, 12, 14, 12)
        icond2 = IconD2()
        icond2.read_file(start_datetime=start_datetime, directory='data', forecast_hours=2)
        icond2.crop(lon_west=11.7, lon_east=15.2, lat_south=50.1, lat_north=51.8)
        nc_filename = 'data/IconD2_example_export.nc'
        icond2.export_netcdf(nc_filename, data_format='f8')

        with nc.Dataset(nc_filename, mode='r') as nc_file:
            icond2_from_nc = np.squeeze(nc_file.variables[icond2.nc_desc.var_data])
        self.assertEqual(icond2.gr_data.data.shape, icond2_from_nc.shape)

        # append to existing netcdf
        icond2_false_cropped = IconD2()
        icond2_false_cropped.read_file(start_datetime=start_datetime, directory='data', forecast_hours=2)
        icond2_false_cropped.crop(lon_west=11.0, lon_east=15.2, lat_south=50.1, lat_north=51.8)
        with self.assertRaises(Exception):
            icond2_false_cropped.export_netcdf_append(nc_filename[1:-2])
        with self.assertRaises(Exception):
            icond2_false_cropped.export_netcdf_append(nc_filename)
        icond2 = IconD2()
        icond2.read_file(start_datetime=start_datetime, directory='data', forecast_hours=2)
        icond2.crop(lon_west=11.7, lon_east=15.2, lat_south=50.1, lat_north=51.8)
        icond2.export_netcdf_append(nc_filename)
        with nc.Dataset(nc_filename, mode='r') as nc_file:
            icond2_from_nc = np.squeeze(nc_file.variables[icond2.nc_desc.var_data])
        self.assertEqual(icond2_from_nc.shape[0], 2)
        self.assertEqual(icond2.gr_data.data.shape, icond2_from_nc.shape[1:])


class TestIconD2EPS(unittest.TestCase):
    def test_read_file(self):
        """
        Test of proper reading of data and the correct usage of scaling and filling.
        """
        # read with start datetime
        start_datetime = dt.datetime(2022, 12, 14, 15)
        icond2eps = IconD2EPS()
        icond2eps.read_file(start_datetime=start_datetime, directory='data', forecast_hours=0)
        self.assertEqual(len(icond2eps.gr_data), 20)
        self.assertEqual(icond2eps.gr_data[0].data.shape, (542040, 4))
        self.assertAlmostEqual(icond2eps.gr_data[0].data[30000, 3], 0.1293945, places=6)

        # read single ensemble members
        icond2eps_member = IconD2EPS()
        icond2eps_member.read_file(start_datetime=start_datetime, directory='data', forecast_hours=0,
                                   eps_member=[4, 7, 15])
        self.assertEqual(len(icond2eps_member.gr_data), 3)
        self.assertEqual(icond2eps_member.gr_data[0].data_description.eps_note, 'eps member 4')
        self.assertEqual(icond2eps_member.gr_data[1].data_description.eps_note, 'eps member 7')
        self.assertEqual(icond2eps_member.gr_data[2].data_description.eps_note, 'eps member 15')

        # read with scale_factor
        scale_factor = 0.1
        icond2eps_scaled = IconD2EPS()
        icond2eps_scaled.read_file(start_datetime=start_datetime, directory='data', forecast_hours=0,
                                   scale_factor=scale_factor)
        self.assertAlmostEqual(icond2eps_scaled.gr_data[0].data[30000, 2] * scale_factor, icond2eps.gr_data[0]
                               .data[30000, 2])

        # read with fill value
        fill_value = -9999
        icond2eps_fill_value = IconD2EPS()
        icond2eps_fill_value.read_file(start_datetime=start_datetime, directory='data', forecast_hours=0,
                                       fill_value=fill_value)
        self.assertEqual(icond2eps_fill_value.gr_data[0].data[0, 0], fill_value)

    def test_export_netcdf(self):
        """
        Test of export to netcdf functionality.
        """
        # export to netcdf
        start_datetime = dt.datetime(2022, 12, 14, 15)
        icond2eps = IconD2EPS()
        icond2eps.read_file(start_datetime=start_datetime, directory='data', forecast_hours=0)
        icond2eps.crop(lon_west=11.7, lon_east=15.2, lat_south=50.1, lat_north=51.8)
        nc_filename = 'data/IconD2EPS_example_export.nc'
        icond2eps.export_netcdf(nc_filename, data_format='f8', version='integrated')

        with nc.Dataset(nc_filename, mode='r') as nc_file:
            icond2eps_from_nc = np.squeeze(nc_file.variables[icond2eps.nc_desc.var_data])
        self.assertEqual(icond2eps.gr_data[0].data.shape, icond2eps_from_nc.shape[0:-1])
        self.assertEqual(icond2eps_from_nc.shape[-1], 20)

        # append to existing netcdf
        icond2eps_false_cropped = IconD2EPS()
        icond2eps_false_cropped.read_file(start_datetime=start_datetime, directory='data', forecast_hours=0)
        icond2eps_false_cropped.crop(lon_west=11.0, lon_east=15.2, lat_south=50.1, lat_north=51.8)
        with self.assertRaises(Exception):
            icond2eps_false_cropped.export_netcdf_append(nc_filename[1:-2])
        with self.assertRaises(Exception):
            icond2eps_false_cropped.export_netcdf_append(nc_filename)
        icond2eps = IconD2EPS()
        icond2eps.read_file(start_datetime=start_datetime, directory='data', forecast_hours=0)
        icond2eps.crop(lon_west=11.7, lon_east=15.2, lat_south=50.1, lat_north=51.8)
        icond2eps.export_netcdf_append(nc_filename)
        with nc.Dataset(nc_filename, mode='r') as nc_file:
            icond2eps_from_nc = np.squeeze(nc_file.variables[icond2eps.nc_desc.var_data])
        self.assertEqual(icond2eps_from_nc.shape[0], 2)
        self.assertEqual(icond2eps.gr_data[0].data.shape, icond2eps_from_nc.shape[1:-1])
        self.assertEqual(icond2eps_from_nc.shape[-1], 20)


class TestIconEU(unittest.TestCase):
    def test_read_file(self):
        """
        Test of proper reading of data and the correct usage of scaling and filling.
        """
        # read with start datetime
        start_datetime = dt.datetime(2023, 3, 23, 0)
        iconeu = IconEU()
        iconeu.read_file(start_datetime=start_datetime, forecast_hours=1, directory='data')
        self.assertEqual(iconeu.gr_data.data.shape, (657, 1377, 2))
        self.assertAlmostEqual(iconeu.gr_data.data[574, 1316, 1], 14.8027344, places=7)

        # read with scale factor
        scale_factor = 0.1
        iconeu_scaled = IconEU()
        iconeu_scaled.read_file(start_datetime=start_datetime, directory='data', forecast_hours=1,
                                scale_factor=scale_factor)
        self.assertAlmostEqual(iconeu_scaled.gr_data.data[574, 1316, 1] * scale_factor, iconeu.gr_data.data[574, 1316, 1])

    def test_export_netcdf(self):
        """
        Test of export to netcdf functionality.
        """
        # export to netcdf
        start_datetime = dt.datetime(2023, 3, 23, 0)
        iconeu = IconEU()
        iconeu.read_file(start_datetime=start_datetime, directory='data', forecast_hours=1)
        iconeu.crop(lon_west=11.7, lon_east=15.2, lat_south=50.1, lat_north=51.8)
        nc_filename = 'data/IconEU_example_export.nc'
        iconeu.export_netcdf(nc_filename, data_format='f8')

        with nc.Dataset(nc_filename, mode='r') as nc_file:
            iconeu_from_nc = np.squeeze(nc_file.variables[iconeu.nc_desc.var_data])
        self.assertEqual(iconeu.gr_data.data.shape, iconeu_from_nc.shape)

        # append to existing netcdf
        iconeu_false_cropped = IconEU()
        iconeu_false_cropped.read_file(start_datetime=start_datetime, directory='data', forecast_hours=1)
        iconeu_false_cropped.crop(lon_west=11.0, lon_east=15.2, lat_south=50.1, lat_north=51.8)
        with self.assertRaises(Exception):
            iconeu_false_cropped.export_netcdf_append(nc_filename[1:-2])
        with self.assertRaises(Exception):
            iconeu_false_cropped.export_netcdf_append(nc_filename)
        iconeu = IconEU()
        iconeu.read_file(start_datetime=start_datetime, directory='data', forecast_hours=1)
        iconeu.crop(lon_west=11.7, lon_east=15.2, lat_south=50.1, lat_north=51.8)
        iconeu.export_netcdf_append(nc_filename)
        with nc.Dataset(nc_filename, mode='r') as nc_file:
            iconeu_from_nc = np.squeeze(nc_file.variables[iconeu.nc_desc.var_data])
        self.assertEqual(iconeu_from_nc.shape[0], 2)
        self.assertEqual(iconeu.gr_data.data.shape, iconeu_from_nc.shape[1:])


class TestIconEUEPS(unittest.TestCase):
    def test_read_file(self):
        """
        Test of proper reading of data and the correct usage of scaling and filling.
        """
        # read with start datetime
        start_datetime = dt.datetime(2023, 3, 23, 0)
        iconeueps = IconEUEPS()
        iconeueps.read_file(start_datetime=start_datetime, directory='data', forecast_hours=1)
        self.assertEqual(len(iconeueps.gr_data), 40)
        self.assertEqual(iconeueps.gr_data[0].data.shape, (164984, 2))
        self.assertAlmostEqual(iconeueps.gr_data[0].data[44681, 1], 8.532471, places=6)

        # read single ensemble members
        iconeueps_member = IconEUEPS()
        iconeueps_member.read_file(start_datetime=start_datetime, directory='data', forecast_hours=1,
                                   eps_member=[4, 7, 25])
        self.assertEqual(len(iconeueps_member.gr_data), 3)
        self.assertEqual(iconeueps_member.gr_data[0].data_description.eps_note, 'eps member 4')
        self.assertEqual(iconeueps_member.gr_data[1].data_description.eps_note, 'eps member 7')
        self.assertEqual(iconeueps_member.gr_data[2].data_description.eps_note, 'eps member 25')

        # read with scale_factor
        scale_factor = 0.1
        iconeueps_scaled = IconEUEPS()
        iconeueps_scaled.read_file(start_datetime=start_datetime, directory='data', forecast_hours=1,
                                   scale_factor=scale_factor)
        self.assertAlmostEqual(iconeueps_scaled.gr_data[0].data[44681, 1] * scale_factor, iconeueps.gr_data[0]
                               .data[44681, 1])

    def test_export_netcdf(self):
        """
        Test of export to netcdf functionality.
        """
        # export to netcdf
        start_datetime = dt.datetime(2023, 3, 23, 0)
        iconeueps = IconEUEPS()
        iconeueps.read_file(start_datetime=start_datetime, directory='data', forecast_hours=1)
        iconeueps.crop(lon_west=11.7, lon_east=15.2, lat_south=50.1, lat_north=51.8)
        nc_filename = 'data/IconEUEPS_example_export.nc'
        iconeueps.export_netcdf(nc_filename, data_format='f8', version='integrated')

        with nc.Dataset(nc_filename, mode='r') as nc_file:
            iconeueps_from_nc = np.squeeze(nc_file.variables[iconeueps.nc_desc.var_data])
        self.assertEqual(iconeueps.gr_data[0].data.shape, iconeueps_from_nc.shape[0:-1])
        self.assertEqual(iconeueps_from_nc.shape[-1], 40)

        # append to existing netcdf
        iconeueps_false_cropped = IconEUEPS()
        iconeueps_false_cropped.read_file(start_datetime=start_datetime, directory='data', forecast_hours=1)
        iconeueps_false_cropped.crop(lon_west=11.0, lon_east=15.2, lat_south=50.1, lat_north=51.8)
        with self.assertRaises(Exception):
            iconeueps_false_cropped.export_netcdf_append(nc_filename[1:-2])
        with self.assertRaises(Exception):
            iconeueps_false_cropped.export_netcdf_append(nc_filename)
        iconeueps = IconEUEPS()
        iconeueps.read_file(start_datetime=start_datetime, directory='data', forecast_hours=1)
        iconeueps.crop(lon_west=11.7, lon_east=15.2, lat_south=50.1, lat_north=51.8)
        iconeueps.export_netcdf_append(nc_filename)
        with nc.Dataset(nc_filename, mode='r') as nc_file:
            iconeueps_from_nc = np.squeeze(nc_file.variables[iconeueps.nc_desc.var_data])
        self.assertEqual(iconeueps_from_nc.shape[0], 2)
        self.assertEqual(iconeueps.gr_data[0].data.shape, iconeueps_from_nc.shape[1:-1])
        self.assertEqual(iconeueps_from_nc.shape[-1], 40)


class TestCosmoD2(unittest.TestCase):
    def test_read_file(self):
        """
        Test of proper reading of data and the correct usage of scaling and filling.
        """
        # read with start datetime
        start_datetime = dt.datetime(2021, 1, 17, 0)
        cosmod2 = CosmoD2()
        cosmod2.read_file(start_datetime=start_datetime, forecast_hours=1, directory='data')
        self.assertEqual(cosmod2.gr_data.data.shape, (716, 651, 8))
        self.assertAlmostEqual(cosmod2.gr_data.data[86, 75, 3], 5.419556, places=5)

        # read with scale factor
        scale_factor = 0.1
        cosmod2_scaled = CosmoD2()
        cosmod2_scaled.read_file(start_datetime=start_datetime, directory='data', forecast_hours=0,
                                 scale_factor=scale_factor)
        self.assertAlmostEqual(cosmod2_scaled.gr_data.data[86, 75, 3] * scale_factor, cosmod2.gr_data.data[86, 75, 3])

    def test_export_netcdf(self):
        """
        Test of export to netcdf functionality.
        """
        # export to netcdf
        start_datetime = dt.datetime(2021, 1, 17, 0)
        cosmod2 = CosmoD2()
        cosmod2.read_file(start_datetime=start_datetime, directory='data', forecast_hours=1)
        cosmod2.crop(lon_west=11.7, lon_east=15.2, lat_south=50.1, lat_north=51.8)
        nc_filename = 'data/CosmoD2_example_export.nc'
        cosmod2.export_netcdf(nc_filename, data_format='f8')

        with nc.Dataset(nc_filename, mode='r') as nc_file:
            cosmod2_from_nc = np.squeeze(nc_file.variables[cosmod2.nc_desc.var_data])
        self.assertEqual(cosmod2.gr_data.data.shape, cosmod2_from_nc.shape)

        # append to existing netcdf
        cosmod2_false_cropped = CosmoD2()
        cosmod2_false_cropped.read_file(start_datetime=start_datetime, directory='data', forecast_hours=1)
        cosmod2_false_cropped.crop(lon_west=11.0, lon_east=15.2, lat_south=50.1, lat_north=51.8)
        with self.assertRaises(Exception):
            cosmod2_false_cropped.export_netcdf_append(nc_filename[1:-2])
        with self.assertRaises(Exception):
            cosmod2_false_cropped.export_netcdf_append(nc_filename)
        cosmod2 = CosmoD2()
        cosmod2.read_file(start_datetime=start_datetime, directory='data', forecast_hours=1)
        cosmod2.crop(lon_west=11.7, lon_east=15.2, lat_south=50.1, lat_north=51.8)
        cosmod2.export_netcdf_append(nc_filename)
        with nc.Dataset(nc_filename, mode='r') as nc_file:
            cosmod2_from_nc = np.squeeze(nc_file.variables[cosmod2.nc_desc.var_data])
        self.assertEqual(cosmod2_from_nc.shape[0], 2)
        self.assertEqual(cosmod2.gr_data.data.shape, cosmod2_from_nc.shape[1:])


class TestCosmoD2EPS(unittest.TestCase):
    def test_read_file(self):
        """
        Test of proper reading of data and the correct usage of scaling and filling.
        """
        # read with start datetime
        start_datetime = dt.datetime(2021, 1, 17, 0)
        cosmod2eps = CosmoD2EPS()
        cosmod2eps.read_file(start_datetime=start_datetime, directory='data', forecast_hours=0)
        self.assertEqual(len(cosmod2eps.gr_data), 20)
        self.assertEqual(cosmod2eps.gr_data[0].data.shape, (716, 651, 4))
        self.assertAlmostEqual(cosmod2eps.gr_data[0].data[86, 84, 3], 3.647949, places=5)

        # read single ensemble members
        cosmod2eps_member = CosmoD2EPS()
        cosmod2eps_member.read_file(start_datetime=start_datetime, directory='data', forecast_hours=0,
                                   eps_member=[5, 10, 19])
        self.assertEqual(len(cosmod2eps_member.gr_data), 3)
        self.assertEqual(cosmod2eps_member.gr_data[0].data_description.eps_note, 'eps member 5')
        self.assertEqual(cosmod2eps_member.gr_data[1].data_description.eps_note, 'eps member 10')
        self.assertEqual(cosmod2eps_member.gr_data[2].data_description.eps_note, 'eps member 19')

        # read with scale_factor
        scale_factor = 0.1
        cosmod2eps_scaled = CosmoD2EPS()
        cosmod2eps_scaled.read_file(start_datetime=start_datetime, directory='data', forecast_hours=0,
                                   scale_factor=scale_factor)
        self.assertAlmostEqual(cosmod2eps_scaled.gr_data[0].data[86, 84, 3] * scale_factor, cosmod2eps.gr_data[0]
                               .data[86, 84, 3])

    def test_export_netcdf(self):
        """
        Test of export to netcdf functionality.
        """
        # export to netcdf
        start_datetime = dt.datetime(2021, 1, 17, 0)
        cosmod2eps = CosmoD2EPS()
        cosmod2eps.read_file(start_datetime=start_datetime, directory='data', forecast_hours=0)
        cosmod2eps.crop(lon_west=11.7, lon_east=15.2, lat_south=50.1, lat_north=51.8)
        nc_filename = 'data/CosmoD2EPS_example_export.nc'
        cosmod2eps.export_netcdf(nc_filename, data_format='f8', version='integrated')

        with nc.Dataset(nc_filename, mode='r') as nc_file:
            cosmod2eps_from_nc = np.squeeze(nc_file.variables[cosmod2eps.nc_desc.var_data])
        self.assertEqual(cosmod2eps.gr_data[0].data.shape, cosmod2eps_from_nc.shape[0:-1])
        self.assertEqual(cosmod2eps_from_nc.shape[-1], 20)

        # append to existing netcdf
        cosmod2eps_false_cropped = CosmoD2EPS()
        cosmod2eps_false_cropped.read_file(start_datetime=start_datetime, directory='data', forecast_hours=0)
        cosmod2eps_false_cropped.crop(lon_west=11.0, lon_east=15.2, lat_south=50.1, lat_north=51.8)
        with self.assertRaises(Exception):
            cosmod2eps_false_cropped.export_netcdf_append(nc_filename[1:-2])
        with self.assertRaises(Exception):
            cosmod2eps_false_cropped.export_netcdf_append(nc_filename)
        cosmod2eps = CosmoD2EPS()
        cosmod2eps.read_file(start_datetime=start_datetime, directory='data', forecast_hours=0)
        cosmod2eps.crop(lon_west=11.7, lon_east=15.2, lat_south=50.1, lat_north=51.8)
        cosmod2eps.export_netcdf_append(nc_filename)
        with nc.Dataset(nc_filename, mode='r') as nc_file:
            cosmod2eps_from_nc = np.squeeze(nc_file.variables[cosmod2eps.nc_desc.var_data])
        self.assertEqual(cosmod2eps_from_nc.shape[0], 2)
        self.assertEqual(cosmod2eps.gr_data[0].data.shape, cosmod2eps_from_nc.shape[1:-1])
        self.assertEqual(cosmod2eps_from_nc.shape[-1], 20)


class TestWeatherData(unittest.TestCase):
    lon = np.arange(4, 15, .3)
    lat = np.arange(54, 48, -.3)
    lon_target, lat_target = np.meshgrid(lon, lat)
    lon_target = LonLatTime(data=lon_target)
    lat_target = LonLatTime(data=lat_target)
    regrid_description = {
        'radolanrw': RegridDescription(lon_target=lon_target, lat_target=lat_target,
                                       file_nearest='data/radrw_regridding.npz'),
        'radvorrq': RegridDescription(lon_target=lon_target, lat_target=lat_target,
                                      file_nearest='data/radrq_regridding.npz'),
        'radolanrv': RegridDescription(lon_target=lon_target, lat_target=lat_target,
                                       file_nearest='data/radrv_regridding.npz'),
        'icond2': RegridDescription(lon_target=lon_target, lat_target=lat_target,
                                    file_nearest='data/icond2_regridding.npz'),
        'icond2eps': RegridDescription(lon_target=lon_target, lat_target=lat_target,
                                       file_nearest='data/icond2eps_regridding.npz')}

    def test_collect_radolanrw(self):
        """
        Test of collecting RadolanRW data with harmonizing features.
        """
        fill_value = np.nan  # for easier content comparison with division
        time_start = dt.datetime(2022, 12, 14, 12, tzinfo=dt.timezone.utc)
        time_now = dt.datetime(2022, 12, 14, 13, tzinfo=dt.timezone.utc)

        # direct loading for content comparison
        rd_1 = RadolanRW()
        rd_1.read_file(start_datetime=time_start, directory='data', fill_value=fill_value)
        rd_1.regrid(regrid_description=self.regrid_description['radolanrw'])

        rd_2 = RadolanRW()
        rd_2.read_file(start_datetime=time_now, directory='data', fill_value=fill_value)
        rd_2.regrid(regrid_description=self.regrid_description['radolanrw'])

        # test 5 min resolution
        wd = WeatherData(time_now=time_now, delta_t=dt.timedelta(minutes=5), fill_value=fill_value,
                         regrid_description=self.regrid_description)
        wd.collect_radolanrw(time_start=time_start, directory='data')
        self.assertEqual(len(wd.data), 11)
        self.assertEqual(wd.data[0].data.shape, self.lon_target.data.shape)
        self.assertEqual(wd.time_now, time_now - dt.timedelta(minutes=10))  # changed time_now due to representable xx:50 observation time of RadolanRW products
        self.assertTrue(np_all_ignore_nan(wd.data[0].data, rd_2.gr_data.data / 12))  # 5 min must 1/12 from the original
        self.assertTrue(np_all_ignore_nan(wd.data[-1].data, rd_2.gr_data.data / 12))

        # test 15 min resolution
        wd = WeatherData(time_now=time_now, delta_t=dt.timedelta(minutes=15), fill_value=fill_value,
                         regrid_description=self.regrid_description)
        wd.collect_radolanrw(time_start=time_start, directory='data')
        self.assertEqual(len(wd.data), 5)
        self.assertTrue(np_all_ignore_nan(wd.data[0].data, rd_1.gr_data.data / 4))  # 15 min must be a quarter from the original
        self.assertTrue(np_all_ignore_nan(wd.data[1].data, rd_2.gr_data.data / 4))
        self.assertTrue(np_all_ignore_nan(wd.data[2].data, rd_2.gr_data.data / 4))
        self.assertTrue(np_all_ignore_nan(wd.data[3].data, rd_2.gr_data.data / 4))
        self.assertTrue(np_all_ignore_nan(wd.data[4].data, rd_2.gr_data.data / 4))

        # test 60 min resolution
        wd = WeatherData(time_now=time_now, delta_t=dt.timedelta(minutes=60), fill_value=fill_value,
                         regrid_description=self.regrid_description)
        wd.collect_radolanrw(time_start=time_start, directory='data')
        self.assertEqual(len(wd.data), 2)
        self.assertTrue(np_all_ignore_nan(wd.data[0].data, rd_1.gr_data.data))
        self.assertTrue(np_all_ignore_nan(wd.data[1].data, rd_2.gr_data.data))

        # test 120 min resolution
        wd = WeatherData(time_now=time_now, delta_t=dt.timedelta(minutes=120), fill_value=fill_value,
                         regrid_description=self.regrid_description)
        with self.assertRaises(Exception):
            wd.collect_radolanrw(time_start=time_start, directory='data')
        wd.collect_radolanrw(time_start=dt.datetime(2022, 12, 14, 11, tzinfo=dt.timezone.utc), directory='data')
        self.assertEqual(len(wd.data), 2)
        self.assertTrue(np_all_ignore_nan(wd.data[1].data, rd_1.gr_data.data + rd_2.gr_data.data))


class TestReadNetcdf(unittest.TestCase):
    def test_read_netcdf(self):
        """
        Test of reading in netcdf functionality. Read in a netcdf file previously exported by a MetEntities subclass
        export_netcdf method. The resulting class shall have all necessary entries and the data must not change.
        """
        # load RadolanRW data and export to netcdf
        rd = RadolanRW()
        rd.read_file(start_datetime=dt.datetime(2022, 12, 14, 11, 50), directory='data', scale_factor=0.1,
                     fill_value=-1)
        rd.crop(lon_west=11.7, lon_east=15.2, lat_south=50.1, lat_north=51.8)
        nc_filename = 'data/RadolanRW_example_export2.nc'
        rd.export_netcdf(nc_filename, data_format='f8')

        # read netcdf file and build all variables
        rd_in = read_netcdf(nc_filename)
        self.assertTrue(type(rd) is type(rd_in))
        self.assertEqual(rd.time_value.data_description.units, rd_in.time_value.data_description.units)
        self.assertEqual(rd.gr_data.data_description.scale_factor, rd_in.gr_data.data_description.scale_factor)
        self.assertEqual(rd.gr_data.data_description.time_note, rd_in.gr_data.data_description.time_note)
        self.assertEqual(rd.gr_data.data_description.units, rd_in.gr_data.data_description.units)
        self.assertTrue(np.all(rd.gr_data.data == rd_in.gr_data.data))

        # load RadolanRW data and export to netcdf - scale_undo
        rd = RadolanRW()
        rd.read_file(start_datetime=dt.datetime(2022, 12, 14, 11, 50), directory='data', scale_factor=0.1,
                     fill_value=-1)
        rd.crop(lon_west=11.7, lon_east=15.2, lat_south=50.1, lat_north=51.8)
        nc_filename = 'data/RadolanRW_example_export2.nc'
        rd.export_netcdf(nc_filename, data_format='f8', scale_undo=True)

        # read netcdf file and build all variables
        rd_in = read_netcdf(nc_filename)
        self.assertTrue(type(rd) is type(rd_in))
        self.assertEqual(rd.time_value.data_description.units, rd_in.time_value.data_description.units)
        self.assertEqual(rd.gr_data.data_description.time_note, rd_in.gr_data.data_description.time_note)
        idx_fill_value = rd.gr_data.data == rd.gr_data.data_description.fill_value
        data = rd.gr_data.data * rd.gr_data.data_description.scale_factor
        data[idx_fill_value] = rd.gr_data.data_description.fill_value
        self.assertTrue(np.all(data == rd_in.gr_data.data))

        # load RadolanRW data and export to netcdf - scale_undo and scale_factor_nc
        rd = RadolanRW()
        rd.read_file(start_datetime=dt.datetime(2022, 12, 14, 11, 50), directory='data', scale_factor=0.1,
                     fill_value=-1)
        rd.crop(lon_west=11.7, lon_east=15.2, lat_south=50.1, lat_north=51.8)
        nc_filename = 'data/RadolanRW_example_export2.nc'
        rd.export_netcdf(nc_filename, data_format='f8', scale_undo=True, scale_factor_nc=0.1)

        # read netcdf file and build all variables
        rd_in = read_netcdf(nc_filename)
        self.assertTrue(type(rd) is type(rd_in))
        self.assertEqual(rd.time_value.data_description.units, rd_in.time_value.data_description.units)
        self.assertEqual(rd.gr_data.data_description.scale_factor, rd_in.gr_data.data_description.scale_factor)
        self.assertEqual(rd.gr_data.data_description.time_note, rd_in.gr_data.data_description.time_note)
        self.assertEqual(rd.gr_data.data_description.units, rd_in.gr_data.data_description.units)
        self.assertTrue(np.all(rd.gr_data.data == rd_in.gr_data.data))


class TestGeoReferencedData(unittest.TestCase):
    def test_regrid_idw(self):
        """
        Test of regridding an instance of GeoReferencedData. The final shape must meet the target coordinates.
        """
        rd = RadolanRW()
        rd.read_file(start_datetime=dt.datetime(2022, 12, 14, 13, 0), directory='data')
        # transfer RadolanRW class gr_data to GeoReferencedData for testing low level
        grd = GeoReferencedData(lon=rd.gr_data.lon, lat=rd.gr_data.lat, data=rd.gr_data.data,
                                data_description=rd.gr_data.data_description)
        lon = np.arange(10, 12, 1)
        lat = np.arange(50, 48, -1)
        lon_target, lat_target = np.meshgrid(lon, lat)
        grd.regrid_idw(lon_target=lon_target, lat_target=lat_target)
        self.assertEqual(grd.data.shape, (lat.size, lon.size))
        self.assertEqual(grd.lon.data.shape, (lat.size, lon.size))
        self.assertEqual(grd.lat.data.shape, (lat.size, lon.size))

    def test_find_nearest(self):
        """
        Test of finding the nearest points from source grid to target grid and the according lengths. The values are
        compared with values found on the map.
        """
        rd = RadolanRW()
        rd.read_file(start_datetime=dt.datetime(2022, 12, 14, 13, 0), directory='data')
        grd = GeoReferencedData(lon=rd.gr_data.lon, lat=rd.gr_data.lat, data=rd.gr_data.data,
                                data_description=rd.gr_data.data_description)
        lon = 13
        lat = 51
        idx_nearest_via_map = [444, 743]
        length_nearest_direct = np.sqrt((grd.lon.data[444, 743] - 13)**2 + (grd.lat.data[444, 743] - 51)**2)
        lon_target, lat_target = np.meshgrid(lon, lat)
        idx_nearest, length_nearest = grd.find_nearest(lon_target=lon_target, lat_target=lat_target, neighbors=1)
        self.assertEqual(idx_nearest[0, 0, 0, 0], idx_nearest_via_map[0])
        self.assertEqual(idx_nearest[0, 0, 1, 0], idx_nearest_via_map[1])
        self.assertAlmostEqual(length_nearest_direct, length_nearest[0, 0, 0], places=5)

    def test_crop(self):
        """
        Test the crop functionality for proper exception raising and the correct size of a cropped dataset with cut
        first and last column/row.
        """
        rd = RadolanRW()
        rd.read_file(start_datetime=dt.datetime(2022, 12, 14, 13, 0), directory='data')
        grd = GeoReferencedData(lon=rd.gr_data.lon, lat=rd.gr_data.lat, data=rd.gr_data.data,
                                data_description=rd.gr_data.data_description)
        lon_west = 11.7
        lon_east = 15.2
        lat_south = 50.1
        lat_north = 51.8
        idx_west = 635
        idx_east = 899
        idx_south = 556
        idx_north = 324
        with self.assertRaises(Exception):
            grd.crop()
        with self.assertRaises(Exception):
            # changed north and south latitudes
            grd.crop(lon_west=lon_west, lon_east=lon_east, lat_south=lat_north, lat_north=lat_south)
        with self.assertRaises(Exception):
            # changed west and east indexes
            grd.crop(idx_west=idx_east, idx_east=idx_west, idx_south=idx_south, idx_north=idx_north)
        # whole matrix except first and last columns must be returned if the target bounding boxes west end is between
        # 2. and 3. lon min column-wise and the east end is between 898. and 899. max lon column-wise
        grd.crop(lon_west=3.61, lon_east=14.598, lat_south=47.21, lat_north=54.57)
        self.assertEqual(grd.data.shape, (898, 898))


def np_all_ignore_nan(mat0, mat1):
    """
    Compare two matrices while ignoring nan entries.

    :param mat0: first matrix to be compared
    :type mat0: numpy.ndarray
    :param mat1: second matrix to be compared
    :type mat1: numpy.ndarray
    :return: True if everything is the same and False if something differs
    :rtype: bool
    """
    mask = ~(np.isnan(mat0) | np.isnan(mat1))  # mask nan data

    return np.all(mat0[mask] == mat1[mask])
