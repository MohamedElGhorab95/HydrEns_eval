#   @author: Michael Wagner
#   @organization: TU Dresden
#   @contact: michael.wagner@tu-dresden.de
import unittest
import numpy as np
import copy

from aux_tools import grib2_tools as gt
from aux_tools import scaling_tools as st
from met_entities.GeoReferencedData import GeoReferencedData
from met_entities.VariableDescription import DataDescription

class TestGrib2Tools(unittest.TestCase):
    def test_get_keys(self):
        """
        Pick the name key in a grib file message and test its value.
        """
        keys = gt.get_keys('data/icon-d2_germany_icosahedral_single-level_2022121412_000_2d_tot_prec.grib2.bz2')
        self.assertEqual(keys['name'], 'Total Precipitation')

    def test_accum_to_instantaneous_flattened(self):
        """
        Test proper de-accumulation function. It shall work with multiple data. If a list with only one element is
        given, it shall return the same list reduced by one (list-) dimension.
        """
        data0 = DataExample(data=np.array([0., 1, 2, 3, 4]))
        data1 = DataExample(data=np.array([0., 1, 2, 3, 4]) + 1)
        data2 = DataExample(data=np.array([0., 1, 2, 3, 4]) + 2)
        data3 = DataExample(data=np.array([0., 1, 2, 3, 4]) + 4)

        # case 1 with only one element
        data_list1 = copy.deepcopy([[data0]])
        res_list1 = gt.accum_to_instantaneous_flattened(data_list1, np.nan)
        self.assertTrue(np.all(res_list1[0].data == data0.data))

        # case 2 with multiple elements
        data_list2 = copy.deepcopy([[data0, data1, data2, data3]])
        res_list2 = gt.accum_to_instantaneous_flattened(data_list2, np.nan)

        self.assertTrue(np.all(res_list2[0].data == data0.data))
        self.assertTrue(np.all(res_list2[1].data == data1.data - data0.data))
        self.assertTrue(np.all(res_list2[2].data == data2.data - data1.data))
        self.assertTrue(np.all(res_list2[3].data == data3.data - data2.data))

    def test_average_to_instantaneous_flattened(self):
        """
        Test proper de-averaging function. It shall work with multiple data. If a list with only one element is
        given, it shall return the same list reduced by one (list-) dimension.
        """
        data0 = DataExample(data=np.array([0, 0, 0, 0]))
        data1 = DataExample(data=np.array([2., 4, 6, 8]))
        data2 = DataExample(data=np.array([2., 4, 6, 8]) + 1)
        data3 = DataExample(data=np.array([2, 4, 6, 8]) + 0.5)

        # case 1 with only one element
        data_list1 = copy.deepcopy([[data1]])
        res_list1 = gt.average_to_instantaneous_flattened(data_list1, np.nan, 4)
        self.assertTrue(np.all(res_list1[0].data == data1.data))

        # case 2 with multiple elements
        data_list2 = copy.deepcopy([[data0, data1, data2, data3]])
        res_list2 = gt.average_to_instantaneous_flattened(data_list2, np.nan, 0)

        self.assertTrue(np.all(res_list2[0].data == data0.data))
        self.assertTrue(np.all(res_list2[1].data == data1.data))
        self.assertTrue(np.all(res_list2[2].data == 2 * data2.data - 1 * data1.data))
        self.assertTrue(np.all(res_list2[3].data == 3 * data3.data - 2 * data2.data))


class TestScalingTools(unittest.TestCase):
    def test_get_fill_value_unpack(self):
        """
        Only matrix elements with fill_value shall get the scaling factor.
        """
        fill_value = -1
        scale_factor_nc = 0.1
        gr_data = GeoReferencedData(data=np.array([[1, 2, fill_value], [4, fill_value, 6]]),
                                    data_description=DataDescription(fill_value=fill_value))
        idx_fill_value = gr_data.data == fill_value
        fill_value_unpack = st.get_fill_value_unpack(gr_data, scale_factor_nc)

        self.assertTrue(np.all(fill_value_unpack[idx_fill_value] == scale_factor_nc))
        self.assertTrue(np.all(fill_value_unpack[~idx_fill_value] != scale_factor_nc))

    def test_gr_data_scaling(self):
        """
        Test proper scaling of GeoReferencedData.
        """
        fill_value = -1
        scale_factor = 0.1
        scale_factor_nc = 0.01
        scale_undo = True
        gr_data = GeoReferencedData(data=np.array([[10, 25, fill_value], [42, fill_value, 67]]),
                                    data_description=DataDescription(fill_value=fill_value))
        idx_fill_value = gr_data.data == fill_value
        result = st.gr_data_scaling(gr_data, scale_undo, scale_factor, scale_factor_nc)

        self.assertTrue(np.all(result[idx_fill_value] == fill_value * scale_factor_nc))
        self.assertTrue(np.all(result[~idx_fill_value] == gr_data.data[~idx_fill_value] * scale_factor))


class DataExample:
    """
    Simple example class with data element.
    """
    def __init__(self,
                 data):
        """
        Instantiate DataExample class.

        :param data: some arbitrary data
        :type data: numpy.ndarray
        """
        self.data = data
