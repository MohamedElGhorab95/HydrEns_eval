#   @author: Michael Wagner
#   @organization: TU Dresden
#   @contact: michael.wagner@tu-dresden.de


class Metadata:
    """
    Metadata class contains relevant metadata of cosmo formatted data.
    """
    def __init__(self,
                 centre=None,
                 centre_description=None,
                 datum=None,
                 datum_end_of_overall_time_interval=None,
                 datum_iso=None,
                 prediction_time=None,
                 step_type=None,
                 grid_type=None,
                 name=None,
                 units=None,
                 missing_value=None,
                 number_of_missing=None,
                 number_of_data_points=None,
                 fill_value=None):
        """
        Initialize Metadata class.

        :param centre: centre responsible for data creation
        :type centre: str
        :param centre_description: further description of the centre
        :type centre_description: str
        :param datum: UTC datum/time of simulation start
        :type datum: datetime.datetime
        :param datum_end_of_overall_time_interval: UTC datum/time of current simulation end
        :type datum_end_of_overall_time_interval: datetime.datetime
        :param datum_iso: ISO 8601 compliant datetime of current simulation end
        :type datum_iso: str
        :param prediction_time: timedelta from datum to datum_end_of_overall_time_interval
        :type prediction_time: datetime.timedelta
        :param step_type: type of step, e.g. 'accum' for accumulated data
        :type step_type: str
        :param grid_type: type of grid, e.g. 'unstructured grid' for non-orthogonal grids
        :type grid_type: str
        :param name: name of data
        :type name: str
        :param units: unit of data
        :type units: str
        :param missing_value: original missing value indicator in grib file
        :type missing_value: float
        :param number_of_missing: number of missing values
        :type number_of_missing: int
        :param number_of_data_points: total number of data points
        :type number_of_data_points: int
        :param fill_value: missing value indicator in corresponding read_cosmo.CosmoData.CosmoData object
        :type fill_value: float
        """
        self.centre = centre
        self.centre_description = centre_description
        self.datum = datum
        self.datum_end_of_overall_time_interval = datum_end_of_overall_time_interval
        self.datum_iso = datum_iso
        self.prediction_time = prediction_time
        self.step_type = step_type
        self.grid_type = grid_type
        self.name = name
        self.units = units
        self.missing_value = missing_value
        self.number_of_missing = number_of_missing
        self.number_of_data_points = number_of_data_points
        self.fill_value = fill_value
