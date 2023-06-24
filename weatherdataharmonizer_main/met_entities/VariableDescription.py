#   @author: Michael Wagner
#   @organization: TU Dresden
#   @contact: michael.wagner@tu-dresden.de


class DataDescription:
    """
    DataDescription class contains metadata to further describe climate data to facilitate the export to netCDF.
    """
    def __init__(self,
                 fill_value=None,
                 units=None,
                 scale_factor=None,
                 coordinate_system=None,
                 long_name=None,
                 standard_name=None,
                 time_note=None,
                 eps_note=None):
        """
        Initialize DataDescription class.

        :param fill_value: value used to fill missing data
        :type fill_value: float, optional
        :param units: unit of data
        :type units: str, optional
        :param scale_factor: the value the data is multiplied with
        :type scale_factor: float, optional
        :param coordinate_system: description of the coordinate system (e.g. 'WGS 84, EPSG:4326')
        :type: str, optional
        :param long_name: long name of the data
        :type long_name: str, optional
        :param standard_name: a shorter standard name of the data
        :type standard_name: str, optional
        :param time_note: notes about temporal aspects, e.g. 'start at ...'
        :type time_note: str, optional
        :param eps_note: notes about ensemble prediction system members, e.g. 'eps member ...'
        :type eps_note: str, optional
        """
        self.fill_value = fill_value
        self.units = units
        self.scale_factor = scale_factor
        self.coordinate_system = coordinate_system
        self.long_name = long_name
        self.standard_name = standard_name
        self.time_note = time_note
        self.eps_note = eps_note


class TimeDescription:
    """
    TimeDescription class contains variables to further describe time aspects of climate data to facilitate the
    export to netCDF.
    """
    def __init__(self,
                 calendar=None,
                 units=None):
        """
        Initialize TimeDescription class.

        :param calendar: description of the calendar of the time data (e.g. 'standard')
        :type calendar: str, optional
        :param units: standard name of netCDF time unit (e.g. 'hours since 2018-01-01 01:00:00')
        :type units: str, optional
        """
        self.calendar = calendar
        self.units = units


class RegridDescription:
    """
    RegridDescription class contains variables do describe the regridding of weather data.
    """
    def __init__(self,
                 lon_target,
                 lat_target,
                 neighbors=3,
                 file_nearest=None):
        """
        Instantiate RegridDescription class.

        :param lon_target: longitudes of the target raster, given as raster with shape (num_lat, num_lon)
        :type lon_target: LonLatTime.LonLatTime
        :param lat_target: latitudes of the target raster, given as raster with shape (num_lat, num_lon)
        :type lat_target: LonLatTime.LonLatTime
        :param neighbors: number of neighbors for IDW
        :type neighbors: int, optional
        :param file_nearest: npz file with indexes and lengths for the current setup
        :type file_nearest: str, optional
        """
        self.lon_target = lon_target
        self.lat_target = lat_target
        self.neighbors = neighbors
        self.file_nearest = file_nearest


class CropDescription:
    """
    CropDescription class contains variables to describe the cropping of weather data.
    """
    def __init__(self,
                 lon_west=None,
                 lon_east=None,
                 lat_south=None,
                 lat_north=None,
                 idx_west=None,
                 idx_east=None,
                 idx_south=None,
                 idx_north=None,
                 idx_array=None):
        """
        Instantiate CropDescription class.

        :param lon_west: western longitude limit
        :type lon_west: float, optional
        :param lon_east: eastern longitude limit
        :type lon_east: float, optional
        :param lat_south: southern latitude limit
        :type lat_south: float, optional
        :param lat_north: northern latitude limit
        :type lat_north: float, optional
        :param idx_west: western index limit
        :type idx_west: int, optional
        :param idx_east: eastern index limit
        :type idx_east: int, optional
        :param idx_south: southern index limit
        :type idx_south: int, optional
        :param idx_north: northern index limit
        :type idx_north: int, optional
        :param idx_array: index for 1D array in lon and lat (e.g. original icond2 data)
        :type idx_array: np.ndarray, optional
        """
        self.lon_west = lon_west
        self.lon_east = lon_east
        self.lat_south = lat_south
        self.lat_north = lat_north
        self.idx_west = idx_west
        self.idx_east = idx_east
        self.idx_south = idx_south
        self.idx_north = idx_north
        self.idx_array = idx_array


class NcDimVarDescription:
    """
    NcVarDimDescription class describes detailed the dimension and variable names of a netcdf file.
    """
    def __init__(self,
                 dim_time=None,
                 dim_forecast=None,
                 dim_lon=None,
                 dim_lat=None,
                 dim_coord=None,
                 dim_eps=None,
                 var_time=None,
                 var_forecast=None,
                 var_lon=None,
                 var_lat=None,
                 var_eps=None,
                 var_data=None):
        """
        Initialize NcVarDimDescription class.

        :param dim_time: name of time dimension
        :type dim_time: str, optional
        :param dim_forecast: name of forecast dimension
        :type dim_forecast: str, optional
        :param dim_lon: name of longitude dimension
        :type dim_lon: str, optional
        :param dim_lat: name of latitude dimension
        :type dim_lat: str, optional
        :param dim_coord: name of coordinate dimension
        :type dim_coord: str, optional
        :param dim_eps: name of ensemble member dimension
        :type dim_eps: str, optional
        :param var_time: name of time variable
        :type var_time: str, optional
        :param var_forecast: name of forecast variable
        :type var_forecast: str, optional
        :param var_lon: name of longitude variable
        :type var_lon: str, optional
        :param var_lat: name of latitude variable
        :type var_lat: str, optional
        :param var_eps: name of ensemble member variable
        :type var_eps: str, optional
        :param var_data: name of data variable
        :type var_data: str, optional
        """
        self.dim_time = dim_time
        self.dim_forecast = dim_forecast
        self.dim_lon = dim_lon
        self.dim_lat = dim_lat
        self.dim_coord = dim_coord
        self.dim_eps = dim_eps
        self.var_time = var_time
        self.var_forecast = var_forecast
        self.var_lon = var_lon
        self.var_lat = var_lat
        self.var_eps = var_eps
        self.var_data = var_data
