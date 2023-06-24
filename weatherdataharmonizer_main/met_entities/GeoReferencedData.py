#   @author: Michael Wagner
#   @organization: TU Dresden
#   @contact: michael.wagner@tu-dresden.de

import os
import numpy as np
from tqdm import tqdm

import met_entities.LonLatTime as LonLatTime


class GeoReferencedData:
    """
    GeoReferencedData serves as a container for spatial data including coordinates, data and their description.
    """

    def __init__(self,
                 lon: LonLatTime = None,
                 lat: LonLatTime = None,
                 data=None,
                 data_description=None,
                 regridded=False,
                 cropped=False):
        """
        Initialize GeoReferencedData object.

        :param lon: longitudes of the spatial data
        :type lon: LonLatTime.LonLatTime, optional
        :param lat: latitudes of the spatial data
        :type lat: LonLatTime.LonLatTime, optional
        :param data: data with shape (num_lat, num_lon, num_forecast)
        :type data: ndarray, optional
        :param data_description: metadata for the data
        :type data_description: VariableDescription.DataDescription, optional
        :param regridded: indicator whether this instance of GeoReferencedData was already regridded
        :type regridded: bool, optional
        :param cropped: indicator whether this instance of GeoReferencedData was already cropped
        :type cropped: bool, optional
        """
        self.lon = lon
        self.lat = lat
        self.data = data
        self.data_description = data_description
        self.regridded = regridded
        self.cropped = cropped

    def regrid_idw(self, lon_target, lat_target, neighbors=3, file_nearest=None, short=None):
        """
        Interpolate the gridded data onto a new raster in the same coordinate system. It uses an inverse distance
        weighted (IDW) method with an arbitrary number of neighbors. The original class instance is changed. The
        function handles 1D/2D data without forecast as well as 2D/3D data including forecasts.

        :param lon_target: longitudes of the target raster, given as raster with shape (num_lat, num_lon)
        :type lon_target: LonLatTime.LonLatTime
        :param lat_target: latitudes of the target raster, given as raster with shape (num_lat, num_lon)
        :type lat_target: LonLatTime.LonLatTime
        :param neighbors: number of neighbors for IDW
        :type neighbors: int, optional
        :param file_nearest: npz file with indexes and lengths for the current setup; if the file does not exist, it is
            built, otherwise the content of the file is checked if it is suitable for the actual data
        :type file_nearest: str, optional
        :param short: if int16 or int32, the large data variables are cast to numpy.int16/int32 to minimize memory
            usage; the user must pay attention to the scale_factor to include the necessary precision
        :type short: str, optional
        """
        if self.cropped:
            print('the dataset is already cropped; consider cropping in one step with regridding\n')
        if self.regridded:
            print('the dataset is already regridded; consider regridding in one step\n')

        if self.lon.data.ndim == 2:
            # coordinates given as matrix
            coord_matrix = True
        else:
            # coordinates given as array
            coord_matrix = False

        # find the nearest points
        if file_nearest is None:
            # no filename is given
            idx_nearest, length_nearest = self.find_nearest(lon_target, lat_target, neighbors)
        elif file_nearest is not None and file_nearest[-4:] != '.npz':
            # the filename must have '.npz' extension
            raise Exception("the filename must end with '.npz'")
        elif file_nearest is not None and not os.path.exists(file_nearest):
            # filename is given but file does not exist -> save arrays to file
            idx_nearest, length_nearest = self.find_nearest(lon_target, lat_target, neighbors)
            # calculate meta parameters to be compared in case of array usage afterwards
            if coord_matrix:
                meta_param = np.array([self.lon.data[-1, 0] * 1e5, self.lat.data[-1, 0] * 1e5,
                                       self.lon.data[0, -1] * 1e5, self.lat.data[0, -1] * 1e5, self.lon.data.size,
                                       lon_target.data[-1, 0] * 1e5, lat_target.data[-1, 0] * 1e5,
                                       lon_target.data[0, -1] * 1e5, lat_target.data[0, -1] * 1e5,
                                       lon_target.data.size], dtype=int)
            else:
                meta_param = np.array([self.lon.data[0] * 1e5, self.lat.data[0] * 1e5, self.lon.data[-1] * 1e5,
                                       self.lat.data[-1] * 1e5, self.lon.data.size,
                                       lon_target.data[-1, 0] * 1e5, lat_target.data[-1, 0] * 1e5,
                                       lon_target.data[0, -1] * 1e5, lat_target.data[0, -1] * 1e5,
                                       lon_target.data.size], dtype=int)
            np.savez(file_nearest, idx_nearest=idx_nearest, length_nearest=length_nearest, meta_param=meta_param)
        else:
            # file with filename exists
            npz_file = np.load(file_nearest)
            meta_param = npz_file['meta_param']
            # compare meta parameters
            if coord_matrix:
                meta_param_act = np.array(
                    [self.lon.data[-1, 0] * 1e5, self.lat.data[-1, 0] * 1e5, self.lon.data[0, -1] * 1e5,
                     self.lat.data[0, -1] * 1e5, self.lon.data.size,
                     lon_target.data[-1, 0] * 1e5, lat_target.data[-1, 0] * 1e5, lon_target.data[0, -1] * 1e5,
                     lat_target.data[0, -1] * 1e5, lon_target.data.size], dtype=int)
            else:
                meta_param_act = np.array([self.lon.data[0] * 1e5, self.lat.data[0] * 1e5, self.lon.data[-1] * 1e5,
                                           self.lat.data[-1] * 1e5, self.lon.data.size,
                                           lon_target.data[-1, 0] * 1e5, lat_target.data[-1, 0] * 1e5,
                                           lon_target.data[0, -1] * 1e5, lat_target.data[0, -1] * 1e5,
                                           lon_target.data.size], dtype=int)
            if np.all(meta_param == meta_param_act):
                # the indexes and lengths can be used here
                idx_nearest = npz_file['idx_nearest']
                length_nearest = npz_file['length_nearest']
            else:
                # the file contains data for a different setup
                raise Exception(f"the file {file_nearest} was built for another setup and cannot be used here")

        # interpolate data at target points
        # IDW method: v_target = sum(weight*val_source) / sum(weight), if dist != 0
        #             v_target = val_source, if dist == 0
        #             weight = 1 / dist^p, with p = 2 (because of squarely growing point density in two dimensions)
        nrows_target, ncols_target = lon_target.data.shape
        if coord_matrix and len(self.data.shape) == 3:
            num_forecast = self.data.shape[2]
        elif not coord_matrix and len(self.data.shape) == 2:
            num_forecast = self.data.shape[1]
        elif coord_matrix and len(self.data.shape) == 2:
            num_forecast = 1
        elif not coord_matrix and len(self.data.shape) == 1:
            num_forecast = 1
        else:
            raise Exception('something went wrong with the dimensionality of input data')

        # cells in idx_nearest with -999 have no nearest neighbors -> missing values
        # temporarily declare with idx_nearest = 0 and fill with fill_value afterwards
        if coord_matrix:
            idx_nearest_missing = idx_nearest[:,:,0,0] == -999
        else:
            idx_nearest_missing = idx_nearest[:,:,0] == -999
        idx_nearest[idx_nearest == -999] = 0
        length_nearest[length_nearest == 0] = 1e5  # zero distances would result in ZeroDivisionError
        if num_forecast > 1:
            # 3D/2D data (including forecast)
            if coord_matrix:
                data_source_nearest = np.empty((nrows_target, ncols_target, neighbors, num_forecast))
                for fc in range(num_forecast):
                    for nb in range(neighbors):
                        data_source_nearest[:, :, nb, fc] = \
                            self.data[idx_nearest[:, :, 0, nb], idx_nearest[:, :, 1, nb], fc]
            else:
                data_source_nearest = self.data[idx_nearest]  # broadcasting applied
            # index for missing values in one of the neighbors -> in this case the interpolation result is not feasible
            idx_data_source_nearest_missing = np.any(data_source_nearest == self.data_description.fill_value, axis=2)
            weights = np.expand_dims(1 / (length_nearest ** 2), axis=3).repeat(num_forecast, axis=3)
            if short == 'int16':
                data_regridded = np.int16(np.around(np.sum(weights * data_source_nearest, axis=2) /
                                                    np.sum(weights, axis=2)))
            elif short == 'int32':
                data_regridded = np.int32(np.around(np.sum(weights * data_source_nearest, axis=2) /
                                                    np.sum(weights, axis=2)))
            else:
                data_regridded = np.sum(weights * data_source_nearest, axis=2) / np.sum(weights, axis=2)
        else:
            # 2D/1D data (no forecast)
            if coord_matrix:
                data_source_nearest = np.empty((nrows_target, ncols_target, neighbors))
                for nb in range(neighbors):
                    data_source_nearest[:, :, nb] = self.data[idx_nearest[:, :, 0, nb], idx_nearest[:, :, 1, nb]]
            else:
                data_source_nearest = self.data[idx_nearest]  # broadcasting applied
            # index for missing values in one of the neighbors -> in this case the interpolation result is not feasible
            idx_data_source_nearest_missing = np.any(data_source_nearest == self.data_description.fill_value, axis=2)
            weights = 1 / (length_nearest ** 2)
            if short == 'int16':
                data_regridded = np.int16(np.around(np.sum((weights * data_source_nearest), axis=2) /
                                                    np.sum(weights, axis=2)))
            elif short == 'int32':
                data_regridded = np.int32(np.around(np.sum((weights * data_source_nearest), axis=2) /
                                                    np.sum(weights, axis=2)))
            else:
                data_regridded = np.sum((weights * data_source_nearest), axis=2) / np.sum(weights, axis=2)

        data_regridded[idx_nearest_missing] = self.data_description.fill_value
        data_regridded[idx_data_source_nearest_missing] = self.data_description.fill_value

        self.data = data_regridded
        self.lon.data = lon_target.data
        self.lat.data = lat_target.data
        self.regridded = True

    def find_nearest(self, lon_target, lat_target, neighbors):
        """
        Find an arbitrary number of the nearest points from the source coordinates for the target grid. It works with
        source coordinates given as 2D matrix and given as 1D array.

        :param lon_target: longitudes of the target raster, given as raster with shape (num_lat, num_lon)
        :type lon_target: LonLatTime.LonLatTime
        :param lat_target: latitudes of the target raster, given as raster with shape (num_lat, num_lon)
        :type lat_target: LonLatTime.LonLatTime
        :param neighbors: number of neighbors for IDW
        :type neighbors: int
        :return: (1) matrix with indexes of the nearest points for source coordinates given as matrix with shape
            (num_lat, num_lon, 2, neighbors), whereas the third dimension differentiates between index of row and index
            of col or for source coordinates given as array with shape (num_lat, num_lon, neighbors); (2) matrix with
            distances of the nearest points with shape (num_lat, num_lon, neighbors)
        :rtype: numpy.ndarray, numpy.ndarray
        """
        nrows_target, ncols_target = lon_target.data.shape
        length_nearest = np.empty((nrows_target, ncols_target, neighbors))
        if self.lon.data.ndim == 2:
            # coordinates given as matrix
            coord_matrix = True
            idx_nearest = np.empty((nrows_target, ncols_target, 2, neighbors), dtype=int)
            lon_max = (self.lon.data[:, 1:] - self.lon.data[:, 0:-1]).max()
            lat_max = (self.lat.data[0:-1, :] - self.lat.data[1:, :]).max()
        else:
            # coordinates given as array
            coord_matrix = False
            idx_nearest = np.empty((nrows_target, ncols_target, neighbors), dtype=int)
            lon_sorted = np.sort(self.lon.data)
            lat_sorted = np.sort(self.lat.data)
            lon_max = (lon_sorted[1:] - lon_sorted[0:-1]).max()
            lat_max = (lat_sorted[1:] - lat_sorted[0:-1]).max()

        for col in tqdm(range(ncols_target), bar_format='{l_bar}{bar}| {n_fmt}/{total_fmt} [{elapsed}<{remaining}]',
                        delay=2):
            for row in range(nrows_target):
                if coord_matrix:
                    idx_source_act1 = np.logical_and(self.lon.data.ravel() > lon_target.data[row, col] - 4 * lon_max,
                                                     self.lon.data.ravel() < lon_target.data[row, col] + 4 * lon_max)
                    idx_source_act2 = np.logical_and(self.lat.data.ravel() > lat_target.data[row, col] - 4 * lat_max,
                                                     self.lat.data.ravel() < lat_target.data[row, col] + 4 * lat_max)
                else:
                    idx_source_act1 = np.logical_and(self.lon.data > lon_target.data[row, col] - 4 * lon_max,
                                                     self.lon.data < lon_target.data[row, col] + 4 * lon_max)
                    idx_source_act2 = np.logical_and(self.lat.data > lat_target.data[row, col] - 4 * lat_max,
                                                     self.lat.data < lat_target.data[row, col] + 4 * lat_max)
                idx_source_act = np.where(np.logical_and(idx_source_act1, idx_source_act2))[0]
                if idx_source_act.size < neighbors:
                    # no or not enough neighbors found, target possibly out of domain -> mark with -999
                    print(f'no or not enough neighbors found for lon={lon_target.data[row, col]} and lat={lat_target.data[row, col]}; target out of domain')
                    if coord_matrix:
                        idx_nearest[row, col, :, :] = -999
                        length_nearest[row, col, :] = -999
                    else:
                        idx_nearest[row, col, :] = -999
                        length_nearest[row, col, :] = -999
                    continue
                # elif idx_source_act.size < neighbors:
                #     raise Exception(f"found less than {neighbors} nearest neighbors; consider a lower number")

                if coord_matrix:
                    lon_source_act = self.lon.data.ravel()[idx_source_act]
                    lat_source_act = self.lat.data.ravel()[idx_source_act]
                else:
                    lon_source_act = self.lon.data[idx_source_act]
                    lat_source_act = self.lat.data[idx_source_act]
                distance = np.sqrt((lon_source_act - lon_target.data[row, col]) ** 2 +
                                   (lat_source_act - lat_target.data[row, col]) ** 2)
                idx_sort = np.argsort(distance)
                length_nearest[row, col, :] = distance[idx_sort[0:neighbors]]
                if coord_matrix:
                    tmp = np.unravel_index(idx_source_act[idx_sort[0:neighbors]], self.lon.data.shape)
                    idx_nearest[row, col, 0, :] = tmp[0]
                    idx_nearest[row, col, 1, :] = tmp[1]
                else:
                    idx_nearest[row, col, :] = idx_source_act[idx_sort[0:neighbors]]

        return idx_nearest, length_nearest

    def find_nearest_slower(self, lon_target, lat_target, neighbors):
        """
        Find an arbitrary number of the nearest points from the source grid for the target grid. The algorithm is
        somewhat easier to understand, but significantly slower than find_nearest. The code is not extended to work with
        source coordinates given as an array (not a matrix).

        :param lon_target: longitudes of the target raster, given as raster with shape (num_lat, num_lon)
        :type lon_target: LonLatTime.LonLatTime
        :param lat_target: latitudes of the target raster, given as raster with shape (num_lat, num_lon)
        :type lat_target: LonLatTime.LonLatTime
        :param neighbors: number of neighbors for IDW
        :type neighbors: int
        :return: (1) matrix with indexes of nearest points with shape (num_lat, num_lon, 2, neighbors), whereas the
            third dimension differentiates between index of row and index of col; (2) matrix with distances of the
            nearest points with shape (num_lat, num_lon, neighbors)
        :rtype: numpy.ndarray, numpy.ndarray
        """
        nrows_target, ncols_target = lon_target.data.shape
        idx_nearest = np.empty((nrows_target, ncols_target, 2, neighbors), dtype=int)
        length_nearest = np.empty((nrows_target, ncols_target, neighbors))

        # find the nearest points
        for col in tqdm(range(ncols_target), bar_format='{l_bar}{bar}| {n_fmt}/{total_fmt} [{elapsed}<{remaining}]'):
            for row in range(nrows_target):
                distance = np.sqrt((self.lon.data - lon_target.data[row, col]) ** 2 +
                                   (self.lat.data - lat_target.data[row, col]) ** 2).ravel()
                idx_sort = np.argsort(distance)
                length_nearest[row, col, :] = distance[idx_sort[0:neighbors]]
                tmp = np.unravel_index(idx_sort[0:neighbors], self.lat.data.shape)
                idx_nearest[row, col, 0, :] = tmp[0]
                idx_nearest[row, col, 1, :] = tmp[1]

        return idx_nearest, length_nearest

    def crop(self, lon_west=None, lon_east=None, lat_south=None, lat_north=None, idx_west=None, idx_east=None,
             idx_south=None, idx_north=None, idx_array=None):
        """
        Cropping of data of this GeoReferencedData class. Usage of indexes directly if given. Otherwise, use lon/lat
        with the guarantee that the whole requested area is within the cropped region. The function handles 1D/2D data
        without forecast as well as 2D/3D data including forecasts.

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
        if self.regridded:
            print('the dataset is already regridded; consider cropping in one step with regridding\n')
        if self.cropped:
            print('the dataset is already cropped; consider cropping in one step\n')

        if self.lon.data.ndim == 2:
            # coordinates given as matrix
            coord_matrix = True
        else:
            # coordinates given as array
            coord_matrix = False

        if idx_west is not None and idx_east is not None and idx_south is not None and idx_north is not None:
            # use indexes directly
            if coord_matrix:
                nrows = np.size(self.data, axis=0)
                ncols = np.size(self.data, axis=1)
                if idx_west < 0 or idx_west >= ncols or \
                        idx_east < 0 or idx_east >= ncols or \
                        idx_south < 0 or idx_south >= nrows or \
                        idx_north < 0 or idx_north >= nrows:
                    raise ValueError('the cropping indexes are out of the relevant domain')
                if idx_east < idx_west:
                    raise ValueError('idx_east must surpass idx_west')
                if idx_south < idx_north:
                    raise ValueError('idx_south must surpass idx_north')

        elif lon_west is not None and lon_east is not None and lat_south is not None and lat_north is not None:
            # calculate the indexes using coordinates
            if lon_east < lon_west:
                raise ValueError('lon_east must surpass lon_west')
            if lat_north < lat_south:
                raise ValueError('lat_north must surpass lat_south')
            if lon_west > np.max(self.lon.data) or lon_east < np.min(self.lon.data) or \
                    lat_south > np.max(self.lat.data) or lat_north < np.min(self.lat.data):
                raise ValueError('the cropping coordinates are out of the relevant domain')

            if coord_matrix:
                # get maxima and minima per col/row and find the corresponding indexes
                max_lon = np.max(self.lon.data, axis=0)
                min_lon = np.min(self.lon.data, axis=0)
                min_lat = np.min(self.lat.data, axis=1)
                max_lat = np.max(self.lat.data, axis=1)

                if np.min(max_lon) >= lon_west:
                    idx_west = 0
                else:
                    # with -1 the col before the first one is chosen (floor of cols)
                    idx_west = np.where(~(max_lon < lon_west))[0][0] - 1
                if np.max(min_lon) < lon_east:
                    idx_east = len(min_lon) - 1
                else:
                    idx_east = np.where(min_lon > lon_east)[0][0]
                if np.min(max_lat) > lat_south:
                    idx_south = len(min_lat) - 1
                else:
                    idx_south = np.where(max_lat < lat_south)[0][0]
                if np.max(min_lat) < lat_north:
                    idx_north = 0
                else:
                    # with -1 the row before the first one is chosen (floor of rows)
                    idx_north = np.where(~(min_lat > lat_north))[0][0] - 1
            else:
                # get indexes for coordinate arrays
                idx_array1 = np.logical_and(self.lon.data >= lon_west, self.lon.data <= lon_east)
                idx_array2 = np.logical_and(self.lat.data >= lat_south, self.lat.data <= lat_north)
                idx_array = np.logical_and(idx_array1, idx_array2)
        else:
            raise Exception('please provide either longitudes and latitudes or indexes in every direction')

        # crop relevant variables in self
        if self.data.ndim == 3:
            # 3D data (included forecast)
            self.data = self.data[idx_north:idx_south + 1, idx_west:idx_east + 1, :]
        elif self.data.ndim == 2 and coord_matrix:
            # 2D data (no forecast)
            self.data = self.data[idx_north:idx_south + 1, idx_west:idx_east + 1]
        elif self.data.ndim == 2 and not coord_matrix:
            # 2D data (included forecast)
            self.data = self.data[idx_array, :]
        else:
            # 1D data (no forecast)
            self.data = self.data[idx_array]

        if coord_matrix:
            self.lon.data = self.lon.data[idx_north:idx_south + 1, idx_west:idx_east + 1]
            self.lat.data = self.lat.data[idx_north:idx_south + 1, idx_west:idx_east + 1]
        else:
            self.lon.data = self.lon.data[idx_array]
            self.lat.data = self.lat.data[idx_array]

        self.cropped = True
