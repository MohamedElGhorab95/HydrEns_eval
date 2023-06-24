#   @author: Michael Wagner
#   @organization: TU Dresden
#   @contact: michael.wagner@tu-dresden.de
import bz2
import datetime as dt
import gzip
import re

import numpy as np
from pyproj import CRS
from pyproj import Transformer

from met_entities.LonLatTime import LonLatTime
from met_entities.VariableDescription import DataDescription
from read_radolan.Metadata import Metadata
from read_radolan.RadolanData import RadolanData


def read_radolan(filename, scale_factor=1, fill_value=np.nan, metadata_only=False, short=None):
    """
    Read in radolan formatted data as bz2 or gz compressed data or as binary.

    :param filename: name of file; distinction is made upon the file extension (.bz2, .gz, or nothing)
    :type filename: str
    :param scale_factor: the final data has to be multiplied with this value
    :type scale_factor: float, optional
    :param fill_value: missing data is filled with that value
    :type fill_value: float, optional
    :param metadata_only: if True, only the metadata of the file is extracted
    :type metadata_only: bool, optional
    :param short: if int16 or int32, the large data variables are cast to numpy.int16/int32 to minimize memory
            usage; the user must pay attention to the scale_factor to include the necessary precision
    :type short: str, optional
    :return: RadolanData class or Metadata class with corresponding content
    :rtype: read_radolan.RadolanData.RadolanData, read_radolan.RadolanData.Metadata
    """
    if filename[-3:] == '.gz':
        # gzip compressed data
        with gzip.open(filename, 'rb') as fid:
            file_content = fid.read()
    elif filename[-4:] == '.bz2':
        # bzip2 compressed data
        with bz2.open(filename, 'rb') as fid:
            file_content = fid.read()
    elif filename[-4:].find('.') == -1:
        # uncompressed data
        with open(filename, 'rb') as fid:
            file_content = fid.read()
    else:
        raise Exception(f'cannot interpret format of {filename}')
    return read_radolan_data(file_content, scale_factor, fill_value, metadata_only, short)


def read_radolan_gz(filename, scale_factor=1, fill_value=np.nan, metadata_only=False, short=None):
    """
    Read in radolan formatted and gzip compressed data.

    :param filename: name of file
    :type filename: str
    :param scale_factor: the final data has to be multiplied with this value
    :type scale_factor: float, optional
    :param fill_value: missing data is filled with that value
    :type fill_value: float, optional
    :param metadata_only: if True, only the metadata of the file is extracted
    :type metadata_only: bool, optional
    :param short: if int16 or int32, the large data variables are cast to numpy.int16/int32 to minimize memory
        usage; the user must pay attention to the scale_factor to include the necessary precision
    :type short: str, optional
    :return: RadolanData class or Metadata class with corresponding content
    :rtype: read_radolan.RadolanData.RadolanData, read_radolan.RadolanData.Metadata
    """
    with gzip.open(filename, 'rb') as fid:
        file_content = fid.read()
    return read_radolan_data(file_content, scale_factor, fill_value, metadata_only, short)


def read_radolan_bz2(filename, scale_factor=1, fill_value=np.nan, metadata_only=False, short=None):
    """
    Read in radolan formatted and bzip2 compressed data.

    :param filename: name of file
    :type filename: str
    :param scale_factor: the final data has to be multiplied with this value
    :type scale_factor: float, optional
    :param fill_value: missing data is filled with that value
    :type fill_value: float, optional
    :param metadata_only: if True, only the metadata of the file is extracted
    :type metadata_only: bool, optional
    :param short: if int16 or int32, the large data variables are cast to numpy.int16/int32 to minimize memory
        usage; the user must pay attention to the scale_factor to include the necessary precision
    :type short: str, optional
    :return: RadolanData class or Metadata class with corresponding content
    :rtype: read_radolan.RadolanData.RadolanData, read_radolan.RadolanData.Metadata
    """
    with bz2.open(filename, 'rb') as fid:
        file_content = fid.read()
    return read_radolan_data(file_content, scale_factor, fill_value, metadata_only, short)


def read_radolan_binary(filename, scale_factor=1, fill_value=np.nan, metadata_only=False, short=None):
    """
    Read in radolan formatted binary data.

    :param filename: name of file
    :type filename: str
    :param scale_factor: the final data has to be multiplied with this value
    :type scale_factor: float, optional
    :param fill_value: missing data is filled with that value
    :type fill_value: float, optional
    :param metadata_only: if True, only the metadata of the file is extracted
    :type metadata_only: bool, optional
    :param short: if int16 or int32, the large data variables are cast to numpy.int16/int32 to minimize memory
        usage; the user must pay attention to the scale_factor to include the necessary precision
    :type short: str, optional
    :return: RadolanData class or Metadata class with corresponding content
    :rtype: read_radolan.RadolanData.RadolanData, read_radolan.RadolanData.Metadata
    """
    with open(filename, 'rb') as fid:
        file_content = fid.read()
    return read_radolan_data(file_content, scale_factor, fill_value, metadata_only, short)


def read_radolan_data(stream, scale_factor=1, fill_value=np.nan, metadata_only=False, short=None):
    """
    Read radolan formatted data from stream.

    :param stream: the data stream
    :type stream: bytes
    :param scale_factor: the final data has to be multiplied with this value
    :type scale_factor: float, optional
    :param fill_value: missing data is filled with that value
    :type fill_value: float, optional
    :param metadata_only: if True, only the metadata of the file is extracted
    :type metadata_only: bool, optional
    :param short: if int16 or int32, the large data variables are cast to numpy.int16/int32 to minimize memory
        usage; the user must pay attention to the scale_factor to include the necessary precision
    :type short: str, optional
    :return: RadolanData class or Metadata class with corresponding content
    :rtype: read_radolan.RadolanData.RadolanData, read_radolan.RadolanData.Metadata
    """
    if short and scale_factor >= 1:
        print(f'the attribute short casts the data to 16 bit signed integer for memory saving; the scale_factor = '
              f'{scale_factor} imposes a low precision')
    if short and np.isnan(fill_value):
        raise Exception(f'if short data type is requested, the fill_value must not be nan')
    if short and (short != 'int16' and short != 'int32'):
        raise Exception(f'datatype {short} not supported; please change to int16 or int32')

    # read header data and fill metadata
    md = Metadata()
    header = (stream[0:300]).decode(errors='ignore')
    md.product = header[0:2]
    supported_products = ['RW', 'RX', 'YW', 'RQ', 'RE', 'RV', 'RY']
    if md.product not in supported_products:
        raise Exception(f"product {md.product} not supported")

    md.datum.append(2000 + int(header[15:17]))  # year
    md.datum.append(int(header[13:15]))  # month
    md.datum.append(int(header[2:4]))  # day
    md.datum.append(int(header[4:6]))  # hour
    md.datum.append(int(header[6:8]))  # minute
    tmp = dt.datetime(md.datum[0], md.datum[1], md.datum[2], md.datum[3], md.datum[4])
    md.datum_iso = tmp.replace(tzinfo=dt.timezone.utc).isoformat()

    idx_by = re.search('BY', header).end()
    idx_vs = re.search('VS', header).end()
    md.product_length = header[idx_by:idx_vs - 2]  # product length in byte
    idx_sw = re.search('SW', header).end()
    md.format_version = int(header[idx_vs:idx_sw - 2])
    idx_pr = re.search('PR', header).end()
    md.version = header[idx_sw:idx_pr - 2]  # software version
    idx_int = re.search('INT', header).end()
    scale_factor_internal = 10 ** int(header[idx_pr + 2:idx_int - 3])  # data precision
    idx_gp = re.search('GP', header).end()
    header_int = header[idx_int:idx_gp - 2]
    if 'U' in header_int:
        # some products have a 'U' for the unit in this element
        header_int = header_int.split('U')
        if int(header_int[1]) == 0:
            # interval length in minutes
            md.interval_length = int(header_int[0])
        else:
            # interval length in days
            md.interval_length = int(header_int[0]) * 60 * 24
    else:
        md.interval_length = int(header[idx_int:idx_gp - 2])  # length of time interval

    nrows = int(header[idx_gp:idx_gp + 4])
    ncols = int(header[idx_gp + 5:idx_gp + 9])
    md.nrows = nrows
    md.ncols = ncols

    if 'VV' in header:
        # prediction time only given for forecast products
        idx_vv = re.search('VV', header).end()
        md.prediction_time = int(header[idx_vv:idx_vv + 4])  # prediction time in min

    if 'MF' in header:
        idx_mf = re.search('MF', header).end()
        md.module_flags = int(header[idx_mf:idx_mf + 9])  # module flags

    if 'QN' in header:
        # kind of quantification only given for RQ and RE products
        idx_qn = re.search('QN', header).end()
        md.quantification_kind = int(header[idx_qn:idx_qn + 4])

    idx_ms = re.search('MS', header).end()
    length_text = int(header[idx_ms:idx_ms + 3])
    md.sites = header[idx_ms + 3:idx_ms + 3 + length_text]

    if metadata_only:
        # if metadata only is enough
        return md
    else:
        # if also the data is demanded
        rd = RadolanData()
        rd.metadata = md  # fill metadata variable in radolan data

        header_length = idx_ms + 3 + length_text  # one empty char before byte block starts

        # read  data block
        if md.product in ['RW', 'YW', 'RQ', 'RE', 'RV', 'RY']:
            # the original 16-bit data is departed into to 8-bit parts
            data_raw = np.frombuffer(stream[header_length + 1:], dtype=np.uint8)
            data_bit = np.unpackbits(data_raw, bitorder='little')
            num_bit = len(data_bit)
            idx_clutter = data_bit[np.arange(15, num_bit, 16)]
            idx_nan = data_bit[np.arange(13, num_bit, 16)]
            idx_missing = np.array(np.logical_or(idx_clutter, idx_nan))

            data_bit[np.arange(12, num_bit, 16)] = 0
            data_bit[np.arange(13, num_bit, 16)] = 0
            data_bit[np.arange(14, num_bit, 16)] = 0
            data_bit[np.arange(15, num_bit, 16)] = 0
            tmp = np.reshape(data_bit, (int(num_bit/16), 16))
            p_tmp = np.dot(tmp, (2**np.arange(0, 16))) * scale_factor_internal / scale_factor
            p_tmp[np.array(idx_missing, dtype=bool)] = fill_value
            data_tmp = p_tmp

        elif md.product == 'RX':
            data_raw = np.frombuffer(stream[header_length + 1:], dtype=np.uint8)
            rvp6 = data_raw  # TODO: is a conversion to float necessary (like in matlab)?
            idx_nan = rvp6 == 250
            rvp6[idx_nan] = fill_value
            idx_clutter = rvp6 == 249
            rvp6[idx_clutter] = fill_value
            dbz = rvp6 / 2 - 32.5
            data_tmp = dbz

        else:
            raise Exception(f"Product {md.product} is not supported yet")

        # convert data block to nrows x ncols in the right orientation
        if short == 'int16':
            rd.data = np.short(np.around(np.flipud(np.reshape(data_tmp, (nrows, ncols)))))
        elif short == 'int32':
            rd.data = np.int32(np.around(np.flipud(np.reshape(data_tmp, (nrows, ncols)))))
        else:
            rd.data = np.flipud(np.reshape(data_tmp, (nrows, ncols)))
        rd.idx_clutter = np.flipud(np.reshape(idx_clutter, (nrows, ncols)))
        rd.idx_nan = np.flipud(np.reshape(idx_nan, (nrows, ncols)))

        return rd


def get_lonlat(format_version, grid_variant=None, num_rows=None, target_epsg=4326, corners=False):
    """
    Calculate longitudes and latitudes for centerpoints of a given Radolan raster. The format version (identifier VS
    in radolan binary files) must be given, as until 4 a sphere earth model is used and from 5 on the WGS84 ellipsoid
    as earth model is applied.

    :param format_version: version from VS identifier in radolan binary files; relevant is only <=4 or >=5
    :type format_version: int
    :param grid_variant: either radolanrx, de1200, de4800, or europecomposite (case-insensitive); either grid variant or
        num_rows must be given
    :type grid_variant: str, optional
    :param num_rows: either 900, 1200, 4800, or 2400; either grid variant or num_rows must be given
    :type num_rows: int, optional
    :param target_epsg: epsg number of target crs (default: 4326)
    :type target_epsg: int, optional
    :param corners: if true, the outer coordinates of all four corners are returned (default: False)
    :type corners: bool, optional
    :return: longitude and latitude objects
    :rtype: met_entities.LonLatTime.LonLatTime, met_entities.LonLatTime.LonLatTime
    """
    grids_supported = {'radolanrx': 900,
                       'de1200': 1200,
                       'de4800': 4800,
                       'europecomposite': 2400}
    rows_supported = {val: key for (key, val) in grids_supported.items()}  # reversed dictionary
    if grid_variant is None and num_rows is not None:
        grid_variant = rows_supported[num_rows]

    # shortcut to older but faster coordinate transformation
    if format_version <= 4 and target_epsg == 4326:
        if grid_variant == 'radolanrx':
            return get_lonlat_sphere(900, corners=corners)
        elif num_rows == 1100:
            return get_lonlat_sphere(1100, corners=corners)

    # regular coordinate transformation
    if num_rows is None and grid_variant is None:
        raise Exception('Please provide either number of rows or the grid variant')
    elif num_rows == 1100:
        raise Exception('The number of rows is given as 1100. This is the extended national composite. This is only '
                        'supported via get_lonlat_sphere function')
    elif num_rows is not None and num_rows not in grids_supported.values():
        raise Exception(f'A number of {num_rows} rows is not supported')
    elif grid_variant is not None and grid_variant.lower() not in grids_supported.keys():
        raise Exception(f'The grid variant {grid_variant.lower()} is not supported. Please choose one of '
                        '{grids_supported.keys()}')

    # lonlat in WGS84 as target crs is provided by default
    crs_target = CRS.from_epsg(target_epsg)

    if grid_variant.lower() == 'radolanrx':
        if format_version <= 4:
            # RadolanRX sphere earth model
            crs_source = CRS.from_proj4('+proj=stere +lat_0=90 +lat_ts=60 +lon_0=10 +a=6370040 +b=6370040 +no_defs '
                                        '+x_0=522962.16692185635 +y_0=3759144.724265574')
        else:
            # RadolanRX WGS84 earth model
            crs_source = CRS.from_proj4('+proj=stere +lat_0=90 +lat_ts=60 +lon_0=10 +a=6378137 +b=6356752.3142451802 '
                                        '+no_defs +x_0=523196.83521777776 +y_0=3772588.861931134')
        if corners:
            x = [-.5 * 1000, 899.5 * 1000]
            y = [-899.5 * 1000, .5 * 1000]
        else:
            x = np.arange(0, 900) * 1000  # coordinates necessary in meters
            y = np.arange(-899, 1) * 1000
    elif grid_variant.lower() == 'de1200':
        if format_version <= 4:
            # DE1200 grid sphere earth model
            crs_source = CRS.from_proj4('+proj=stere +lat_0=90 +lat_ts=60 +lon_0=10 +a=6370040 +b=6370040 +no_defs '
                                        '+x_0=542962.16692185658 +y_0=3609144.7242655745')
        else:
            # DE1200 grid WGS84 earth model
            crs_source = CRS.from_proj4('+proj=stere +lat_0=90 +lat_ts=60 +lon_0=10 +a=6378137 +b=6356752.3142451802 '
                                        '+no_defs +x_0=543196.83521776402 +y_0=3622588.8619310018')
        if corners:
            x = [-.5 * 1000, 1099.5 * 1000]
            y = [-1199.5 * 1000, .5 * 1000]
        else:
            x = np.arange(0, 1100) * 1000
            y = np.arange(-1199, 1) * 1000
    elif grid_variant.lower() == 'de4800':
        if format_version <= 4:
            # DE4800 grid sphere earth model
            crs_source = CRS.from_proj4('+proj=stere +lat_0=90 +lat_ts=60 +lon_0=10 +a=6370040 +b=6370040 +no_defs '
                                        '+x_0=543337.16692185646 +y_0=3608769.7242655735')
        else:
            # DE4800 grid WGS84 earth model
            crs_source = CRS.from_proj4('+proj=stere +lat_0=90 +lat_ts=60 +lon_0=10 +a=6378137 +b=6356752.3142451802 '
                                        '+no_defs +x_0=543571.83521776402 +y_0=3622213.8619310018')
        if corners:
            x = [-.125 * 1000, 1099.875 * 1000]
            y = [-1199.875 * 1000, .125 * 1000]
        else:
            x = np.arange(0, 1100, .25) * 1000
            y = np.arange(-1199.75, .25, .25) * 1000
    elif grid_variant.lower() == 'europecomposite':
        if format_version <= 4:
            # Europe composite sphere earth model
            crs_source = CRS.from_proj4('+proj=stere +lat_0=90 +lon_0=10 +lat_ts=60 +a=6370040 +b=6370040 +no_defs '
                                        '+x_0=1622962.1669218568 +y_0=3059144.724265574')
        else:
            # Europe composite WGS84 earth model
            crs_source = CRS.from_proj4('+proj=stere +lat_0=90 +lat_ts=60 +lon_0=10 +a=6378137 +b=6356752.3142451802 '
                                        '+no_defs +x_0=1623196.8352178331 +y_0=3072588.8619308411')
        if corners:
            x = [-.5 * 1000, 2399.5 * 1000]
            y = [-2399.5 * 1000, .5 * 1000]
        else:
            x = np.arange(0, 2400) * 1000
            y = np.arange(-2399, 1) * 1000
    else:
        crs_source = None
        x = None
        y = None

    # transform coordinates to target crs
    transformer = Transformer.from_crs(crs_from=crs_source, crs_to=crs_target, always_xy=True)

    x_mat, y_mat = np.meshgrid(x, y)
    lonlat = transformer.transform(x_mat, y_mat)
    lon_description = DataDescription(units='degrees_east',
                                      long_name='longitude of center',
                                      coordinate_system='WGS84, EPSG:4326')
    lon = LonLatTime(data=np.flipud(lonlat[0]), data_description=lon_description)

    lat_description = DataDescription(units='degrees_north',
                                      long_name='latitude of center',
                                      coordinate_system='WGS84, EPSG:4326')
    lat = LonLatTime(data=np.flipud(lonlat[1]), data_description=lat_description)

    return lon, lat


def get_lonlat_sphere(num_rows, corners=False):
    """
    Calculate longitudes and latitudes for centerpoints of the currently used Radolan rasters (either 900 or 1100 rows
    long). The eastbound shift of 1100 rows products is respected.

    :param num_rows: number of rows (900 or 1100)
    :type num_rows: int
    :param corners: if true, the outer coordinates of all four corners are returned (default: False)
    :type corners: bool, optional
    :return: longitude and latitude objects
    :rtype: met_entities.LonLatTime.LonLatTime, met_entities.LonLatTime.LonLatTime
    """
    radius = 6370.04
    phi_0 = 60 * np.pi / 180
    lambda_0 = 10 * np.pi / 180

    if num_rows == 900:
        if corners:
            x = [-523.4622, -523.4622 + 900]
            y = [-4658.645, -4658.645 + 900]
        else:
            x = -523.4622 + .5 + np.arange(0, 900)
            y = -4658.645 + .5 + np.arange(0, 900)
    elif num_rows == 1100:
        if corners:
            x = [-523.4622 + 80, -523.4622 + 80 + 900]
            y = [-4658.645 - 100, -4658.645 - 100 + 1100]
        else:
            x = -523.4622 + .5 + 80 + np.arange(0, 900)
            y = -4658.645 + .5 - 100 + np.arange(0, 1100)
    else:
        raise ValueError(f"Coordinates for {num_rows} rows cannot be calculated, please change to either 900 or 1100")

    x_mat, y_mat = np.meshgrid(x, y)
    lambda_r = np.flipud((np.arctan(-x_mat / y_mat) + lambda_0) * 180 / np.pi)
    phi_r = np.flipud((np.arcsin((radius ** 2 * (1 + np.sin(phi_0)) ** 2 - (x_mat ** 2 + y_mat ** 2)) /
                                 (radius ** 2 * (1 + np.sin(phi_0)) ** 2 + (x_mat ** 2 + y_mat ** 2)))) * 180 / np.pi)

    lon_description = DataDescription(units='degrees_east',
                                      long_name='longitude of center',
                                      coordinate_system='WGS84, EPSG:4326')
    lon = LonLatTime(data=lambda_r, data_description=lon_description)

    lat_description = DataDescription(units='degrees_north',
                                      long_name='latitude of center',
                                      coordinate_system='WGS84, EPSG:4326')
    lat = LonLatTime(data=phi_r, data_description=lat_description)

    return lon, lat
