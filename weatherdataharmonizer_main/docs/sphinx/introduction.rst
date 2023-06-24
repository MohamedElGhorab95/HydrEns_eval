Introduction
============

The weatherDataHarmonizer project facilitates the harmonizing of different weather data of the
German Weatherservice (DWD). It incorporates recent data as well as forecast
data from RadolanRW, RadvorRQ, RadolanRV, Icon-D2, Icon-D2-EPS, Icon-EU, and Icon-EU-EPS products.

This chapter includes examples for typical usages. The comprehensive technical description follows in the next chapter.

Content and usage for fast readers
----------------------------------

The weatherDataHarmonizer ist parted into different packages (e.g. ``met_entities``,
``read_radolan``, ``read_icon``), which can be imported into own projects. Various
example implementations are given in https://gitlab.hrz.tu-chemnitz.de/ihmtud/weatherdataharmonizer/-/blob/main/eval_met_packages.py.

The code examples below show low level access using modules in ``read_radolan``, and high level access for advanced
reading and harmonizing features.


Example: low level access of a raw dataset
------------------------------------------

Using the low level access variant, it is possible to directly load a raw file of the appropriate data type. Currently
supported are data of type Radolan, Icon, and Cosmo.

The following piece of code reads in a Radolan type data file and stores the data in a ``RadolanData`` class instance. Further,
the longitudes and latitudes of all centerpoints of the original 900 x 900 grid are given in a tuple of ``LonLatTime`` class instances.
The ``readRadolan.read_radolan`` method allows further attributes, e.g. ``scale_factor`` or ``fill_value`` to further
influence the file reading.

.. code-block:: python

    from read_radolan import readRadolan as readRadolan

    data1 = readRadolan.read_radolan('path/to/radolan/data_file')
    lonlat = readRadolan.get_lonlat(format_version=4, grid_variant='radolanrx')


Example: read/export one raw dataset
------------------------------------

This code example defines a target grid, reads in one time step of RadolanRW data,
regrids the data using IDW method with the nearest three neighbors, and exports the
resulting grid to a netcdf file.

.. code-block:: python

    import datetime as dt
    import numpy as np
    from met_entities.RadolanRW import RadolanRW
    from met_entities.LonLatTime import LonLatTime

    # define target grid
    lon = np.arange(4, 15, .3)
    lat = np.arange(54, 48, -.3)
    lon_target, lat_target = np.meshgrid(lon, lat)
    lon_target = LonLatTime(data=lon_target)
    lat_target = LonLatTime(data=lat_target)

    # read, regrid, and export data
    start_datetime = dt.datetime(2022, 10, 25, 6, 50)
    radrw = RadolanRW()
    radrw.read_file(start_datetime=start_datetime, directory='path/to/data')
    radrw.regrid(lon_target=lon_target, lat_target=lat_target)
    radrw.export_netcdf('path/to/netcdf.nc')


Example: read/export multiple raw datasets
------------------------------------------

The next example shows the collection of data from multiple time steps. Here, only IconD2 data is read and
exported/appended to a netcdf file. At the first time step the netcdf file is created. All other steps will append data
to this file along the time dimension. The data is packed via ``data_format='i2'`` and ``scale_factor_nc=0.01`` (explanation
see below in :ref:`harmonize_example_label`).

.. code-block:: python

    import datetime as dt
    import numpy as np
    from met_entities.LonLatTime import LonLatTime
    from met_entities.IconD2 import IconD2

    # define target grid
    lon = np.arange(4, 15, .3)
    lat = np.arange(54, 48, -.3)
    lon_target, lat_target = np.meshgrid(lon, lat)
    lon_target = LonLatTime(data=lon_target)
    lat_target = LonLatTime(data=lat_target)

    # define forecast time steps to be read in
    start_datetime = dt.datetime(2022, 12, 14, 9)
    end_datetime = dt.datetime(2022, 12, 14, 15)
    period_step = dt.timedelta(hours=3)
    periods = [start_datetime]
    ct = 0
    while periods[-1] < end_datetime:
        ct = ct + 1
        periods.append(start_datetime + ct * period_step)

    # read, regrid, and export/append data
    nc_file = 'path/to/netcdf.nc'
    mode = 'create'
    for period in periods:
        print(period)
        icond2 = IconD2()
        icond2.read_file(period, directory='path/to/data', forecast_hours=2, fill_value=-1)
        icond2.regrid(lon_target=lon_target, lat_target=lat_target, file_nearest='path/to/regridding/rule.npz')
        if mode == 'create':
            icond2.export_netcdf(nc_file, data_format='i2', scale_factor_nc=0.01)
        else:
            icond2.export_netcdf_append(nc_file)
        mode = 'append'


.. _harmonize_example_label:

Example: read/harmonize/export multiple datasets of different types
-------------------------------------------------------------------

The class WeatherData is developed for the harmonizing of different data resources (the same temporal and spatial
resolution, the same scaling and the same filling of missing values). For instance, we start with measured RadolanRW
data from 2 hours ago until the most recent time step, include IconD2 forecast data and afterwards replace the nearest
forecast by nowcast products like RadolanRV. The latter two products are assumed to be stored in monthly separated
directories (e.g. .../icond2/202302/...). All summarized data shall be exported as a netcdf file with a filename
reflecting the actual datetime. The following code shows this example.

.. code-block:: python

    import numpy as np
    import datetime as dt

    from met_entities.VariableDescription import RegridDescription
    from met_entities.LonLatTime import LonLatTime
    from met_entities.WeatherData import WeatherData

    def main_combined_data():
      # define an arbitrary target grid
      lon = np.arange(4, 15, .3)
      lat = np.arange(54, 48, -.3)
      lon_target, lat_target = np.meshgrid(lon, lat)
      lon_target = LonLatTime(data=lon_target)
      lat_target = LonLatTime(data=lat_target)

      # define the time of the supposed last observation
      time_now = dt.datetime.now(tz=dt.timezone.utc)
      time_pivot = dt.datetime(time_now.year, time_now.month, time_now.day, time_now.hour, tzinfo=dt.timezone.utc)

      # define regridding descriptions for all used products
      regrid_description = {'radolanrw': RegridDescription(lon_target=lon_target, lat_target=lat_target,
                                                           file_nearest='data/radrw_regridding.npz'),
                            'radolanrv': RegridDescription(lon_target=lon_target, lat_target=lat_target,
                                                          file_nearest='data/radrv_regridding.npz'),
                            'icond2': RegridDescription(lon_target=lon_target, lat_target=lat_target,
                                                        file_nearest='data/icond2_regridding.npz')}

      # instantiate WeatherData class with central specifications (temporal/spatial resolution, filling, scaling, usage of
      # np.int16 type for memory saving)
      wd = WeatherData(time_now=time_pivot, delta_t=dt.timedelta(minutes=15), fill_value=-1, scale_factor=0.01,
                       regrid_description=regrid_description, short='int16')

      # collect all data, RadolanRW starting two hours ago
      wd.collect_radolanrw(time_start=time_pivot - dt.timedelta(hours=2), directory='path/to/radolanrw')
      wd.collect_icond2(latest_event=time_now, directory='path/to/icond2', dir_time_descriptor=['%Y%m'])
      wd.collect_radolanrv(latest_event=time_now, directory='/mnt/08_hwstore/RadolanRV', dir_time_descriptor=['%Y%m'])

      # export to netcdf with packed data (type i2: short integer, internal netcdf scaling of 0.01, undo the scaling from
      # above); the large data variables are compressed with zlib and compression level 4
      filename_nc = f'data/precipitation_data_example_{wd.time_now.strftime("%Y%m%d%H%M")}.nc'
      wd.export_netcdf(filename=filename_nc, institution='Institution as global attribute', data_format='i2',
                       scale_factor_nc=0.01, scale_undo=True, data_kwargs={'compression': 'zlib', 'complevel': 4})

Please note the ``scale_factor=0.01`` in the WeatherData instance in combination with ``short='int16'``. This combination
guarantees a precision of two floating points with a possible maximum of 327 (reasonable for 15 min precipitation). The
data could be left in original float precision (simply omit ``scale_factor`` and ``short``) if memory consumption is not
limiting. The possibility of ``short`` is specifically included for IconD2EPS data that consumes a lot more memory if used
for the whole prediction region and time.

Please also note the exporting to netcdf method with ``data_format='i2'``, ``scale_factor_nc=0.01``, and ``scale_undo=True``.
This combination takes back the scaling from data import and defines an internal scaling for netcdf and short integer
type for the data. The resulting netcdf file will be much smaller but ensures the same precision of two floating points.
This method is known as packing data values. All typical libraries to access netcdf content consider the ``scale_factor``
variable attribute and de-pack the data automatically. The accompanying ``add_offset`` attribute for packing is not
supported here as the data is typically extended over time and a perfect packing (recommendations here:
https://docs.unidata.ucar.edu/nug/current/best_practices.html) cannot be done. Moreover, it would introduce more
complexity with missing values.

Further the data variable is compressed to save storage.


Example: Download and collect Radolan data
------------------------------------------

To facilitate retrieving the data a download helper is incorporated in the small package ``download_data``. It enables
the download of recent data of currently RadolanRW, RadvorRQ, RadolanRV, IconD2, IconD2EPS, IconEU, and IconEUEPS data.
Please note, that nowcast and forecast data in particular typically do not stay on the website longer than 24 hours.

In the code block below the downloader is used to obtain the data from https://opendata.dwd.de/. Further a the data is
read in, cropped, and exported to netCDF. All time steps after the first one are appended to the existing netCDF file.

.. code-block:: python

    import datetime as dt
    from download_data.DownloadJob import DownloadJob
    from met_entities.RadolanRW import RadolanRW

    start_datetime = dt.datetime(2023, 4, 28, 2, tzinfo=dt.timezone.utc)
    end_datetime = dt.datetime(2023, 4, 28, 12, tzinfo=dt.timezone.utc)

    dj = DownloadJob(product='RadolanRW', directory='tmp', date_start=start_datetime, date_end=end_datetime)
    dj.download_files()

    period_step = dt.timedelta(hours=1)
    periods = [start_datetime]
    ct = 0
    while periods[-1] < end_datetime:
        ct = ct + 1
        periods.append(start_datetime + ct * period_step)
    mode = 'create'
    for period in periods:
        radrw = RadolanRW()
        radrw.read_file(start_datetime=period, directory='tmp', fill_value=-1)
        radrw.crop(lon_west=11.7, lon_east=15.2, lat_south=50.1, lat_north=51.8)
        if mode == 'create':
            radrw.export_netcdf('radrw_example.nc', data_format='i2', scale_factor_nc=0.1)
        else:
            radrw.export_netcdf_append('radrw_example.nc')
        mode = 'append'

    dj.delete_files()



Structure
---------

In the image below the UML class diagram of the weatherDataHarmonizer is depicted. A detailed version can be found in
https://gitlab.hrz.tu-chemnitz.de/ihmtud/weatherdataharmonizer/-/blob/main/weatherDataHarmonizer_UML.png.

.. image:: ../../weatherDataHarmonizer_UML_compact.png
    :alt: UML class diagram of weatherDataHarmonizer


License
-------

The weatherDataHarmonizer is licensed under Apache-2.0. You may obtain a copy
of the License in the project's LICENSE at https://gitlab.hrz.tu-chemnitz.de/ihmtud/weatherdataharmonizer/-/blob/main/LICENSE
or at http://www.apache.org/licenses/LICENSE-2.0.


Get Code and Contact
--------------------

The most recent version of weatherDataHarmonizer package can be obtained at https://gitlab.hrz.tu-chemnitz.de/ihmtud/weatherdataharmonizer/.
The author can be contacted via email: michael.wagner@tu-dresden.de.