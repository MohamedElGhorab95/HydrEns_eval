"""
Created on Sun Apr  9 19:31:52 2023

@author: M Elghorab
"""
import copy
import pandas as pd
import geopandas as gp
import rioxarray as rio
from shapely.geometry import mapping
import shapefile
import pyproj
import xarray as xr
import numpy as np
import datetime as dt
import matplotlib.pyplot as plt
from shapely.geometry import mapping
import rasterio
from dask.distributed import Client
import pandas as pd

class Rainfall(object):

    def __init__(self, array):

        self.arr = array
        self.average = 0

    def rtrn_arr(self):
        """
        Returns
        -------
        array : xarray object

        """

        return self.arr

    def get_title(self):
        
        return self.rtrn_arr().name
        
        
    def aggr_spatial(self, scale_factor):
        """
        This method changes the spatial resolution of the Xarray type object.

        Parameters
        ----------

        scale_factor : int
            new grid cell size / old grid cell size.

        """
        
        # this prevents overwriting the original data object
        newself = copy.deepcopy(self)

        newself.arr = self.arr.coarsen(lat=scale_factor, lon=scale_factor, boundary="trim").mean()

        return newself

    def aggr_temporal(self, resolution):
        """
        This method changes the temporal resolution of the Xarray type object.

        Parameters
        ----------

        resolution : int
            hours to sum over.

        """
        
        # this prevents overwriting the original data object
        newself = copy.deepcopy(self)
        

        #newself.arr = newself.arr.resample(time="{}H".format(resolution)
         #                                  , loffset=pd.Timedelta(hours=resolution)).sum()

        vals = []
        times = []
        
        if type(self) == Ensemble_run:
            variable = newself.arr
            
        if type(self) == Deterministic_run:
            
            if 'average' in dir(self):
                if isinstance(self.average,int) == False:
                    variable = newself.average
                else:
                    variable = newself.arr
        else:
            variable = newself.arr
        
        for i in range(0,len(variable.time),resolution):
            try:
                times.append(variable.isel(time=i+resolution-1).time.values)
                vals.append(variable.isel(time=slice(i, i+resolution)).sum(dim='time'))
            except:
               break
        
        new_arr = xr.concat(vals, dim='time')
        new_arr['time'] = times

        # Assign the new DataArray to the modified object and rename it
        newself.arr = new_arr
        
        return newself

    def avg_areal_prec(self):
        """
        This method calculates The mean areal precipitation (MAP) over the range of 
        given pixels

        NOTE: this method only works with pre generated rainfall fields,
        however it can also work with ensemble object
        """
        
        # this prevents overwriting the original data object
        newself = copy.deepcopy(self)

        if (type(self) == Ensemble_run) or (type(self) == Deterministic_run):
            try:
                newself.average = newself.arr.mean(('lon', 'lat'))
            except:
                newself.average = newself.fr.mean(('lon', 'lat'))


        elif 'arr' in dir(self):

            newself.average = newself.arr.mean(('lon', 'lat'))

        return newself

    def extract_by_coords(self, lat_min, lat_max, lon_min, lon_max):
        """
        This method extracts smaller parts of the grid.
        NOTE:
            this method works with either the objects or the fields, however it is recommended
            to work with the objects for faster calculations

        Parameters
        ----------
        lat_min : float
            minimum latitude (bottom).
        lat_max : float
            maximum latitude (top).
        lon_min : float
            minimum longitude (left).
        lon_max : float
            maximum longitude (right).

        Returns
        -------
        Forecast/ observation Object
            rectangular grid bounded by the indicated coordinates.

        """

        
        # this prevents overwriting the original data object
        newself = copy.deepcopy(self)
        if 'gen_observation_field' in dir(self):
            rf = newself.gen_observation_field()
        elif 'gen_deterministic_field' in dir(self):
            rf = newself.gen_deterministic_field()
        else:
            rf = newself

        # get indexes of the coordinates

        # ilon = [rf.arr.indexes['lon'].get_loc(lon_min, method='nearest'),
        #         rf.arr.indexes['lon'].get_loc(lon_max, method='nearest')]
        # ilat = [rf.arr.indexes['lat'].get_loc(lat_min, method='nearest'),
        #         rf.arr.indexes['lat'].get_loc(lat_max, method='nearest')]

        ilon = rf.arr.indexes['lon'].get_indexer([lon_min, lon_max], method='nearest').tolist()
        ilat = rf.arr.indexes['lat'].get_indexer([lat_min, lat_max], method='nearest').tolist()

        # slicing the data array
        newself.arr = newself.arr.isel(lat=slice(ilat[1], ilat[0]), # lat is reversed to match the convention of xarray
                                       lon=slice(ilon[0], ilon[1]))

        return newself

    def extract_by_shp(self, shapefile):
        """
        This method extracts smaller parts of the grid based on a given shapefile.
        NOTE:
            this method works with either the objects or the fields, however it is recommended
            to work with the objects for faster calculations

        Parameters
        ----------
        lat_min : str
            path for mask shapefile.


        Returns
        -------
        Forecast/ observation Object
            rectangular grid bounded by the polygon given in the shapefile.

        """
        
         

        newself = copy.deepcopy(self)
        
        if isinstance(newself,Observation):
            if 'arr' not in dir(newself):
                rf = newself.gen_observation_field().arr
            else:
                rf = newself.arr
        
        if isinstance(newself,Deterministic_run):
            # if 'arr' not in dir(newself):
            #     rf = newself.gen_deterministic_field().arr
            # else:
                rf = newself.arr
            
            
        # if 'gen_observation_field' in dir(self):
        #     rf = newself.gen_observation_field().arr
        # elif 'gen_deterministic_field' in dir(self) and 'latitude' in dir(self.fr):
        #     rf = newself.gen_deterministic_field().fr
        # else:
        #     rf = newself.original

        # reading the shapefile
        da = rf.rio.write_crs("EPSG:4326")
        da = da.rename({"lat": 'y', 'lon': 'x'})

        shp = gp.read_file(shapefile)

        clipped = da.rio.clip(shp.geometry.apply(mapping), shp.crs)
        clipped = clipped.rename({"y": "lat", "x": "lon"})

        newself.arr = clipped
       
        return newself

    def shp_extents(self, shape_file):
        """
        This method gets the extents of a rectangle surrounding a shapefile.

        Parameters
        ----------
        shape_file : str
            relative path for the shapefile and file name WITHOUT .shp

        Returns
        -------
        lonmax : float
            eastmost point.
        lonmin : float
            westmost point.
        latmax : float
            northmost point.
        latmin : float
            southmost point.

        """

        
        if not shape_file.endswith(".shp"):
            
            # without extension
            shp = shape_file + ".shp"
            prj = shape_file + ".prj"
        else:
            shp = shape_file
            prj = shape_file.replace(".shp",".prj")

        # Open the shapefile
        sf = shapefile.Reader(shp)

        # determine the file projection
        with open(prj) as prj_file:
            prj_text = prj_file.read()
        crs = pyproj.CRS.from_string(prj_text).to_epsg()
        # Define the input projection
        input_proj = pyproj.Proj(init="EPSG:{}".format(crs))

        # Get the first polygon (assuming there is only one)
        polygon = sf.shapes()[0]

        # extract all vertices of the polygon
        lons = [lo[0] for lo in polygon.points]
        lats = [la[1] for la in polygon.points]

        # Get the latitudes of the northmost and southmost points
        latmax = max(lats)
        latmin = min(lats)

        # Get the longitudes of the easthmost and westhmost points

        lonmax = max(lons)
        lonmin = min(lons)

        # Define the output projection (WGS84)
        output_proj = pyproj.Proj(init='EPSG:4326')

        # Transform the coordinates
        lonmax, latmax = pyproj.transform(input_proj, output_proj, x=lonmax, y=latmax)
        lonmin, latmin = pyproj.transform(input_proj, output_proj, x=lonmin, y=latmin)

        return lonmax, lonmin, latmax, latmin

    def limit_to_shp(self, shape_file):
        """
        This method limits the calculation extents to the square 
        bounding the given shapefile.

        Parameters
        ----------
        shape_file : str
            relative path for the shapefile and file name WITHOUT .shp

        """

        import copy
        oper_self = copy.deepcopy(self)

        lonmax, lonmin, latmax, latmin = oper_self.shp_extents(shape_file)

        return oper_self.extract_by_coords(latmin - 0.00625, latmax + 0.00625, lonmin - 0.00625, lonmax + 0.00625)


    def limit_to_time(self, start_datetime, end_datetime):
        '''
        This method extracts a sub part of the data based on the specified time window.
        
        Parameters
        ----------
        start_datetime : iterable (list, tuple)
            iterable with the following format (YYYY, M, D,H).
        end_datetime : iterable (list, tuple)
            iterable with the following format (YYYY, M, D,H).

        '''
        
        oper_self = copy.deepcopy(self)
        
        strt = dt.datetime(start_datetime[0], start_datetime[1], start_datetime[2], start_datetime[3])
        endt = dt.datetime(end_datetime[0], end_datetime[1], end_datetime[2], end_datetime[3])
        
        oper_self.arr = oper_self.arr.sel(time=slice(strt, endt))
        
        return oper_self

    def to_raster(self, output_path_name):
        '''
        This method exports the data array to a raster file
        
        Parameters
        ----------
        output_path_name : str
            path and name of the file without ending with .tif
            .
        '''
        # set the profile for the output file
        profile = {
            'driver': 'GTiff',
            'dtype': self.arr.dtype,
            'nodata': 999999,
            'width': self.arr.sizes['lon'],
            'height': self.arr.sizes['lat'],
            'count': 1,
            'crs': 'EPSG:4326',
            'transform': rasterio.transform.from_bounds(self.arr.lon.values[0], 
                                            self.arr.lat.values[-1], 
                                            self.arr.lon.values[-1], 
                                            self.arr.lat.values[0], 
                                            self.arr.sizes['lon'], 
                                            self.arr.sizes['lat'])}
            # write the data array to the output file
        for t in range(len(self.arr.time.values)):
           # timestep for the file name
            n = self.arr.time.values[t] 
            with rasterio.open(output_path_name + "{}".format(str(n)[:-16]) + ".tif", 'w', **profile) as dst:
                dst.write(self.arr.isel(time=t).values, 1)
                
        return
        

# =================================================================================================================        
# =================================================================================================================
# =================================================================================================================


class Observation(Rainfall):

    def __init__(self, observation_file, horizon = None):
        self.file = observation_file
        self.wnd = horizon
        
        
        
    def gen_observation_field(self):
        """
        This method generates an observation xarray from the netCDF file and 
        limited to the indicated horizon



        """
        
        self.obs = xr.open_dataset(self.file, chunks='auto')
        
        
        # self.obs = self.obs.reindex(lat=self.obs['lat'][::-1])
               

        # adding the actual cooridnates to the xarray
        # creating a new list with the actual coordinates
       # # xarray flips the coordinates when opening the dataset
       # # max to min in the latitude to rearrange the data to its correct location
        l_lat = np.linspace(self.obs.latitude.max(), self.obs.latitude.min(), len(self.obs['lat'])).tolist()
        l_lon = np.linspace(self.obs.longitude.min(), self.obs.longitude.max(), len(self.obs['lon'])).tolist()
        # assigning the actual coordinates
        self.obs = self.obs.assign_coords(lat=l_lat, lon=l_lon)
        # drop unnecessary variables
        self.obs = self.obs.drop_vars(["longitude", "latitude"])

        # extracting the value variable to a dataXarray and structuring the array | observation
        self.arr = self.obs.radolanrw.rename("rainfall radar observation | Radolan-RW")
        self.arr = self.arr.assign_coords({"lat": self.obs.lat, "lon": self.obs.lon})

        # fixing the time coordinates from 00:50 to 01:00 by ceiling to the nearest 30 minutes
        self.arr['time'] = (self.arr['time'] + np.timedelta64(30, 'm')).astype('datetime64[h]')

        if self.wnd != None:        
            # limiting the observations to the chosen horizon
            self.arr = self.arr.isel(time=slice(0, self.wnd))
            
        self.arr = self.arr.load()
            
      

        return self

    
    def plot_all(self):
        
        array = self.rtrn_arr()
        
        for t in array.time:
            capture = array.sel(time=t)
            plt.figure()  # to avoid plotting figures on top of each other
            capture.plot(vmin=0, vmax=14.5)
            plt.title('Observation time\n {}'.format(str(t.values)[:-13]))

        return
    
    def plot(self, date_time):
        """
        This method plots the rainfall field for the extents of the data array 
        at the defined timestep

        

        Returns
        -------
        plot
            plot of data array.

        """

        

        # defining the time variable 
        window = dt.datetime(date_time[0], date_time[1], date_time[2], date_time[3])

        # ==========================================================================
        # plot title building
        # extract the radar time
        # array = self.gen_observation_field().rtrn_arr()
        array = self.rtrn_arr()
        radt = array.time.values
        # convert from numpy datetime array to list with the date as a string
        df = pd.DataFrame({'datetime': radt})
        df.set_index('datetime', inplace=True)
        df.index = df.index.strftime('%Y/%m/%d %H:%M')
        # ==========================================================================

        # capture the frame for the required time step 
        capture = array.sel(time=window)
        plt.figure()  # to avoid plotting figures on top of each other
        capture.plot(vmin=0, vmax=14.5)
        plt.title('Observation time\n {}'.format(df.index[window.hour - 1]))

        return


    def ident_event(self, threshold):
        '''
        This method returns the timesteps where the observation has a value
        that exceeds (or equals to) a certain given threshold

        Parameters
        ----------
        threshold : float
            rainfall threshold in mm/hr.

        Returns
        -------
        formatted_time : list
            list of pandas time stamps.

        '''
        
        # create a boolean array where values > threshold are True and others are False
        ev_trfl = self.arr >= threshold
        times = []
        formatted_time = []
        for t in ev_trfl.time:
            if True in ev_trfl.sel(time =t):
                times.append(t.values)
        for t in times:
            
            ft = pd.to_datetime(str(t)) 
            formatted_time.append(ft)
            
            
        return formatted_time

# =================================================================================================================        
# =================================================================================================================
# =================================================================================================================


class Forecast(Rainfall):

    def __init__(self, forecast_file, forecast_cycle=None, horizon=None):

        
        self.fr = xr.open_dataset(forecast_file)
        self.cyc = forecast_cycle
        self.wnd = horizon
        self.dssid = 0
        
        self.original = self.gen_dataset()
        
    def gen_dataset(self):
        """
        This method generates the dataset built on the indicated forecast cycle and horizon.

        Returns
        -------
        xarray dataset
            DESCRIPTION.

        """

        
        # # adding the actual cooridnates to the xarray
        # # creating a new list with the actual coordinates
        # l_lat = np.linspace(self.fr.latitude.max(), self.fr.latitude.min(), len(self.fr['lat'])).tolist()
        # l_lon = np.linspace(self.fr.longitude.min(), self.fr.longitude.max(), len(self.fr['lon'])).tolist()
        # # assigning the actual coordinates
        # self.fr = self.fr.assign_coords(lat=l_lat, lon=l_lon)
        # # drop unnecessary variables
        # self.fr = self.fr.drop_vars(["longitude", "latitude"])
        
        if 'start_time' in self.fr.dims: # should always give false, except for the case of leadtime netcdfs which are already preset
            # extracting the forecast cycle data
            # adding the actual cooridnates to the xarray
            # creating a new list with the actual coordinates
            l_lat = np.linspace(self.fr.latitude.max(), self.fr.latitude.min(), len(self.fr['lat'])).tolist()
            l_lon = np.linspace(self.fr.longitude.min(), self.fr.longitude.max(), len(self.fr['lon'])).tolist()
            # assigning the actual coordinates
            self.fr = self.fr.assign_coords(lat=l_lat, lon=l_lon)
            # drop unnecessary variables
            self.fr = self.fr.drop_vars(["longitude", "latitude"])
            
            self.fr = self.fr.isel(time=self.cyc - 1)
            # limiting the forecast to the chosen horizon
            self.fr = self.fr.isel(
                forecast_time=slice(0, self.wnd * 4))  # multiply by 4 because for each hour 4 available timesteps
            # skipping the first zero value
            self.fr = self.fr.isel(forecast_time=slice(1, len(self.fr.forecast_time)))
            # temporal aggregation from 15 minutes to one hour for the icond2 values
            self.fr = self.fr.resample(forecast_time="1H").sum()  # mm/15min -> mm/hr
            # correcting dimension name
            self.fr = self.fr.rename({"time": "start_time", "forecast_time": "time"})
    
            # create a list of  date time objects  ending with the forecast time
            time = [dt.timedelta(hours=i) for i in np.arange(1, len(self.fr.time) + 1).tolist()]
            # assigning new times to solve the 1:15, 2:15,.... time problem
            self.fr = self.fr.assign_coords({'time': time})
            # referencing the time coordinate based on the forecast release time as its starting point
            self.fr = self.fr.assign_coords({'time': self.fr.time + self.fr.start_time})
        # else:
        #     l_lat = np.linspace(51.921707, 49.991165, 212).tolist()
        #     l_lon = np.linspace(11.634748, 15.226268, 253).tolist()
        #     # assigning the actual coordinates
        #     self.fr = self.fr.assign_coords(lat=l_lat, lon=l_lon)
        #     # drop unnecessary variables
        #     self.fr = self.fr.drop_vars(["longitude", "latitude"])

        
        return self.fr.load()

    def gen_leadtime(lt, events):
        pass



    def preplot(self, member):

        # if isinstance(self, Ensemble_run):
        if 'gen_quantiles' in dir(self):

            array = xr.DataArray(self.fr.variables['icond2eps_{}'.format(member-1)], coords=self.fr.coords)
        else:

            array = self.gen_deterministic_field().rtrn_arr()

        return array

    def plot(self, date_time, Ens_mem=None):
        """
        This method plots the rainfall field for the extents of the data array 
        at the defined timestep

        Parameters
        ----------

        date_time : iterable
            exact date and time to be plotted
            in the format of (YYYY,M,D,H).

        Returns
        -------
        plot
            plot of data array.

        """

       

        # defining the time variable 
        window = dt.datetime(date_time[0], date_time[1], date_time[2], date_time[3])

        # ==========================================================================
        # plot title building
        # extract the radar time

        # to be able to excute the plot command multiple times without generating again
        try:
            array = self.arr
        except:
            array = self.preplot(Ens_mem)

        # capture = data_array.sel(time = dt.timedelta(hours=plot_date.hour))  
        capture = array.sel(time=window)
        # extract forecast release time
        if Ens_mem == None:
            release = array.start_time.values
            formatted_releasedate = np.datetime_as_string(release, unit='s').replace('T', ' ').replace('-', '/')

        # ==========================================================================    
        # Creating the plot

        plt.figure()  # to avoid plotting figures on top of each other
        capture.plot(vmin=0, vmax=14.5)
        if Ens_mem == None:
            plt.title('Date/ Time: {}  \n Forecast released: {}'.format(window, formatted_releasedate))
        else:
            plt.title('Date/ Time: {} \nEnsemble member: {}  '.format(window, Ens_mem))
        return
    
    
    def plot_all(self, averaging_method=None):
        
        try:
            array = self.arr
        except:
            array = self.preplot(averaging_method)
        
        # extract forecast release time
        release = array.start_time.values
        formatted_releasedate = np.datetime_as_string(release, unit='s').replace('T', ' ').replace('-', '/')

        
        for t in array.time:
            capture = array.sel(time=t)
            plt.figure()  # to avoid plotting figures on top of each other
            capture.plot(vmin=0, vmax=14.5)
            plt.title('Date/ Time: {}  \n Forecast released: {}'.format(str(t.values)[:-13],formatted_releasedate))

        return


# =================================================================================================================


class Deterministic_run(Forecast):

    def __init__(self, forecast_file, forecast_cycle=None, horizon=48):
        Forecast.__init__(self, forecast_file, forecast_cycle=None, horizon=48)

    def gen_deterministic_field(self):
        """
        This method generates the datarray built on the indicated forecast cycle and horizon.


        """
        

        # if self.dssid == 0:
        #     ds = self.gen_dataset()
        # else:
        #     ds = self.fr
        ds = self.fr
        # extracting the variable name
        name = list(ds.data_vars.keys())[0]
        # extracting the value variable to a dataXarray and renaming it       
        self.arr = xr.DataArray(ds[name])
        self.arr = self.arr.rename("rainfall forecast | deterministic run")

        self.dssid = 1
        
        
        return self


# =================================================================================================================


class Ensemble_run(Forecast):

    def __init__(self, forecast_file, forecast_cycle=None, horizon=48):
        Forecast.__init__(self, forecast_file, forecast_cycle=None, horizon=48)
        self.dssid = 2

    

    def get_ensemble_members(self):
        
        return self.fr

    def gen_quantiles(self, averaging):
        """
        This method generates the datarray built on the indicated forecast cycle and horizon.

        Parameters
        ----------


        averaging : str / int
            numerical method to obtain a descriptive value for the ensemble members,
            if left empty a dataset with all ensemble members is returned,
            if averaging = 'mean'; the arithmatic mean of the ensemble members is returned
            if averaging = 'median'; the median of the ensemble members is returned
            if averaging = int between 0 and 100; the percentile of the ensemble members is returned.


        """
        

        # check if the members are areal averages or spatially distribured
        # if len(self.fr.dims) > 1 and self.dssid != 2.1:
        #     ds = self.gen_dataset().fr
        # else:
        #     ds = self.fr
        
                
        # newself = copy.deepcopy(self)
        
        # if self.average ==0:
        #     ds = copy.deepcopy(self.original)
        # else:
        #     ds = self.average
            
        # # parallelizing the computation for acceleration
        

        # # # Set up a Dask distributed client
        # # client = Client()
            
        # # Enable Dask and chunk the dataset along desired dimensions
        # ds = ds.chunk({'time': -1})  # Chunking along the 'time' and 'new' dimensions
        
        
        # # if isinstance(self.average,int) == False:
        # # # if self.average !=0:
        # #     ds = self.average
        # # else:
        # #     ds = copy.deepcopy(self.original)
        
        
        
        # if averaging == 'mean':
        #     # calculate ensemble mean
        #     mean = ds.to_array(dim='new').mean(dim='new', skipna=True).compute()
        #     # add the mean to the dataset
        #     ds = ds.assign(ensemble_val=mean)
        #     # extracting the mean variable to a dataXarray and structuring the array | forecast
        #     newself.arr = ds.ensemble_val.rename("EPS rainfall forecast | mean")
        # elif averaging == 'median':
        #     # calculate ensemble median
        #     median = ds.to_array(dim='new').median(dim='new', skipna=True).compute()
        #     # add the median to the dataset
        #     ds = ds.assign(ensemble_val=median)
        #     # extracting the median variable to a dataXarray and structuring the array | forecast
        #     newself.arr = ds.ensemble_val.rename("EPS rainfall forecast | median")
        # else:
        #     # calculate ensemble percentile
        #     per = ds.to_array(dim='new')  # extract the ensembles to a new dimension and put it a new array
        #     # calculate the percentile over this dimension
        #     per = per.quantile(q=averaging / 100, dim='new', skipna=True).compute()

        #     # assing the ensemble variable in the original data set to the created array
        #     ds = ds.assign(ensemble_val=per)
        #     # extracting the mean variable to a dataXarray and renaming it
        #     newself.arr = ds.ensemble_val.rename("EPS rainfall forecast | {}th percentile".format(averaging))
        
        
        
        newself = copy.deepcopy(self)
        
               
        try:
            if isinstance(self.average,int) == False:
            # if self.average !=0:
                ds = self.average
        except:
            ds = copy.deepcopy(self.original)
        
        # Enable Dask and chunk the dataset along desired dimensions
        ds = ds.chunk({'time': -1})  # Chunking along the 'time' and 'new' dimensions
        
        if averaging == 'mean':
            # calculate ensemble mean
            mean = ds.to_array(dim='new').chunk({'new': -1}).mean(dim='new', skipna=True).compute()
            # add the mean to the dataset
            ds = ds.assign(ensemble_val=mean)
            # extracting the mean variable to a dataXarray and structuring the array | forecast
            newself.arr = ds.ensemble_val.rename("EPS rainfall forecast | mean").load()
        elif averaging == 'median':
            # calculate ensemble median
            median = ds.to_array(dim='new').chunk({'new': -1}).median(dim='new', skipna=True).compute()
            # add the median to the dataset
            ds = ds.assign(ensemble_val=median)
            # extracting the median variable to a dataXarray and structuring the array | forecast
            newself.arr = ds.ensemble_val.rename("EPS rainfall forecast | median").load()
        else:
            # calculate ensemble percentile
            per = ds.to_array(dim='new').chunk({'new': -1})  # extract the ensembles to a new dimension and put it a new array
            # calculate the percentile over this dimension
            per = per.quantile(q=averaging / 100, dim='new', skipna=True).compute()

            # assing the ensemble variable in the original data set to the created array
            ds = ds.assign(ensemble_val=per)
            # extracting the mean variable to a dataXarray and renaming it
            newself.arr = ds.ensemble_val.rename("EPS rainfall forecast | {}th percentile".format(averaging)).load()
        
        
        return newself

    

    def get_ens_quantiles(self):
        
        return self.arr

    def eps_accelerate_by_shp(self, shape_file):
        """
        This method accelerates the ensemble calculations by limiting
        the extents to a square around the shapefile extents

        NOTE:
            this does not limit the calculations to the inside of the shapefile,
            to do that execute eps_extract_by_shp() after acceleration and 
            ensemble field generation.


        Parameters
        ----------
        shape_file : str
            relative path for the shapefile and file name WITHOUT .shp

        """

        

        # this prevents overwriting the original data object
        oper_self = copy.deepcopy(self)

        lonmax, lonmin, latmax, latmin = oper_self.shp_extents(shape_file)

        return oper_self.eps_extract_by_coords(latmin - 0.00625, latmax + 0.00625, lonmin - 0.00625, lonmax + 0.00625)

    def eps_extract_by_coords(self, lat_min, lat_max, lon_min, lon_max):
        """
        This method extracts smaller parts of the grid

        Parameters
        ----------
        lat_min : float
            minimum latitude (bottom).
        lat_max : float
            maximum latitude (top).
        lon_min : float
            minimum longitude (left).
        lon_max : float
            maximum longitude (right).

        Returns
        -------
        Forecast Object
            rectangular grid bounded by the indicated coordinates.

        """
        

        # this prevents overwriting the original data object
        # oper_self = copy.deepcopy(self)
        # ds = oper_self.gen_dataset().fr
        # # ds = oper_self.fr
        newself = copy.deepcopy(self)
        
        oper_self = copy.deepcopy(self.original)

        # get indexes of the coordinates        
        ilon = oper_self.indexes['lon'].get_indexer([lon_min, lon_max], method='nearest').tolist()
        ilat = oper_self.indexes['lat'].get_indexer([lat_min, lat_max], method='nearest').tolist()

        # slicing the data array

        newself.fr = oper_self.isel(lat=slice(ilat[1], ilat[0]),  # reversed to match the xarray convention
                                                lon=slice(ilon[0], ilon[1]))

        # newself.dssid = 2.1
        return newself

    def eps_extract_by_shp(self, shapefile):

        
        # this prevents overwriting the original data object
        oper_self = copy.deepcopy(self.original)
        # ds = oper_self.gen_dataset()
        # ds = oper_self.arr
        ds = oper_self
        newself = copy.deepcopy(self)

        da = ds.rio.write_crs("EPSG:4326")
        da = da.rename({"lat": 'y', 'lon': 'x'})

        # reading the shapefile
        shp = gp.read_file(shapefile)

        clipped = da.rio.clip(shp.geometry.apply(mapping), shp.crs)
        clipped = clipped.rename({"y": "lat", "x": "lon"})

        newself.fr = clipped
        # newself.fr = clipped

        # newself.dssid = 2.1

        return newself

# =================================================================================================================        
# =================================================================================================================
# =================================================================================================================


class leadtime(object):
    def __init__(self, forecast_file, forecast_type ,leadtime, events):
        
        self.lt = leadtime
        self.events = events
        self.typ = forecast_type
    
    def gen(self):
        if self.typ == "deterministic":
            pass

# =================================================================================================================        
# =================================================================================================================
# =================================================================================================================


class Runoff(object):

    def __init__(self, forecast_file):
        
        self.fr = xr.open_dataset(forecast_file)
         
       
    def aggr_temporal(self, resolution):
        """
        This method changes the temporal resolution of the Xarray type object.

        Parameters
        ----------

        resolution : int
            hours to sum over.

        """
        
        # this prevents overwriting the original data object
        newself = copy.deepcopy(self)
        

        #newself.arr = newself.arr.resample(time="{}H".format(resolution)
         #                                  , loffset=pd.Timedelta(hours=resolution)).sum()

        vals = []
        times = []
        
        variable = newself.fr
        
        for i in range(0,len(variable.timeQ),resolution):
            try:
                times.append(variable.isel(timeQ=i+resolution).timeQ.values)
                vals.append(variable.isel(timeQ=slice(i, i+resolution)).sum(dim='timeQ'))
            except:
               break
        
        new_arr = xr.concat(vals, dim='timeQ')
        new_arr['timeQ'] = times

        # Assign the new DataArray to the modified object and rename it
        newself.fr = new_arr
        newself.fr = newself.fr.rename(timeQ='time')
        
        return newself

# =================================================================================================================        
# =================================================================================================================
# =================================================================================================================


class R_Observation(Runoff):

    def __init__(self, observation_file):
        # Read the CSV file using pandas
        df = pd.read_csv(observation_file)
        df['begin'] = pd.to_datetime(df['begin']).dt.to_pydatetime()
        # df = df.rename(columns={'begin': 'time', 'pikobytes$hwims$550940$q-ziel-tw-15m_hi': 'runoff_discharge'})
        df = df.rename(columns={'begin': 'time', 'hi': 'runoff_discharge'})
        
        # # Convert the DataFrame to an xarray DataArray
        self.fr = xr.DataArray(data=df["runoff_discharge"], coords=[df["time"]], dims=["time"])

        
        
    # def gen_obs(self):
        
    #     import copy
    #     newself = copy.deepcopy(self)
    #     ds = newself.fr
    #     newself.fr = ds.Q_obs_value
    #     return newself


# =================================================================================================================        
# =================================================================================================================
# =================================================================================================================


class R_Forecast(Runoff):

    def __init__(self, forecast_file):
        Runoff.__init__(self, forecast_file)
        
       
        # find the common time steps between timeQobs and timeQ
        common_times = set(self.fr['timeQ'].values).intersection(self.fr['timeQobs'].values)
        forecast_time = [time for time in self.fr['timeQ'].values if time not in common_times]
        
        self.fr["Q_for_mean_value"]= self.fr["Q_for_mean_value"].sel(timeQ=forecast_time)
        self.fr["Q_q10"]= self.fr["Q_q10"].sel(timeQ=forecast_time)
        self.fr["Q_q25"]= self.fr["Q_q25"].sel(timeQ=forecast_time)
        self.fr["Q_q50"]= self.fr["Q_q50"].sel(timeQ=forecast_time)
        self.fr["Q_q75"]= self.fr["Q_q75"].sel(timeQ=forecast_time)
        self.fr["Q_q90"]= self.fr["Q_q90"].sel(timeQ=forecast_time)
        
        for e in range(20):
            self.fr['eps_{}'.format(e)] = self.fr.Q_for_eps_value[e].sel(timeQ=forecast_time)
       
        
    def gen_quantiles(self, averaging):
        
        
        newself = copy.deepcopy(self)
        ds = newself.fr
        
        if averaging == 'mean':
            # extracting the mean variable to a dataXarray and structuring the array | forecast
            ds['Q_for_mean_value'] = ds.Q_for_mean_value.rename({'timeQ':'time'})
            
            newself.fr = ds.Q_for_mean_value.dropna(dim='time')
        
        else:
            # calculate ensemble percentile
            per = "Q_q{}".format(averaging)  # extract the ensembles to a new dimension and put it a new array
            
            # extracting the mean variable to a dataXarray and renaming it
            ds[per] = ds[per].rename({'timeQ':'time'})
            
            newself.fr = ds[per].dropna(dim='time')
        
        return newself
    
    def gen_ensmembers(self):
        
        newself = copy.deepcopy(self)
        ds = newself.fr
        
        da = xr.Dataset()
        
        for e in range(20):
            
            da['eps_{}'.format(e)] = ds['eps_{}'.format(e)].dropna('timeQ')
          
            
        da = da.rename_dims({'timeQ':'time'})
        da = da.rename_vars({'timeQ': 'time'})
        newself.fr = da
        
        return newself
        
        


# ==============================================================================
###########################   Testing | Examples  #############################
# ==============================================================================


if __name__ == '__main__':

    def test_Obs():
        """
        Testing the Observation Class

        """

        # generate an observation object
        rad = Observation("Data/NetCDFs/Example/radRW_example.nc", 48)
        # generating a rainfall field
        # only needed for manual low level operations, however if the object is to be 
        # used in other classes (Cont, ROC) it gets generated automatically
        rad.gen_observation_field()
        # aggregating the rainfall field
        rad_agg = rad.aggr_spatial(3)
        # extracting a part of the observation using geographical coordinates
        sub_rad = Observation("Data/NetCDFs/Example/radRW_example.nc", 48).extract_by_coords(50.02, 51.01, 11.66, 12.68)
        # mean areal rainfall
        rad_avg = rad.avg_areal_prec()

        # Weise Elster
        weise = rad.extract_by_shp("shp/WeiseElster/OelsnitzTEZG_DHDN.shp")

        rad_3 = rad.aggr_temporal(3)

        rad.plot([2021,10,4, 23])
        rad_agg.plot([2021,10,4, 23])
        sub_rad.plot([2021,10,4, 23])

        for i in range(1, 23):
            weise.plot([2021,10,4, i])

        return rad, rad_agg, sub_rad, rad_avg, weise, rad_3


    # rad, rad_agg, sub_rad, rad_avg, weise, rad_3 = test_Obs()

    def test_Det():
        """
        Testing the Deterministic run Class

        """

        # generate the forecast object
        icond2 = Deterministic_run("Data/NetCDFs/Example/icond2_example.nc", 1, 12)
        # generating a rainfall field
        # only needed for manual low level operations, however if the object is to be 
        # used in other classes (Cont, ROC) it gets generated automatically
        icond2.gen_deterministic_field()
        # aggregating the rainfall field
        icon_agg = icond2.aggr_spatial(2)
        # extracting a part of the observation using geographical coordinates
        sub_icon = Deterministic_run("Data/NetCDFs/Example/icond2_example.nc", 1, 12).limit_to_shp("shp/Mugliz/mugliz_cats")
        # mean areal rainfall
        icon_avg = icond2.avg_areal_prec()

        mugliz = icond2.extract_by_shp("shp/Mugliz/mugliz_cats.shp")

        mandau = icond2.extract_by_shp("shp/Mandau/egp6741491_Zittau_5_TEZG_neu.shp")

        ico_3 = icond2.aggr_temporal(4)

        icond2.plot([2021, 7, 14, 7])
        icon_agg.plot([2021, 7, 14, 1])
        sub_icon.plot([2021, 7, 14, 1])

        for i in range(1, 12):
            mugliz.plot([2021, 7, 14, i])
            mandau.plot([2021, 7, 14, i])

        return icond2, icon_agg, sub_icon, icon_avg, ico_3


    # icond2, icon_agg, sub_icon, icon_avg, ico_3 = test_Det()

    def test_EPS():
        """
        Testing the Ensemble run Class

        """

        # generate the forecast object
        eps = Ensemble_run("Data/NetCDFs/Example/icond2eps_example.nc", 1)
        # generating a rainfall field for a specific quantile
        # only needed for manual low level operations, however if the object is to be 
        # used in other classes (Cont, ROC) it gets generated automatically
        # eps95 = eps.gen_quantiles(95)  # 95th quantile      
        # epsavg = eps.gen_quantiles("mean")  # mean
        # epsmed = eps.gen_quantiles("median")  # median

        # extracting a part of the forecast using geographical coordinates
        # the method eps_accelerate_by_shp is used instead of the general extract_by_coords
        # to speed up the calculation process for the smaller areas
        sub95 = eps.eps_accelerate_by_shp("shp/WeiseElster/OelsnitzTEZG_DHDN").gen_quantiles(95)
        subavg = eps.eps_accelerate_by_shp("shp/WeiseElster/OelsnitzTEZG_DHDN").gen_quantiles("mean")
        submed = eps.eps_accelerate_by_shp("shp/WeiseElster/OelsnitzTEZG_DHDN").gen_quantiles("median")

        # mean areal rainfall

        m_eps95 = sub95.avg_areal_prec()
        m_epsavg = subavg.avg_areal_prec()
        m_epsmed = submed.avg_areal_prec().rtrn_arr()

        sub_agg = sub95.aggr_spatial(3).aggr_temporal(3)

        # # the plot method handles plotting ensemble objects,
        # however it is recommended to generate the field first then plot
        # eps95.plot([2021,10,4, 23])
        # epsavg.plot([2021,10,4, 23])
        # epsmed.plot([2021,10,4, 23])
        sub95.plot([2021,10,4, 23])
        subavg.plot([2021,10,4, 23])
        submed.plot([2021,10,4, 23])
        sub_agg.plot([2021,10,4, 23])

        return  m_eps95, m_epsavg, m_epsmed, sub95, subavg, submed, eps
        
    
    
    m_eps95, m_epsavg, m_epsmed, sub95, subavg, submed, eps = test_EPS()
    