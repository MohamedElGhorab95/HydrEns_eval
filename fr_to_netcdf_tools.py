
"""
Created on Mon Apr 10 14:44:46 2023

@author: M Elghorab
"""



def target(lat1, lat2, lon1, lon2, product, version):
    """
    # defining target extents and resolution based on the Radolan-RW raster

    Parameters
    ----------
    lat1 : minimum latitude.
    lat2 : maximum latitude.
    lon1 : minimum longtude.
    lon2 : maximum longitude.
    product : version of radar product.

    Returns
    -------
    lon_target : object containing target longitude data.
    lat_target : object containing target Latitude data.

    """
    import numpy as np
    from read_radolan import readRadolan
    # creating two matrices with coordinates similar to the product
    lon_target, lat_target = readRadolan.get_lonlat(version, product)
    # getting a matrix of booleans corressponding to datapoints covering the extents enclosed by the given cooedinates 
    idx = np.all((lon_target.data >= lon1, lon_target.data <= lon2,
                  lat_target.data >= lat1, lat_target.data <= lat2),
                 axis=0)
    # getting an array of indices along the vertical direction (axis=0) that are both TRUE and !=0
    idx_lon = np.any(idx, axis=0).nonzero()[0]
    # getting an array of indices along the horizontal direction (axis=1) that are both TRUE and !=0
    idx_lat = np.any(idx, axis=1).nonzero()[0]
    # modifying the matrices based on the index of TRUE values created before
    lon_target.data = lon_target.data[idx_lat[0]:idx_lat[-1], idx_lon[0]:idx_lon[-1]]
    lat_target.data = lat_target.data[idx_lat[0]:idx_lat[-1], idx_lon[0]:idx_lon[-1]]
    return lon_target, lat_target, (np.min(idx_lon), np.max(idx_lon)), (np.min(idx_lat), np.max(idx_lat))


# =========================================================================================================================


def timeframe(start_datetime, end_datetime, product):
    """
    this function creates a list of the time steps of the files (forecast / radar) to be read

    Parameters
    ----------
    start_datetime : iterable (list, tuple)
        iterable with the following format (YYYY, M, D,H).
    end_datetime : iterable (list, tuple)
        iterable with the following format (YYYY, M, D,H).
    product : str
        type of product to create the list for i.e. "forecast" or "radar".

    Returns
    -------
    times : List of datetime objects
        DESCRIPTION.

    """

    import datetime as dt

    if product == "forecast":

        # defining start times
        # define forecast time steps to be read in

        start_datetime = dt.datetime(start_datetime[0], start_datetime[1], start_datetime[2], start_datetime[3])
        end_datetime = dt.datetime(end_datetime[0], end_datetime[1], end_datetime[2], end_datetime[3])
        period_step = dt.timedelta(hours=3)
        periods = [start_datetime]
        ct = 0
        while periods[-1] < end_datetime:
            ct = ct + 1
            periods.append(start_datetime + ct * period_step)

        times = periods

    else:

        start_datetime = dt.datetime(start_datetime[0], start_datetime[1], start_datetime[2], start_datetime[3], 50)
        end_datetime = dt.datetime(end_datetime[0], end_datetime[1], end_datetime[2], end_datetime[3], 50)
        radolan_times = [start_datetime]
        hour = dt.timedelta(hours=1)
        ct = 0
        while radolan_times[-1] < end_datetime:
            ct += 1
            radolan_times.append(start_datetime + ct * hour)

        times = radolan_times

    return times


# =========================================================================================================================
def IconD2toNetCDF(ST, datafolder, longitude, latitude, nearestpoints, outputfile):
    """
    this function creates a netCDF file from IconD2 forecast files

    Parameters
    ----------
    ST : python dictionary.
        a dictionary containing forecast cycle starting times 
        with the following format {Year:YYYY,Month:M,Day:DD}
    datafolder : string
        Data folder location and name.
    forecast_issue: int or a python list
        [0,3,6,9,12,15,18]
    nearestpoints : string 
        file nearest .npz location and name.
    outputfile : string
        output file location and name .nc.
    longitude : container object
          target longitude object.
    latitude : container object
          target latitude object.

    Returns
    -------
    None.

    """

    from met_entities.IconD2 import IconD2

    # creating an instance of Icon object
    icond2 = IconD2()
    # looping over forecast cycle starting times
    for x in ST:
        # read Icon
        icond2.read_file(x, datafolder)
        # gridding the data
        icond2.regrid(lon_target=longitude, lat_target=latitude, file_nearest=nearestpoints)
        # exporting to netcdf if no file exists i.e. first forecast cycle
        if x == ST[0]:
            icond2.export_netcdf(filename=outputfile)
        # appending to the netcdf file
        else:
            icond2.export_netcdf_append(filename=outputfile)


# ====================================================================================================================================

def IconD2EPStoNetCDF(ST, datafolder, longitude, latitude, nearestpoints, outputfile):
    """
    this function creates a netCDF file from IconD2EPS forecast files

    Parameters
    ----------
   ST : python dictionary.
       a dictionary containing forecast cycle starting times 
       with the following format {Year:YYYY,Month:M,Day:DD}
    datafolder : string
        Data folder location and name.

    nearestpoints : string 
        file nearest .npz location and name.
    outputfile : string
        output file location and name .nc.
    longitude : container object
          target longitude object.
    latitude : container object
          target latitude object.

    Returns
    -------
    None.

    """
    from met_entities.IconD2EPS import IconD2EPS

    # creating an instance of IconEPS object
    icond2eps = IconD2EPS()
    # looping over forecast cycle starting times
    for x in ST:
        # read Icon ensemble files
        icond2eps.read_file(x, datafolder)
        # gridding the data
        icond2eps.regrid(lon_target=longitude, lat_target=latitude, file_nearest=nearestpoints)
        # exporting to netcdf if no file exists i.e. first forecast cycle
        if x == ST[0]:
            icond2eps.export_netcdf(filename=outputfile)
        else:
            icond2eps.export_netcdf_append(filename=outputfile)


# =================================================================================================================

def radolantoNetCDF(radolan_times, datafolder, idx_lon, idx_lat, outputfile):
    """
    this function creates a netCDF file from Radolan observation files

    Parameters
    ----------
    ST : python dictionary.
        a dictionary containing forecast cycle starting times 
        with the following format {Year:YYYY,Month:M,Day:DD}
    datafolder : string
        Data folder location and name.
    nearestpoints : string 
        file nearest .npz location and name.
    outputfile : string
        output file location and name .nc.
    idx_lon : tuple
          indexes of target max/min longitude.
    idx_lat : tuple
          indexes of target max/min latitude.

    Returns
    -------
    None.

    """
    from met_entities.RadolanRW import RadolanRW

    # creating an instance of RadolanRW object
    rad = RadolanRW()
    # looping over forecast cycle starting times
    for x in radolan_times:

        # read Icon ensemble files
        rad.read_file(start_datetime=x, directory=datafolder)
        # gridding the data
        rad.crop(idx_west=idx_lon[0], idx_east=idx_lon[1] - 1,
                 idx_south=idx_lat[1] - 1, idx_north=idx_lat[0])
        # rad.regrid(lon_target=lon, lat_target=lat,file_nearest=nearestpoints,neighbors=1)
        # exporting to netcdf if no file exists i.e. first forecast cycle
        if x == radolan_times[0]:
            rad.export_netcdf(filename=outputfile)
        else:
            rad.export_netcdf_append(filename=outputfile)

# =================================================================================================================


# ==============================================================================
###########################   Testing  ########################################
# ==============================================================================


if __name__ == '__main__':


# reading the files and exporting 

# getting all datapoint indices covering the extents of Radolan-RW raster for
# Sachsen after WGS84 > lat 50.1 - 51.8 | long 11.7 - 15.2 <
    lon, lat, id_lon, id_lat = target(50.1,51.8,11.7,15.2,"radolanrx",version=4)

    fortime = timeframe((2022,9,9,0), (2022,9,10,0), "forecast")

    radtime = timeframe((2022,9,5,0), (2022,9,10,0), "radar")


    # radolantoNetCDF(radtime, datafolder="//vs-grp08.zih.tu-dresden.de/hwstore/RadolanRW/202209",
    #                 idx_lon=id_lon, idx_lat=id_lat, 
    #                 outputfile="D:/Erasmus_FRM/05.Masterarbeit/03.Bearbeitung/02.netCDFs/RadolanRw/radRW_05_10_09_22.nc")
    
    # IconD2toNetCDF(ST=fortime,datafolder="//vs-grp08.zih.tu-dresden.de/hwstore/IconD2/tot_prec/202209",
    #                 longitude=lon,
    #                 latitude=lat,
    #                 nearestpoints="D:/Erasmus_FRM/05.Masterarbeit/03.Bearbeitung/01.Code/WorkspaceSachsen_nearestpoints.npz",
    #                 outputfile="D:/Erasmus_FRM/05.Masterarbeit/03.Bearbeitung/02.netCDFs/Icond2/icond2__05_10_09_22.nc")
    
    # IconD2EPStoNetCDF(fortime, "//vs-grp08.zih.tu-dresden.de/hwstore/IconD2eps/tot_prec/202209",
    #                   lon, lat,
    #                   "D:/Erasmus_FRM/05.Masterarbeit/03.Bearbeitung/01.Code/WorkspaceSachsen_nearestpoints.npz",
    #                   "D:/Erasmus_FRM/05.Masterarbeit/03.Bearbeitung/02.netCDFs/Icond2EPS/icond2eps__09_10_09_22.nc")