# -*- coding: utf-8 -*-
"""
Created on Sat May  6 13:58:28 2023

@author: M Elghorab
"""


import time
import os
from engine import fr_to_netcdf_tools
from datetime import datetime, timedelta






# reading the files and exporting 

def modify_date_dict(dict_A):
    dict_B = {}

    for key, value in dict_A.items():
        start_date_A, end_date_A = value
        start_date_A = datetime(*start_date_A)
        end_date_A = datetime(*end_date_A)

        start_date_B = start_date_A - timedelta(days=1)
        end_date_B = end_date_A 

        dict_B[key] = ((start_date_B.year, start_date_B.month, start_date_B.day, 0),
                       (end_date_B.year, end_date_B.month, end_date_B.day, 0))

    return dict_B

def event_file(source_file):
    '''
    this function reads the text file containing the event start and end dates

    Parameters
    ----------
    source_file : str
        path to txt file.

    Returns
    -------
    lib : dict
        dictionary of all events in the text file.

    '''
    with open(source_file,'r') as f:
        # lines = f.readlines()[1:]
        alist = [line.rstrip() for line in f]
        lines = alist[1:]
    ################################################################################################
    # creating a dictionary of the events
    events = [[],[]]
    evid = [idi[0:2] for idi in lines]  
    for l in lines:
        
        events[0].append((int((l[9:13])),int((l[6:8])),int((l[3:5])),0))
        
        events[1].append((int(float(l[20:24])),int(float(l[17:19])),int(float(l[14:16])),0))

    # event library dictionary
    lib = dict(zip(evid,zip(events[0],events[1])))

    # garbage collection
    del f, l, lines, events, evid
    return lib






def create_radar(file):
    
    ###############################################################################################
    lib = event_file(file)
    t = []
    # iterate over the event library
    for ev in lib:

        radtime = fr_to_netcdf_tools.timeframe(lib[ev][0], lib[ev][1], "radar")
        for i in radtime:
            # creating a list of radolan timestamps to be loaded
            t.append(i)
        
       
        #=============================================================================================================================
    if 'Cos' in file:
        out = 'radRW_COS.nc'
    elif 'Ico' in file:
        out = 'radRW_ICO.nc'
        
    fr_to_netcdf_tools.radolantoNetCDF(t, datafolder="//vs-grp08.zih.tu-dresden.de/hwstore/RadolanRW/",
                    idx_lon=id_lon, idx_lat=id_lat, 
                    # outputfile="//vs-grp07.zih.tu-dresden.de/howa/work/students/Mohamed_Elghorab/netCDFs/radRW_ICO.nc")
                    outputfile="{}/Data/{}".format(os.getcwd(),out))
        
       
    return t, lib



def create_forecast(file, product, mode= None):
    
    lib = event_file(file)
    lib = modify_date_dict(lib)
    ###############################################################################################
    x = [] 
    # iterate over the event library
    for ev in lib:
       
        fortime = [fr_to_netcdf_tools.timeframe(lib[ev][0], lib[ev][1], "forecast")]
        for t in fortime:
            for z in t:
                # creating a list of icond2 timestamps to be loaded
                x.append(z)
       
        #=============================================================================================================================
    
    if product == 'ICON':

        if mode == "D":
            fr_to_netcdf_tools.IconD2toNetCDF(x, "//vs-grp08.zih.tu-dresden.de/hwstore/IconD2/tot_prec/",
                                              
                            lon, lat,
                            "{}/in/Sachsen_nearestpoints.npz".format(os.getcwd()),
                            "{}/Data/icond2_ev.nc".format(os.getcwd()))
        elif mode == "eps":
            fr_to_netcdf_tools.IconD2EPStoNetCDF(x, "//vs-grp08.zih.tu-dresden.de/hwstore/IconD2eps/tot_prec/",
                              lon, lat,
                              "{}/in/Sachsen_nearestpoints.npz".format(os.getcwd()),
                              "{}/Data/icond2eps_ev.nc".format(os.getcwd()))
    elif product == 'COSMO':
        
        if mode == "D":
            fr_to_netcdf_tools.CosmoD2toNetCDF(x, "//vs-grp08.zih.tu-dresden.de/hwstore/CosmoD2/",
                            lon, lat,
                            "{}/in/Sachsen_nearestpoints_cosmo.npz".format(os.getcwd()),
                            "{}/Data/cosmod2_ev.nc".format(os.getcwd()))
        elif mode == "eps":
            fr_to_netcdf_tools.CosmoD2EPStoNetCDF(x, "//vs-grp08.zih.tu-dresden.de/hwstore/CosmoD2eps/",
                              lon, lat,
                              "{}/in/Sachsen_nearestpoints_cosmo.npz".format(os.getcwd()),
                              "{}/Data/cosmod2eps_ev.nc".format(os.getcwd()))
    else:
        print('product should be either ICON or COSMO')
  
    
    
    
################################################################################################

if __name__ == '__main__':
################################################################################################
# getting all datapoint indices covering the extents of Radolan-RW raster for
# Sachsen after WGS84 > lat 50.1 - 51.8 | long 11.7 - 15.2 <
    
    lon, lat, id_lon, id_lat = fr_to_netcdf_tools.target(50.1,51.8,11.7,15.2,"radolanrx",version=4)

################################################################################################

    tic = time.time()
    create_radar('in/Cosmo_events.txt')
    print((time.time()-tic)/3600) # operation time in hours
    
    
    tic = time.time()
    create_forecast('in/Cosmo_events.txt','COSMO',"D")
    print((time.time()-tic)/3600) # operation time in hours
    
    
    tic = time.time()
    create_forecast('in/Cosmo_events.txt','COSMO',"eps")
    print((time.time()-tic)/3600) # operation time in hours
    
    
    tic = time.time()
    create_radar('in/Icon_events.txt')
    print((time.time()-tic)/3600) # operation time in hours
    
    
    tic = time.time()
    create_forecast('in/Icon_events.txt','ICON',"D")
    print((time.time()-tic)/3600) # operation time in hours
    
    
    tic = time.time()
    create_forecast('in/Icon_events.txt','ICON',"eps")
    print((time.time()-tic)/3600) # operation time in hours

