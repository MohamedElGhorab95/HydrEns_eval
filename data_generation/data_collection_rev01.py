# -*- coding: utf-8 -*-
"""
Created on Sat May  6 13:58:28 2023

@author: M Elghorab
"""


import time
from engine import fr_to_netcdf_tools

# reading the files and exporting 


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
   
        
    fr_to_netcdf_tools.radolantoNetCDF(t, datafolder="//vs-grp08.zih.tu-dresden.de/hwstore/RadolanRW/",
                    idx_lon=id_lon, idx_lat=id_lat, 
                    # outputfile="//vs-grp07.zih.tu-dresden.de/howa/work/students/Mohamed_Elghorab/netCDFs/radRW_ICO.nc")
                    outputfile="D:/Erasmus_FRM/05.Masterarbeit/03.Bearbeitung/01.Code/radRW_cosmod2eps.nc")
        
       
    return t, lib



def create_icond2(file, mode= None):
    
    lib = event_file(file)
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
             
    if mode == "D":
        fr_to_netcdf_tools.IconD2toNetCDF(x, "//vs-grp08.zih.tu-dresden.de/hwstore/IconD2/tot_prec/",
                        lon, lat,
                        "data_generation/Sachsen_nearestpoints.npz",
                        "//vs-grp07.zih.tu-dresden.de/howa/work/students/Mohamed_Elghorab/netCDFs/icond2_ev.nc")
    elif mode == "eps":
        fr_to_netcdf_tools.IconD2EPStoNetCDF(x, "//vs-grp08.zih.tu-dresden.de/hwstore/IconD2eps/tot_prec/",
                          lon, lat,
                          "data_generation/Sachsen_nearestpoints.npz",
                          "//vs-grp07.zih.tu-dresden.de/howa/work/students/Mohamed_Elghorab/netCDFs/icond2eps_ev11.nc")
  

def create_Cosmo(file, mode= None):
    
    lib = event_file(file)
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
             
    if mode == "D":
        fr_to_netcdf_tools.CosmoD2toNetCDF(x, "//vs-grp08.zih.tu-dresden.de/hwstore/CosmoD2/",
                        lon, lat,
                        "Sachsen_nearestpoints_cosmo.npz",
                        "//vs-grp07.zih.tu-dresden.de/howa/work/students/Mohamed_Elghorab/netCDFs/fertig/cosmod2.nc")
    elif mode == "eps":
        fr_to_netcdf_tools.CosmoD2EPStoNetCDF(x, "//vs-grp08.zih.tu-dresden.de/hwstore/CosmoD2eps/",
                          lon, lat,
                          "Sachsen_nearestpoints_cosmo.npz",
                          "//vs-grp07.zih.tu-dresden.de/howa/work/students/Mohamed_Elghorab/netCDFs/fertig/cosmod2eps.nc")    
    
    
    
    ###############################################################################################


################################################################################################
# getting all datapoint indices covering the extents of Radolan-RW raster for
# Sachsen after WGS84 > lat 50.1 - 51.8 | long 11.7 - 15.2 <

lon, lat, id_lon, id_lat = fr_to_netcdf_tools.target(50.1,51.8,11.7,15.2,"radolanrx",version=4)

################################################################################################

tic = time.time()
create_radar('data_generation/cosmo_events_erad.txt')
print((time.time()-tic)/3600) # operation time in hours


# tic = time.time()
# create_icond2('data_generation/Icon_events_strt.txt',"D")
# print((time.time()-tic)/3600) # operation time in hours


# tic = time.time()
# create_icond2('data_generation/Icon_events_strt.txt',"eps")
# print((time.time()-tic)/3600) # operation time in hours


# tic = time.time()
# create_Cosmo('data_generation/cosmo_events_d.txt',"D")
# print((time.time()-tic)/3600) # operation time in hours


# tic = time.time()
# create_Cosmo('data_generation/cosmo_events_e.txt',"eps")
# print((time.time()-tic)/3600) # operation time in hours






