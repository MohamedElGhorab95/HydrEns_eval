# -*- coding: utf-8 -*-
"""
Created on Tue May 30 17:48:07 2023

@author: M Elghorab
"""
# importing the required packages

import numpy as np
import math
import concurrent.futures
import xarray as xr
import datetime as dt



        

def load_events():
    # reading the events record
    with open('data_generation/Icon_events_strt.txt','r') as f:
        lines = f.readlines()[1:]
    ################################################################################################
    # creating a dictionary of the events
    events = [[],[],[]]
    evid = [i for i in range(1,21)]  
    for l in lines:   
        # events[0].append((int(l[9:13]),int(l[6:8]),int(l[3:5]),0))
        date_st = l[3:13].split('.')
        date_en = l[14:24].split('.')
        date_st = np.datetime64(f'{date_st[2]}-{date_st[1]}-{date_st[0]}T{0:02d}:00:00')
        date_en = np.datetime64(f'{date_en[2]}-{date_en[1]}-{date_en[0]}T{0:02d}:00:00')
        events[0].append(date_st)
        events[1].append(date_en)
        # adding event time span to the dictionary
        timedelta = date_en - date_st
        timedelta = timedelta / np.timedelta64(1, 'D')
        events[2].append(int(math.ceil((timedelta).astype(float)))+1)
    # event library dictionary
    evs =   dict(zip(evid,zip(events[0],events[1],events[2])))
    
    return evs



def process_lead(lead, netcdf_file, temp_res, run_name, events):
    # Your existing code inside the for loop for leadt
    # ...
    
    
    for ev in events.keys():
        data_low_level = xr.open_dataset(netcdf_file,chunks='auto').isel(forecast_time=slice(0,100))
        data_low_level = data_low_level.rename({"time": "start_time", "forecast_time": "time"})
        print('leadtime {}hrs \ngenerating event............ {}'.format(lead,ev))
        
        t1 = data_low_level.sel(start_time=slice(events[ev][0],events[ev][1]))
        l_lat = np.linspace(51.921707, 49.991165, 212).tolist()
        l_lon = np.linspace(11.634748, 15.226268, 253).tolist()
        # assigning the actual coordinates
        t1 = t1.assign_coords(lat=l_lat, lon=l_lon)
        # drop unnecessary variables
        t1 = t1.drop_vars(["longitude", "latitude"])
        del data_low_level
        t1 = t1.isel(time = slice(1, len(t1.time))) 
        t1 = t1.resample(time="1H").sum()  # mm/15min -> mm/hr
        t1['time'] = t1.time + np.timedelta64(dt.timedelta(minutes=45)) 
        day0 = t1.isel(start_time=7)
        
        for cyc in range(1,len(t1.start_time)-7):
            t2 = t1.isel(start_time= cyc +7).resample(time="{}H".format(temp_res),offset=dt.timedelta(hours=temp_res-1)).sum()    
            t2['time'] = t2.start_time + t2.time
            
            
            if cyc == 1:
                if lead >3:
                    val = day0.resample(time="{}H".format(temp_res),offset=dt.timedelta(hours=temp_res-1)).sum()
                    val['time'] = val.start_time + val.time
                    lead_ds = val.isel(time=slice(1,int(lead/3)))
                    lead_ds = xr.concat([lead_ds,t2.isel(time=int(lead/3)-1)], dim='time')
                else:
                    lead_ds = t2.isel(time=int(lead/3)-1)
            else:
                lead_ds = xr.concat([lead_ds,t2.isel(time=int(lead/3)-1)], dim='time')   
        del t1, t2                        
        if ev == 1:
            dataset_final = lead_ds
        else:
            dataset_final = xr.concat([dataset_final,lead_ds], dim='time')
    
    # Save the lead dataset
    dataset_final.to_netcdf('//vs-grp07.zih.tu-dresden.de/howa/work/students/Mohamed_Elghorab/netCDFs/{}/{}hour_{}.nc'.format(lead, lead, run_name))


def gen_leadtime_netcdf(events, netcdf_file, run_name):
    # Your existing code ...
       
    
    leadt= np.arange(3,25,3)
    temp_res = 3
    
    # Create a ThreadPoolExecutor with the desired number of threads
    num_threads = 2
    executor = concurrent.futures.ProcessPoolExecutor(max_workers=num_threads)
    
    # Submit the tasks to the executor
    futures = [executor.submit(process_lead, lead, netcdf_file, temp_res, run_name, events) for lead in leadt]
    
    # Wait for all tasks to complete
    concurrent.futures.wait(futures)













events = load_events()

if __name__ == '__main__':

    # generating lead time netcdfs for icond2
    gen_leadtime_netcdf(events, '//vs-grp07.zih.tu-dresden.de/howa/work/students/Mohamed_Elghorab/netCDFs/fertig/icond2eps.nc','icond2eps')































