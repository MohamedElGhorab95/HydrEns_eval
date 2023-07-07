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
    # splitting to smaller sub events for faster calculation and memory usage
    # split_events = {}
    # for key, (start, end, duration) in evs.items():
    #     num_splits = int(np.ceil(duration / 2))  # Calculate the number of splits
    #     split_duration = int(np.ceil(duration / num_splits))  # Calculate the duration for each split
    
    #     for i in range(num_splits):
    #         split_start = start + np.timedelta64(i * 2, 'D')  # Calculate the split start time
    #         split_end = split_start + np.timedelta64(split_duration, 'D')  # Calculate the split end time
        
    #         split_events[key * 100 + i] = (split_start, split_end, split_duration)
    # renumbered_events = {i: split_events[key] for i, key in enumerate(split_events, start=1)}
    # return renumbered_events
    return evs



def process_lead(lead, data_low_level, temp_res, run_name, events):
    # Your existing code inside the for loop for leadt
    # ...
    for ev in events.keys():
        print('leadtime {}hrs \ngenerating event............ {}'.format(lead,ev))
        
        t1 = data_low_level.sel(start_time=slice(events[ev][0],events[ev][1]))
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
            
        if ev == 1:
            dataset_final = lead_ds
        else:
            dataset_final = xr.concat([dataset_final,lead_ds], dim='time')
    
    # Save the lead dataset
    dataset_final.to_netcdf('C:/netCDFs/{}/{}hour_{}.nc'.format(lead, lead, run_name))


def gen_leadtime_netcdf(events, netcdf_file, run_name):
    # Your existing code ...
    data_low_level = xr.open_dataset(netcdf_file,chunks='auto').isel(forecast_time=slice(0,100))
    data_low_level = data_low_level.rename({"time": "start_time", "forecast_time": "time"})
    
    leadt= np.arange(3,25,3)
    temp_res = 3
    
    # Create a ThreadPoolExecutor with the desired number of threads
    num_threads = 4
    executor = concurrent.futures.ThreadPoolExecutor(max_workers=num_threads)
    
    # Submit the tasks to the executor
    futures = [executor.submit(process_lead, lead, data_low_level, temp_res, run_name, events) for lead in leadt]
    
    # Wait for all tasks to complete
    concurrent.futures.wait(futures)













events = load_events()



# generating lead time netcdfs for icond2
gen_leadtime_netcdf(events, '//vs-grp07.zih.tu-dresden.de/howa/work/students/Mohamed_Elghorab/netCDFs/fertig/icond2eps.nc','icond2eps')































