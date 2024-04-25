# -*- coding: utf-8 -*-
"""
Created on Tue Aug  1 18:47:26 2023

@author: M Elghorab
"""

from engine.fr_entities_tools import R_Forecast
import xarray as xr
import os
import numpy as np
from tqdm import tqdm

def find_folders(directory, date_list):
    folder_list = []
    for root, dirs, files in os.walk(directory):
        for folder in dirs:
            for date in date_list:
                if date in folder:
                    folder_path = os.path.join(root, folder)
                    if len(os.listdir(folder_path)) > 0:
                        folder_list.append(folder)
                        # No need to check other dates once a match is found
                        break
    return folder_list


def gen_runoff_datasets(dates, catchment):
    
    path = "//vs-grp07.zih.tu-dresden.de/howa/work/students/Mohamed_Elghorab/Mohamed Elghorab_MSc_Working_files/forecast_system/ForData"
    
    for lead in np.arange(3,25,3):
        print('generating for leadtime: {}hrs'.format(lead))
        print('finding folders...........')
        fs = find_folders(path, dates)
        # Create a progress bar
        progress_bar = tqdm(total=len(fs), desc='Generating arrays.................')
        for f in fs:
                      
            try:
                run = R_Forecast(path+'/'+f+'/'+f+catchment+'.nc')
              
            
                q10 = run.gen_quantiles(10).fr.isel(time=lead*4) ; q25 = run.gen_quantiles(25).fr.isel(time=lead*4)
                q50 = run.gen_quantiles(50).fr.isel(time=lead*4) ; q75 = run.gen_quantiles(75).fr.isel(time=lead*4)
                q90 = run.gen_quantiles(90).fr.isel(time=lead*4) ; mean = run.gen_quantiles('mean').fr.isel(time=lead*4)
                eps = run.gen_ensmembers().fr.isel(time=lead*4)
                if f == fs[0]:
                    Q_q10 = q10 ; Q_q25 = q25 ; Q_q50 = q50
                    Q_q75 = q75 ; Q_q90 = q90 ; Qm = mean
                    EPS = eps
                    
                    
                else:
                    Q_q10 = xr.concat([Q_q10,q10],dim='time') ; Q_q25 = xr.concat([Q_q25,q25],dim='time')
                    Q_q50 = xr.concat([Q_q50,q50],dim='time') ; Q_q75 = xr.concat([Q_q75,q75],dim='time')
                    Q_q90 = xr.concat([Q_q90,q90],dim='time') ; Qm = xr.concat([Qm,mean],dim='time')
                    EPS = xr.concat([EPS,eps], dim='time')
            except:
                continue 
            
            progress_bar.update(1)
        
        # Close the progress bar
        progress_bar.close()
        
        dataset_final = xr.Dataset({'Q_q10':Q_q10.sortby('time').drop_duplicates('time')
                                    ,'Q_q25':Q_q25.sortby('time').drop_duplicates('time')
                                    , 'Q_q50':Q_q50.sortby('time').drop_duplicates('time')
                                    , 'Q_q75':Q_q75.sortby('time').drop_duplicates('time')
                                    , 'Q_q90':Q_q90.sortby('time').drop_duplicates('time')
                                    , 'Q_mean':Qm.sortby('time').drop_duplicates('time')})
        
        print('exporting dataset.........')
        # dataset_final.to_netcdf('//vs-grp07.zih.tu-dresden.de/howa/work/students/Mohamed_Elghorab/Mohamed Elghorab_MSc_Working_files/netCDFs/Runoff/{}hrs_leadtime_{}.nc'.format(lead,catchment))
        
        dataset_final.to_netcdf('E:/Data2/NetCDFs/Runoff/{}hrs_leadtime_{}.nc'.format(lead,catchment))
        
        
        EPS = EPS.sortby('time')
        
        ENS = EPS.drop_duplicates('time')
        
        # ENS.to_netcdf('//vs-grp07.zih.tu-dresden.de/howa/work/students/Mohamed_Elghorab/Mohamed Elghorab_MSc_Working_files/netCDFs/Runoff/{}hrs_leadtime_{}_EPS.nc'.format(lead,catchment))
        ENS.to_netcdf('E:/Data2/NetCDFs/Runoff/{}hrs_leadtime_{}_EPS.nc'.format(lead,catchment))




# cat = ['_5661371_Oelsnitz_WeisseElster_data', '_56611313_BadElster_WeisseElster_data', '_5661311_Adorf_WeisseElster_data']
# cat = ['_53718979_Dohna_Mueglitz_data', '_5371831_Lauenstein4_Mueglitz_data', '_5371823_Geising1_WeisseMueglitz_data' ]
cat = ['_67414799_Zittau_Mandau_data', '_67414651_Niederoderwitz_Mandau_data', '_67414511_Grossschoenau_Mandau_data', '_67414311_Seifhennersdorf_Mandau_data']

for c in cat:
    gen_runoff_datasets(['2021','2022','2023'], c)
    

# cat = ['_5661371_Oelsnitz_WeisseElster_data', '_56611313_BadElster_WeisseElster_data', '_5661311_Adorf_WeisseElster_data']    
# for c in cat:
#     gen_runoff_datasets(['2021','2022','2023'], c)
  
    
# cat = ['_53718979_Dohna_Mueglitz_data', '_5371831_Lauenstein4_Mueglitz_data', '_5371823_Geising1_WeisseMueglitz_data' ]    
# for c in cat:
#     gen_runoff_datasets(['2021','2022','2023'], c)    