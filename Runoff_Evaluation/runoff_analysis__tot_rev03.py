# -*- coding: utf-8 -*-
"""
Created on Fri Aug  4 19:19:49 2023

@author: M Elghorab
"""

import xarray as xr
from engine.fr_ROCtools import *
from engine.fr_Conttools import *
from engine.fr_Fullens import *
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator
import seaborn as sns





def calc_metric(observation, forecast , threshold ,metric):
    print("calculating performance metric..........",flush= True)
    
    
    if metric == 'Area under ROC curve':
        return ROC_totens(observation, forecast, threshold).area
    elif metric == 'CRPS':
        return Full_Ens(observation.fr, forecast).crps()
    elif metric == "RMSE":
        return Full_Ens(observation.fr, forecast).rmse()
    

def create_for_regions( metric):
    
    region1 = {"Oelsnitz":[328,25.6], 'Adorf':[170,13.3] ,'Bad Elster':[48,3.72]}
    
    region2 = {'Dohna':[200,30.4], 'Lauenstein':[76,12.3], 'Geising':[26,6.6] }
    
    # region3 = { 'Zittau':[279,63.5], 'Grossschoenau': [162,41.6],'Seifhennersdorf':[75,20.2] ,'Niederoderwitz':[29,12.6]}
    # region3 = { 'Zittau':[279,3.35], 'Grossschoenau': [162,2.34],'Seifhennersdorf':[75,0.928] ,'Niederoderwitz':[29,0.282]}
    region3 = { 'Zittau':[279,3.35], 'Grossschoenau': [162,2.34] ,'Niederoderwitz':[29,0.282]}
   
    
        
    data1   = {key: [None,region1[key]]  for key in region1.keys()}
    data2   = {key: [None,region2[key]]  for key in region2.keys()}
    data3   = {key: [None,region3[key]]  for key in region3.keys()}
    
    
    for a in data1.keys():
        res = []
        
        for lead in range(3,25,3):
            if a == 'Oelsnitz':
                obs = R_Observation('E:/HydrEns_eval/Runoff_Evaluation/Oelsnitz.csv')
                # full = xr.open_dataset('//vs-grp07.zih.tu-dresden.de/howa/work/students/Mohamed_Elghorab/Mohamed Elghorab_MSc_Working_files/netCDFs/Runoff/{}hrs_leadtime__5661371_Oelsnitz_WeisseElster_data_EPS.nc'.format(lead))
                
                full = xr.open_dataset('//vs-grp07.zih.tu-dresden.de/howa/work/students/Mohamed_Elghorab/Mohamed Elghorab_MSc_Working_files/netCDFs/Runoff/{}hrs_leadtime__5661371_Oelsnitz_WeisseElster_data_EPS.nc'.format(lead))
                
                res.append(calc_metric(obs, full, region1[a][1]*0.25, metric))
            elif a == 'Adorf':
                obs = R_Observation('E:/HydrEns_eval/Runoff_Evaluation/Adorf.csv')
                # full = xr.open_dataset('//vs-grp07.zih.tu-dresden.de/howa/work/students/Mohamed_Elghorab/Mohamed Elghorab_MSc_Working_files/netCDFs/Runoff/{}hrs_leadtime__5661311_Adorf_WeisseElster_data_EPS.nc'.format(lead))
                
                full = xr.open_dataset('//vs-grp07.zih.tu-dresden.de/howa/work/students/Mohamed_Elghorab/Mohamed Elghorab_MSc_Working_files/netCDFs/Runoff/{}hrs_leadtime__5661311_Adorf_WeisseElster_data_EPS.nc'.format(lead))
                
                res.append(calc_metric(obs, full, region1[a][1]*0.25, metric))
            elif a == 'Bad Elster':
                obs = R_Observation('E:/HydrEns_eval/Runoff_Evaluation/Bad_Elster.csv')
                # full = xr.open_dataset('//vs-grp07.zih.tu-dresden.de/howa/work/students/Mohamed_Elghorab/Mohamed Elghorab_MSc_Working_files/netCDFs/Runoff/{}hrs_leadtime__56611313_BadElster_WeisseElster_data_EPS.nc'.format(lead))
                
                full = xr.open_dataset('//vs-grp07.zih.tu-dresden.de/howa/work/students/Mohamed_Elghorab/Mohamed Elghorab_MSc_Working_files/netCDFs/Runoff/{}hrs_leadtime__56611313_BadElster_WeisseElster_data_EPS.nc'.format(lead))
                
                res.append(calc_metric(obs, full, region1[a][1]*0.25, metric))
        
        data1[a][0]=res   
        
        
    for a in data2.keys():
        res = []
        
        for lead in range(3,25,3):
            if a == 'Dohna':
                obs = R_Observation('E:/HydrEns_eval/Runoff_Evaluation/Dohna.csv')
                # full = xr.open_dataset('//vs-grp07.zih.tu-dresden.de/howa/work/students/Mohamed_Elghorab/Mohamed Elghorab_MSc_Working_files/netCDFs/Runoff/{}hrs_leadtime__53718979_Dohna_Mueglitz_data_EPS.nc'.format(lead))
               
                full = xr.open_dataset('//vs-grp07.zih.tu-dresden.de/howa/work/students/Mohamed_Elghorab/Mohamed Elghorab_MSc_Working_files/netCDFs/Runoff/{}hrs_leadtime__53718979_Dohna_Mueglitz_data_EPS.nc'.format(lead))
                
                res.append(calc_metric(obs, full, region2[a][1]*0.25, metric))
            elif a == 'Lauenstein':
                obs = R_Observation('E:/HydrEns_eval/Runoff_Evaluation/Lauenstein.csv')
                # full = xr.open_dataset('//vs-grp07.zih.tu-dresden.de/howa/work/students/Mohamed_Elghorab/Mohamed Elghorab_MSc_Working_files/netCDFs/Runoff/{}hrs_leadtime__5371831_Lauenstein4_Mueglitz_data_EPS.nc'.format(lead))
                
                full = xr.open_dataset('//vs-grp07.zih.tu-dresden.de/howa/work/students/Mohamed_Elghorab/Mohamed Elghorab_MSc_Working_files/netCDFs/Runoff/{}hrs_leadtime__5371831_Lauenstein4_Mueglitz_data_EPS.nc'.format(lead))
                
                res.append(calc_metric(obs, full, region2[a][1]*0.25, metric))
            elif a == 'Geising':
                obs = R_Observation('E:/HydrEns_eval/Runoff_Evaluation/Geising.csv')
                # full = xr.open_dataset('//vs-grp07.zih.tu-dresden.de/howa/work/students/Mohamed_Elghorab/Mohamed Elghorab_MSc_Working_files/netCDFs/Runoff/{}hrs_leadtime__5371823_Geising1_WeisseMueglitz_data_EPS.nc'.format(lead))
                
                full = xr.open_dataset('//vs-grp07.zih.tu-dresden.de/howa/work/students/Mohamed_Elghorab/Mohamed Elghorab_MSc_Working_files/netCDFs/Runoff/{}hrs_leadtime__5371823_Geising1_WeisseMueglitz_data_EPS.nc'.format(lead))
                
                res.append(calc_metric(obs, full, region2[a][1]*0.25, metric))
        
        data2[a][0]=res  
        
    for a in data3.keys():
        res = []
        
        for lead in range(3,25,3):
            if a == 'Zittau':
                obs = R_Observation('E:/HydrEns_eval/Runoff_Evaluation/Zittau.csv')
                # full = xr.open_dataset('//vs-grp07.zih.tu-dresden.de/howa/work/students/Mohamed_Elghorab/Mohamed Elghorab_MSc_Working_files/netCDFs/Runoff/{}hrs_leadtime__67414799_Zittau_Mandau_data_EPS.nc'.format(lead))
                full = xr.open_dataset('//vs-grp07.zih.tu-dresden.de/howa/work/students/Mohamed_Elghorab/Mohamed Elghorab_MSc_Working_files/netCDFs/Runoff/{}hrs_leadtime__67414799_Zittau_Mandau_data_EPS.nc'.format(lead))
                res.append(calc_metric(obs, full, region3[a][1], metric))
            elif a == 'Grossschoenau':
                obs = R_Observation('E:/HydrEns_eval/Runoff_Evaluation/Grossschoenau.csv')
                # full = xr.open_dataset('//vs-grp07.zih.tu-dresden.de/howa/work/students/Mohamed_Elghorab/Mohamed Elghorab_MSc_Working_files/netCDFs/Runoff/{}hrs_leadtime__67414511_Grossschoenau_Mandau_data_EPS.nc'.format(lead))
                full = xr.open_dataset('//vs-grp07.zih.tu-dresden.de/howa/work/students/Mohamed_Elghorab/Mohamed Elghorab_MSc_Working_files/netCDFs/Runoff/{}hrs_leadtime__67414511_Grossschoenau_Mandau_data_EPS.nc'.format(lead))
                res.append(calc_metric(obs, full, region3[a][1], metric))
            # elif a == 'Seifhennersdorf':
            #     obs = R_Observation('E:/HydrEns_eval/Runoff_Evaluation/Seifhennersdorf.csv')
            #     # full = xr.open_dataset('//vs-grp07.zih.tu-dresden.de/howa/work/students/Mohamed_Elghorab/Mohamed Elghorab_MSc_Working_files/netCDFs/Runoff/{}hrs_leadtime__67414311_Seifhennersdorf_Mandau_data_EPS.nc'.format(lead))
            #     full = xr.open_dataset('//vs-grp07.zih.tu-dresden.de/howa/work/students/Mohamed_Elghorab/Mohamed Elghorab_MSc_Working_files/netCDFs/Runoff/{}hrs_leadtime__67414311_Seifhennersdorf_Mandau_data_EPS.nc'.format(lead))
            #     res.append(calc_metric(obs, full, region3[a][1], metric))
            elif a == 'Niederoderwitz':
                obs = R_Observation('E:/HydrEns_eval/Runoff_Evaluation/Niederoderwitz.csv')
                # full = xr.open_dataset('//vs-grp07.zih.tu-dresden.de/howa/work/students/Mohamed_Elghorab/Mohamed Elghorab_MSc_Working_files/netCDFs/Runoff/{}hrs_leadtime__67414651_Niederoderwitz_Mandau_data_EPS.nc'.format(lead))
                full = xr.open_dataset('//vs-grp07.zih.tu-dresden.de/howa/work/students/Mohamed_Elghorab/Mohamed Elghorab_MSc_Working_files/netCDFs/Runoff/{}hrs_leadtime__67414651_Niederoderwitz_Mandau_data_EPS.nc'.format(lead))
                res.append(calc_metric(obs, full, region3[a][1], metric))
        
        data3[a][0]=res   
    
    return [data1, data2, data3]
    


def plot_spatial(list_of_dics, metric):
    
    max_lead = 24
    leads = list(range(3,max_lead+1,3))
    
    region = 0
    regions = ['Weiße Elster region', 'Müglitz region', 'Mandau region']
    for l in list_of_dics:
        sns.set_theme(style="whitegrid")
        fig, ax = plt.subplots(dpi=1000)
        
        l = dict(sorted(l.items()))
        
        for u in l.keys():
            
           
            plt.plot(leads, l[u][0],label='{} |{}km\u00b2'.format(u,l[u][1][0]) )
        
        # TODO change metric name in y axis and folder    
        ax.legend()
        plt.xticks(leads)
        ax.set_title("Performance of runoff forecast \n{}".format(regions[region]))
        ax.set_xlabel('Leadtime (hrs)')
        ax.set_ylabel('{}\n Full Ensemble'.format(metric))
        plt.xlim(3, max_lead)
        plt.ylim(0,6)
        # plt.savefig("//vs-grp07.zih.tu-dresden.de/howa/work/students/Mohamed_Elghorab/Mohamed Elghorab_MSc_Working_files/results/spatial extent analysis/{}/{}".format(metric, regions[region]), bbox_inches = 'tight')        
        region +=1
        

    return plt.show()



# mets = [ 'CRPS', "RMSE", 'Area under ROC curve'  ]
mets = ['RMSE'  ]
for m in mets:

    
    df = create_for_regions(m)
  
    plot_spatial(df,m)



plot_spatial(df,'Normalized RMSE')
