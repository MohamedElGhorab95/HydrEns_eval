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
    
    
    
        
    data1   = {key: [None,region1[key]]  for key in region1.keys()}
    
    for a in data1.keys():
        res = []
        if a == 'Oelsnitz':
            obs = R_Observation('C:/Users/User/Downloads/Oelsnitz.csv')
        elif a == 'Adorf':
            obs = R_Observation('C:/Users/User/Downloads/Adorf.csv')
        elif a == 'Bad Elster':
            obs = R_Observation('C:/Users/User/Downloads/Bad_Elster.csv')
            
        
        for lead in range(3,25,3):
            if a == 'Oelsnitz':
                full = xr.open_dataset('//vs-grp07.zih.tu-dresden.de/howa/work/students/Mohamed_Elghorab/netCDFs/Runoff/{}hrs_leadtime__5661371_Oelsnitz_WeisseElster_data_EPS.nc'.format(lead))
                res.append(calc_metric(obs, full, region1[a][1]*0.25, metric))
            elif a == 'Adorf':
                full = xr.open_dataset('//vs-grp07.zih.tu-dresden.de/howa/work/students/Mohamed_Elghorab/netCDFs/Runoff/{}hrs_leadtime__5661311_Adorf_WeisseElster_data_EPS.nc'.format(lead))
                res.append(calc_metric(obs, full, region1[a][1]*0.25, metric))
            elif a == 'Bad Elster':
                full = xr.open_dataset('//vs-grp07.zih.tu-dresden.de/howa/work/students/Mohamed_Elghorab/netCDFs/Runoff/{}hrs_leadtime__56611313_BadElster_WeisseElster_data_EPS.nc'.format(lead))
                res.append(calc_metric(obs, full, region1[a][1]*0.25, metric))
        
        data1[a][0]=res   
    return [data1]
    


def plot_spatial(list_of_dics, metric):
    
    max_lead = 24
    leads = list(range(3,max_lead+1,3))
    
    region = 0
    regions = ['Weiße Elster region']
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
        # plt.savefig("//vs-grp07.zih.tu-dresden.de/howa/work/students/Mohamed_Elghorab/results/spatial extent analysis/{}/{}".format(metric, regions[region]), bbox_inches = 'tight')        
        region +=1
        

    return plt.show()



mets = [ 'CRPS', "RMSE", 'Area under ROC curve'  ]
for m in mets:

    
    df = create_for_regions(m)
  
    plot_spatial(df,m)


    

