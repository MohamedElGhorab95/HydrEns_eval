# -*- coding: utf-8 -*-
"""
Created on Fri Jul  7 10:28:06 2023

@author: M Elghorab
"""

import xarray as xr
import xskillscore as xs
import numpy as np
import seaborn as sns
from HydrEns_eval.engine.fr_entities_tools import *
from HydrEns_eval.engine.fr_Conttools import *
from HydrEns_eval.engine.fr_ROCtools import *


rad = Observation("C:/netCDFs/fertig/radRW_ICO.nc").gen_observation_field().extract_by_shp("shp/Mugliz/mugliz_cats.shp").aggr_temporal(3).avg_areal_prec()
# icond2 = Deterministic_run("C:/netCDFs/3/3hour_icond2.nc").gen_deterministic_field().extract_by_shp("shp/Mugliz/mugliz_cats.shp").avg_areal_prec()
icond2eps = Ensemble_run("C:/netCDFs/18/18hour_icond2eps.nc").eps_extract_by_shp("shp/Mugliz/mugliz_cats.shp").avg_areal_prec()

members = list(icond2eps.average.variables)[3:]

mem2 = icond2eps.average[members[1]]


pod = []
pofd = []
auc = []

for mem in members:
    r = ROC(rad,icond2eps.average[mem])
    auc.append(float(r.roc_auc()))
    pod.append(r.pod.values)
    pofd.append(r.pofd.values)
    # r.plot_roc(mem)


np.quantile(auc, 0.9)
np.mean(auc)
pod90 = [np.quantile([arr[i] for arr in pod], q=0.1) for i in range(len(pod[1]))]
pofd90 = [np.quantile([arr[i] for arr in pofd], q=0.1) for i in range(len(pofd[1]))]

def plot_roc(lead):
    '''
    Returns
    -------
    a single ROC plot for the entire extents of the data array.

    '''
    
   
    sns.set_theme(style="whitegrid")
    
    fig, ax = plt.subplots(dpi=750)
    
    ax.plot(pofd90,pod90,color='r')
    # ax.scatter(self.pofd,self.pod,color='r', marker='*')
    
    # Add a 45-degree line
    ax.plot([0, 1], [0, 1], transform=ax.transAxes, ls='--', color='k')
        
    ax.set_xlabel('PROBABILITY OF FALSE DETECTION')
    ax.set_ylabel('PROBABILITY OF DETECTION')
    
    ax.set_xlim(0,1)
    ax.set_ylim(0,1)
    
    plt.title("ROC at lead time: {}hrs\nagainst ICOND2EPS 90% Quantile (of members)".format(lead))

    # Show the plot
    plt.show()
    
plot_roc(18)
