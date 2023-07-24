# -*- coding: utf-8 -*-
"""
Created on Sat Jul 15 18:44:13 2023

@author: M Elghorab
"""

from engine.fr_entities_tools import *
from engine.fr_Fullens import *
from engine.fr_Conttools import *
from engine.fr_ROCtools import *
import matplotlib.pyplot as plt
import seaborn as sns

def load_event_file(leadtime):
    source_file_path = "data_generation/nc_paths.txt"
    with open(source_file_path, 'r') as f:
        # Read all the lines in the file
        lines = f.readlines()
        # TODO change  to 0 to read from server
        netcdf = lines[1] 
        
        
     
        radar_file = netcdf+ 'fertig' + '/' + 'radRW_ICO.nc'
        
        ensem_file = netcdf+ str(leadtime) + '/' + '{}hour_icond2eps.nc'.format(leadtime)
        return [radar_file, ensem_file]
    
   


def load_catchment(name):
    # loading the catchment of interest
    source_file_path = "data_generation/paths.txt"
    with open(source_file_path, 'r') as f:
        # Read all the lines in the file
        lines = f.readlines()
        if name == 'Zittau':
            return lines[2].strip()
        if name == 'Dohna':
            return lines[4].strip()
        if name == 'Oelsnitz':
            return lines[6].strip()
        if name == 'Adorf':
            return lines[8].strip()
        if name == 'Bad Elster':
            return lines[10].strip()
        if name == 'Lauenstein':
            return lines[12].strip()
        if name == 'Niederoderwitz':
            return lines[14].strip()
        if name == 'Seifhennersdorf':
            return lines[16].strip()
        if name == 'Geising':
            return lines[18].strip()
        if name == 'Grossschoenau':
            return lines[20].strip()


lib = load_event_file(3)

# creating instances 
rad = Observation(lib[0]).gen_observation_field().aggr_temporal(3)

fore = Ensemble_run(lib[1])

# fore.plot([2021, 8, 5, 15],20)

rad = rad.extract_by_shp(load_catchment('Zittau'))

fore = fore.eps_extract_by_shp(load_catchment('Zittau'))

# rad.plot([2021,8,5,15])

rad = rad.avg_areal_prec()

fore_ar = fore.avg_areal_prec()


# sns.set_theme(style="whitegrid")
# fig, ax = plt.subplots(dpi=750)

# thr = [0,2.5,5,7.5,10,12.5,15,17.5,20, 22.5]

# rad.average.plot(label='Radar')
# fore_ar.average['icond2eps_19'].plot(label='ICOND2EPS_m20')

# # Add horizontal lines for thresholds
# for t in thr:
#     ax.axhline(y=t, color='red', linestyle='--', linewidth=1)


# ax.legend()
# ax.set_title("MAP | Zittau catchment \nLeadtime = {}hrs".format(3))
# ax.set_ylabel('Rainfall (mm/3hrs)')
# plt.show()

# print(Full_Ens(rad, fore_ar).crps())

# print(Full_Ens(rad, fore_ar).rmse())


# print(CONT(rad, fore_ar, 2.5).acc(10))

qs = [10,25,50,75,90, 'mean']

for q in qs:
    r=ROC(rad, fore_ar)
    print(r.roc_auc(quantile=q))
    r.plot_roc(3, q)

