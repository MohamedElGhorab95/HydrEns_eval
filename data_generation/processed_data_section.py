# -*- coding: utf-8 -*-
"""
Created on Sun Jul 16 15:53:28 2023

@author: M Elghorab
"""


from engine.fr_entities_tools import *
import matplotlib.pyplot as plt
import seaborn as sns


def load_event_file(leadtime):
    source_file_path = "data_generation/nc_paths.txt"
    with open(source_file_path, 'r') as f:
        # Read all the lines in the file
        lines = f.readlines()
        
        netcdf = lines[1] 
        
        
     
        radar_file = netcdf+ 'fertig' + '/' + 'radRW_ICO.nc'
        deter_file = netcdf+ str(leadtime) + '/' + '{}hour_icond2.nc'.format(leadtime)
        ensem_file = netcdf+ str(leadtime) + '/' + '{}hour_icond2eps.nc'.format(leadtime)
        
    files = [radar_file, deter_file, ensem_file]
    return files




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


for cat in ['Oelsnitz', 'Dohna', 'Zittau']:
    
    for lt in [3,12]:
        
        lib = load_event_file(lt)
        
        det = Deterministic_run(lib[1]).gen_deterministic_field().extract_by_shp(load_catchment(cat)).avg_areal_prec()
        
        eps = Ensemble_run(lib[2]).eps_extract_by_shp(load_catchment(cat)).avg_areal_prec()
        
        sns.set_theme(style="whitegrid")
        fig, ax = plt.subplots(dpi=750)

       
        det.average.plot(label='Main run')
        eps.average['icond2eps_0'].plot(label='ICOND2EPS_m1')
        eps.average['icond2eps_9'].plot(label='ICOND2EPS_m10')
        eps.average['icond2eps_19'].plot(label='ICOND2EPS_m20')

        ax.legend()
        ax.set_title("MAP | {} \nLeadtime = {}hrs".format(cat,lt))
        ax.set_ylabel('Rainfall (mm/3hrs)')
        plt.show()
        plt.close()