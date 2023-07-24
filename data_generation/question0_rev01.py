# -*- coding: utf-8 -*-
"""
Created on Fri Jul 14 17:17:36 2023

@author: M Elghorab
"""

from engine.fr_entities_tools import *
from engine.fr_Fullens import *
import matplotlib.pyplot as plt
import seaborn as sns
import time


def load_event_file(leadtime):
    source_file_path = "data_generation/nc_paths.txt"
    with open(source_file_path, 'r') as f:
        # Read all the lines in the file
        lines = f.readlines()
        
        netcdf = lines[1] 
        
        
     
        radar_file = netcdf+ 'fertig' + '/' + 'radRW_ICO.nc'
        # deter_file = netcdf+ str(leadtime) + '/' + '{}hour_icond2.nc'.format(leadtime)
        ensem_file = netcdf+ str(leadtime) + '/' + '{}hour_icond2eps.nc'.format(leadtime)
        
    # files = [radar_file, deter_file, ensem_file]
    files = [radar_file,  ensem_file]
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



def gen_for_catchment(catchment, locations):
    
    print('generating instances..........',flush= True)
    # print(f"generating instances..........")
    # progress bar graphics
   
    
    # creating instances 
    rad = Observation(locations[0]).gen_observation_field().aggr_temporal(3)
    eps = Ensemble_run(locations[1])
   
    print('cropping to: {}'.format(catchment),flush= True)
    # print(f"cropping to: {catchment}")
    # cropping for selected catchment
    rad_c = rad.extract_by_shp(load_catchment(catchment))

    eps_c = eps.eps_extract_by_shp(load_catchment(catchment))
    # garbage collection
    del rad, eps
   
    # average areal precipitation
    print('calculaing areal averages',flush= True)
    # print(f"calculaing areal averages")
    rad_c = rad_c.avg_areal_prec()
    eps_c = eps_c.avg_areal_prec()
   
    library = {'radar': rad_c,
               'ensemble_run': eps_c
               }
   
    del rad_c, eps_c
    
    return library




def evaluate_for_catchment(cats):

   
    # TODO change
    max_lead = 21

    crps = []
    rmse = []
   
    for lead in range(3,max_lead+1,3):
        print('Catchment: {}\ncalculating for lead time: {}hrs'.format(cats,lead),flush= True)
        
        # loading the files
        files = load_event_file(lead)
        # creating the instances for a specific catchment
        library = gen_for_catchment(cats,files)

        crps.append(Full_Ens(library['radar'], library['ensemble_run']).crps())
        rmse.append(Full_Ens(library['radar'], library['ensemble_run']).rmse())
            
        del library
    
    results = {'CRPS':crps, 'RMSE':rmse}
   
    return results




def create_for_regions():
    
    region1 = {"Oelsnitz":328, 'Bad Elster':48 ,'Adorf':170 }
    
    region2 = {'Dohna':200, 'Lauenstein':76, 'Geising':26 }
    
    region3 = {'Niederoderwitz':29, 'Zittau':279, 'Seifhennersdorf':75 , 'Grossschoenau': 162}
    
    
    data1   = {key: [None,region1[key]]  for key in region1.keys()}
    for a in data1.keys():
            data1[a][0]=evaluate_for_catchment(a)   
    
    
    data2   = {key: [None,region2[key]]  for key in region2.keys()}
    for a in data2.keys():
            data2[a][0]=evaluate_for_catchment(a)
    
    
    data3   = {key: [None,region3[key]]  for key in region3.keys()}
    for a in data3.keys():
            data3[a][0]=evaluate_for_catchment(a)
    
    
    return [data1, data2, data3]



def plot_curves(results):
    # PLOTTING
    # TODO  change
    max_lead = 21
    leads = list(range(3,max_lead+1,3))
    
     
    region = 0
    regions = ['Weiße Elster region', 'Müglitz region', 'Eastern Saxony region']
    
    for res in results:
        
        for cat in res.keys():
            # create the figure
            sns.set_theme(style="whitegrid")
            fig, ax = plt.subplots(dpi=750)
            # plotting the results
            plt.plot(leads, cat['CRPS'],label = cat)
            ax.legend()
            plt.xticks(leads)
            ax.set_xlim(3,max_lead)        
            ax.set_title("Performance of ICOND2 Ensemble forecast  \n{}".format(regions[region]))
            ax.set_xlabel('Leadtime (hrs)')
            ax.set_ylabel('CRPS')
            
            plt.savefig("//vs-grp07.zih.tu-dresden.de/howa/work/students/Mohamed_Elghorab/results/CRPS_{}".format(regions[region]), bbox_inches = 'tight')        
        for cat in res.keys():
            # create the figure
            sns.set_theme(style="whitegrid")
            fig, ax = plt.subplots(dpi=750)
            # plotting the results
            plt.plot(leads, cat['RMSE'],label = cat)
            ax.legend()
            plt.xticks(leads)
            ax.set_xlim(3,max_lead)        
            ax.set_title("Performance of ICOND2 Ensemble forecast  \n{}".format(regions[region]))
            ax.set_xlabel('Leadtime (hrs)')
            ax.set_ylabel('RMSE')
            
            plt.savefig("//vs-grp07.zih.tu-dresden.de/howa/work/students/Mohamed_Elghorab/results/RMSE_{}".format(regions[region]), bbox_inches = 'tight')        
        
        region+=1
        
    return plt.show()
    
###############################################################################


    
df = create_for_regions()
plot_curves(df)




