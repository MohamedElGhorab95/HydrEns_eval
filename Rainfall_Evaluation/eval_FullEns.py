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
import os


def load_event_file(leadtime):
    # source_file_path = "in/nc_paths.txt"
    # with open(source_file_path, 'r') as f:
    #     # Read all the lines in the file
    #     lines = f.readlines()
        
    #     netcdf = lines[0] 
        
    netcdf = os.getcwd()+'/Data/NetCDFs'    
     
    radar_file = netcdf+ '/' + 'radRW_ICO.nc'
    # deter_file = netcdf+ str(leadtime) + '/' + '{}hour_icond2.nc'.format(leadtime)
    ensem_file = netcdf+ '/' + str(leadtime) + '/' + '{}hour_icond2eps.nc'.format(leadtime)
    
    # files = [radar_file, deter_file, ensem_file]
    files = [radar_file,  ensem_file]
    return files



def load_catchment(name):
    # loading the catchment of interest
    source_file_path = "in/paths.txt"
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
    max_lead = 24

    crps = []
    rmse = []
    disc = []
   
    for lead in range(3,max_lead+1,3):
        print('Catchment: {}\ncalculating for lead time: {}hrs'.format(cats,lead),flush= True)
        
        # loading the files
        files = load_event_file(lead)
        # creating the instances for a specific catchment
        library = gen_for_catchment(cats,files)

        crps.append(Full_Ens(library['radar'], library['ensemble_run']).crps())
        rmse.append(Full_Ens(library['radar'], library['ensemble_run']).rmse())
        disc.append(Full_Ens(library['radar'], library['ensemble_run']).disc_dia())    
        del library
    
    results = {'CRPS':crps, 'RMSE':rmse, 'Discrimination':disc}
    # results = {'CRPS':crps}
   
    return results




def plot_curves(results):
    # PLOTTING
    # TODO  change
    max_lead = 24
    leads = list(range(3,max_lead+1,3))
    
     
    region = 0
    regions = ['Weiße Elster region', 'Müglitz region', 'Mandau region']
    # regions = [ 'Weiße Elster region', 'Müglitz region']
    
    for res in results:
        
        # create the figure
        sns.set_theme(style="whitegrid")
        fig, ax = plt.subplots(dpi=750)
        for cat in res.keys():  
            # plotting the results
            plt.plot(leads, res[cat][0]['CRPS'],label = cat)
        ax.legend()
        plt.xticks(leads)
        ax.set_xlim(3,max_lead)    
        ax.set_ylim(0.15, 0.45) 
        ax.set_title("Performance of ICOND2 Ensemble forecast  \n{}".format(regions[region]))
        ax.set_xlabel('Leadtime (hrs)')
        ax.set_ylabel('CRPS')
        
        # plt.savefig("E:/results/CRPS_{}".format(regions[region]), bbox_inches = 'tight')        
        
        # create the figure
        sns.set_theme(style="whitegrid")
        fig, ax = plt.subplots(dpi=750)
        for cat in res.keys():
            
            # plotting the results
            plt.plot(leads, res[cat][0]['RMSE'],label = cat)
        ax.legend()
        plt.xticks(leads)
        ax.set_xlim(3,max_lead) 
        ax.set_ylim(1.0,2.8)
        ax.set_title("Performance of ICOND2 Ensemble forecast  \n{}".format(regions[region]))
        ax.set_xlabel('Leadtime (hrs)')
        ax.set_ylabel('RMSE')
        
        # plt.savefig("E:/results/RMSE_{}".format(regions[region]), bbox_inches = 'tight')        
        
        
        # create the figure
        sns.set_theme(style="whitegrid")
        fig, ax = plt.subplots(dpi=750)
        for cat in res.keys():
            
            # plotting the results
            plt.plot(leads, res[cat][0]['Discrimination'],label = cat)
        ax.legend()
        plt.xticks(leads)
        ax.set_xlim(3,max_lead) 
        ax.set_ylim(1.2,1.65)
        ax.set_title("Performance of ICOND2 Ensemble forecast  \n{}".format(regions[region]))
        ax.set_xlabel('Leadtime (hrs)')
        ax.set_ylabel('Discrimination')
        
        # plt.savefig("E:/results/Disc_{}".format(regions[region]), bbox_inches = 'tight')        
        
        
        region+=1
        
    return plt.show()
    
###############################################################################

################################################################################################

if __name__ == '__main__':
################################################################################################
    
    

    def create_for_regions():
        
        region1 = {"Oelsnitz":328 ,'Adorf':170 , 'Bad Elster':48}
        
        region1 = dict(sorted(region1.items()))
        
        region2 = {'Dohna':200, 'Lauenstein':76, 'Geising':26 }
        region2 = dict(sorted(region2.items()))
        
        region3 = { 'Zittau':279, 'Grossschoenau': 162, 'Seifhennersdorf':75,'Niederoderwitz':29 }
        region3 = dict(sorted(region3.items()))
        
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



    
    df = create_for_regions()
    plot_curves(df)




