# -*- coding: utf-8 -*-
"""
Created on Sun Jul  9 23:05:00 2023

@author: M Elghorab
"""

import time
import matplotlib.pyplot as plt
import seaborn as sns
from engine.fr_entities_tools import *
from engine.fr_Conttools import *
from engine.fr_ROCtools import *
from tqdm import tqdm


def load_event_file(leadtime):
    source_file_path = "data_generation/nc_paths.txt"
    with open(source_file_path, 'r') as f:
        # Read all the lines in the file
        lines = f.readlines()
        
        netcdf = lines[1] 
        
        
        radar_det_cos = netcdf + 'fertig' + '/' + 'radRW_cosmod2.nc'
        radar_ens_cos = netcdf + 'fertig' + '/' + 'radRW_cosmod2eps.nc'
        det_cos = netcdf+ str(leadtime) + '/' + '{}hour_cosmod2.nc'.format(leadtime)
        ens_cos = netcdf+ str(leadtime) + '/' + '{}hour_cosmod2eps.nc'.format(leadtime)
        
    files = [radar_det_cos, radar_ens_cos, det_cos, ens_cos]
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
   
    
   
    
    # instantiating the cosmo fields
    rad_det = Observation(locations[0]).gen_observation_field()
    rad_ens = Observation(locations[1]).gen_observation_field()
    det_cos = Deterministic_run(locations[2]).gen_deterministic_field()
    for_cos = Ensemble_run(locations[3])
   
    print('cropping to: {}\ncalculaing areal averages'.format(catchment),flush= True)
   
    
    
    rad_det_c = rad_det.extract_by_shp(load_catchment(catchment)).avg_areal_prec().aggr_temporal(3)
    rad_ens_c = rad_ens.extract_by_shp(load_catchment(catchment)).avg_areal_prec().aggr_temporal(3)
    det_cos_c = det_cos.extract_by_shp(load_catchment(catchment)).avg_areal_prec()
    for_cos_c = for_cos.eps_extract_by_shp(load_catchment(catchment)).avg_areal_prec()
    
    # garbage collection
    del rad_det, rad_ens, det_cos, for_cos   
    
    
    
    library_cosmo = {'radar': rad_det_c, 'COSMOD2':det_cos_c}
    
    
    library_cosmo_eps = {'radar': rad_ens_c, 'COSMOD2EPS': for_cos_c}
    
    
    
    
    return library_cosmo, library_cosmo_eps



      




def evaluate_for_catchment(cats):
    
    variables = {
               'main_run':['COSMOD2',None],
               'ensemble_q10':['COSMOD2EPS',10],
               'ensemble_q25':['COSMOD2EPS',25],
               'ensemble_q75':['COSMOD2EPS',75],
               'ensemble_q90':['COSMOD2EPS',90],
               'ensemble_mean':['COSMOD2EPS','mean']} 
   
    # TODO change
    max_lead = 24
    
    roc_auc_dic = {key: [None]*int(max_lead/3)  for key in variables}   
    
    
    c = 0 # leadtime index counter
    
    for lead in range(3,max_lead+1,3):
        print('Catchment: {}\ncalculating for lead time: {}hrs'.format(cats,lead),flush= True)
        
        # loading the files
        files = load_event_file(lead)
        # creating the instances for a specific catchment
        lib_cos, lib_coseps = gen_for_catchment(cats,files)
        
        
        
        roc_auc_dic['main_run'][c] = ROC(lib_cos['radar'], lib_cos['COSMOD2']).roc_auc(2.5)
        
        for x in list(variables.keys())[1:]:
            roc_auc_dic[x][c] = ROC(lib_coseps['radar'], lib_coseps[variables[x][0]],variables[x][1]).roc_auc(2.5)
            
            
            
            
       
        c+=1
    
    return {'Area under ROC curve':roc_auc_dic}



def plot_curves(results, catchment_name):
    # PLOTTING
    # TODO  change
    max_lead = 24
    leads = list(range(3,max_lead+1,3))
    
    for res in list(results.keys()):
        # create the figure
        sns.set_theme(style="whitegrid")
        fig, ax = plt.subplots(dpi=750)
        
        
        
        
        # plotting the results
        plt.plot(leads, results[res]['main_run'],label='Main run' )
        plt.plot(leads, results[res]['ensemble_mean'],label='Ensemble mean' )
        
        plt.fill_between(leads, results[res]['ensemble_q10'], 
                         results[res]['ensemble_q90'],
                         alpha=0.4, color='mediumturquoise',label='Quantile 10-90%')
        
        plt.fill_between(leads, results[res]['ensemble_q25'],
                         results[res]['ensemble_q75'],
                         alpha=0.6, color='darkcyan',label='Quantile 25-75%')
        

        ax.legend()
        plt.xticks(leads)
        ax.set_xlim(3,max_lead)
        ax.set_title("Performance of COSMOD2 forecast \n{} catchment".format(catchment_name))
        ax.set_xlabel('Leadtime (hrs)')
        ax.set_ylabel(res)
        
        print('exporting results to file')
        plt.savefig("//vs-grp07.zih.tu-dresden.de/howa/work/students/Mohamed_Elghorab/results/product analysis/{}/ROC_eval".format(catchment_name), bbox_inches = 'tight')        
    
    return plt.show()


###############################################################################

def evaluate_single(cname):
    tic = time.time()
    results = evaluate_for_catchment(cname)
    print((time.time()-tic)/3600) # operation time in hours
    plot_curves(results, cname)


cats = ['Dohna', "Oelsnitz", 'Zittau', 'Adorf', 'Geising', 'Grossschoenau',
        'Bad Elster', 'Lauenstein', 'Niederoderwitz','Seifhennersdorf']

for cat in cats:
    evaluate_single(cat)