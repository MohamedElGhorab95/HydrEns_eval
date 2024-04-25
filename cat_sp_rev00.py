# -*- coding: utf-8 -*-
"""
Created on Tue Jun  6 17:36:05 2023

@author: M Elghorab
"""
from engine.fr_entities_tools import *
from engine.fr_Conttools import *
from engine.fr_ROCtools import *
import matplotlib.pyplot as plt
import seaborn as sns
import xarray as xr
import numpy as np
import time
from tqdm import tqdm


def load_event_file(leadtime):
    source_file_path = "in/nc_paths.txt"
    with open(source_file_path, 'r') as f:
        # Read all the lines in the file
        lines = f.readlines()
        
        netcdf = lines[0] 
        
        
     
        radar_file = netcdf+ 'fertig' + '/' + 'radRW_ICO.nc'
        deter_file = netcdf+ str(leadtime) + '/' + '{}hour_icond2.nc'.format(leadtime)
        ensem_file = netcdf+ str(leadtime) + '/' + '{}hour_icond2eps.nc'.format(leadtime)
        
    files = [radar_file, deter_file, ensem_file]
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
    rad = Observation(locations[0]).gen_observation_field()
    det = Deterministic_run(locations[1]).gen_deterministic_field()
    eps = Ensemble_run(locations[2])
   
    print('cropping to: {}'.format(catchment),flush= True)
    # print(f"cropping to: {catchment}")
    # cropping for selected catchment
    rad_c = rad.extract_by_shp(load_catchment(catchment)).aggr_temporal(3)
    det_c = det.extract_by_shp(load_catchment(catchment))
    eps_c = eps.eps_extract_by_shp(load_catchment(catchment))
    # garbage collection
    del rad, det, eps
   
    # average areal precipitation
    print('calculaing areal averages',flush= True)
    # print(f"calculaing areal averages")
    rad_c = rad_c.avg_areal_prec()
    det_c = det_c.avg_areal_prec()
    eps_c = eps_c.avg_areal_prec()
   
    print('generating quantiles',flush= True)
    # print(f"generating quantiles")
    # generating quantiles
    qs = [10,25,75,90,'mean']
    eps_q= []
    for i in qs:
        x = eps_c.gen_quantiles(i)    
        eps_q.append(x)
        del x
    
    
    library = {'radar': rad_c,
               'deterministic_run': det_c,
               'ensemble_q10': eps_q[0],
               'ensemble_q25': eps_q[1],
               'ensemble_q75': eps_q[2],
               'ensemble_q90': eps_q[3],
               'ensemble_mean':eps_q[4]}
   
    del rad_c, det_c, eps_q
    return library



def calc_metric(event_dictionary, metric, threshold = 2.5):
    print("calculating {}..........".format(metric),flush= True)
    
    results = []
 
    
    env = event_dictionary
    
    names = list(env.keys())
    names.remove('radar')
    
    for t in names:
        
            
            if metric == "cont.hits":
                results.append(CONT(env['radar'], env[t], threshold).hits())
            elif metric == "cont.misses":
                results.append(CONT(env['radar'], env[t], threshold).misses())
            elif metric == "cont.fal":
                results.append(CONT(env['radar'], env[t], threshold).false_alarms())
            elif metric == 'accuracy':
                results.append(CONT(env['radar'], env[t], threshold).acc()*100)
            elif metric == 'fbias':
                results.append(CONT(env['radar'], env[t], threshold).fbias())
    if metric == 'cont':
        return res_cont
    else:
        return results







def evaluate_for_catchment(cats):

    variables = [
               'main_run',
               'ensemble_q10',
               'ensemble_q25',
               'ensemble_q75',
               'ensemble_q90',
               'ensemble_mean']
    # TODO change
    max_lead = 24
    
    # initialize the result dictionaries
    con_hit_dic = {key: [None]*int(max_lead/3)  for key in variables}
    con_mis_dic = {key: [None]*int(max_lead/3)  for key in variables}
    con_fal_dic = {key: [None]*int(max_lead/3)  for key in variables}
    con_acc_dic = {key: [None]*int(max_lead/3)  for key in variables}
    con_fbias_dic = {key: [None]*int(max_lead/3)  for key in variables}
    
    
    
    
    c = 0 # leadtime index counter
    
    for lead in range(3,max_lead+1,3):
        print('Catchment: {}\ncalculating for lead time: {}hrs'.format(cats,lead),flush= True)
        
        # loading the files
        files = load_event_file(lead)
        # creating the instances for a specific catchment
        library = gen_for_catchment(cats,files)
        
        # calculate metrics
        
        con_hit = calc_metric(library, "cont.hits")
        con_mis = calc_metric(library, "cont.misses")
        con_fal = calc_metric(library, "cont.fal")
        con_acc = calc_metric(library, "accuracy")
        con_fbias = calc_metric(library, "fbias")
        
        del library, files
        
        # fill the dictionaries with values
        for var in variables:
            
            con_hit_dic[var][c]= con_hit[variables.index(var)]
            con_mis_dic[var][c]= con_mis[variables.index(var)]
            con_fal_dic[var][c]= con_fal[variables.index(var)]
            con_acc_dic[var][c]= con_acc[variables.index(var)]
            con_fbias_dic[var][c]= con_fbias[variables.index(var)]
        
        
        c+=1
    
    return {'Hits':con_hit_dic,
            'Misses': con_mis_dic, 'False Alarms': con_fal_dic, 
            'Accuracy':con_acc_dic, 'Frequency bias':con_fbias_dic}














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
        plt.plot(leads, results[res]['ensemble_q10'],ls='--',label='Q10%' )
        plt.plot(leads, results[res]['ensemble_q90'],ls='--',label='Q90%' )
        plt.plot(leads, results[res]['ensemble_q25'],ls='-.',label='Q25%' )
        plt.plot(leads, results[res]['ensemble_q75'],ls='-.',label='Q75%' )
        
        

        ax.legend()
        plt.xticks(leads)
        ax.set_xlim(3,max_lead)
        
        
        if res == 'Accuracy'            : ax.set_ylim(90,100)
        elif res == 'Frequency bias'      : ax.set_ylim(0.5,1.4)
        
        
        
        
        ax.set_title("Performance of ICOND2 forecast  \nGauge name:{}".format(catchment_name))
        ax.set_xlabel('Leadtime (hrs)')
        ax.set_ylabel(res)
        
        print('exporting results to file')
        plt.savefig("E:/results/{}/{}".format(catchment_name,res))        
    
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