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
    source_file_path = "data_generation/nc_paths.txt"
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
    source_file_path = "data_generation/paths.txt"
    with open(source_file_path, 'r') as f:
        # Read all the lines in the file
        lines = f.readlines()
        if name == 'mandau':
            return lines[2].strip()
        if name == 'mugliz':
            return lines[4].strip()
        if name == 'weisse':
            return lines[6].strip()



def gen_for_catchment(catchment, locations):
    
    print('generating instances..........')
    # progress bar graphics
   
    
    # creating instances 
    rad = Observation(locations[0]).gen_observation_field().aggr_temporal(3)
    det = Deterministic_run(locations[1]).gen_deterministic_field()
    # eps = Ensemble_run(locations[2])
   
    print('cropping to: {}'.format(catchment))
    # cropping for selected catchment
    rad_c = rad.extract_by_shp(load_catchment(catchment))
    det_c = det.extract_by_shp(load_catchment(catchment))
    # eps_c = eps.eps_extract_by_shp(load_catchment(catchment))
    # garbage collection
    # del rad, det, eps
   
    # average areal precipitation
    print('calculaing areal averages')
    rad_c = rad_c.avg_areal_prec()
    det_c = det_c.avg_areal_prec()
    # eps_c = eps_c.avg_areal_prec()
   
    print('generating quantiles')
    # generating quantiles
    # qs = [10,25,75,90,'mean']
    # eps_q = [eps_c.gen_quantiles(i) for i in qs]
    
    
    # library = {'radar': rad_c,
    #            'deterministic_run': det_c,
    #            'ensemble_q10': eps_q[0],
    #            'ensemble_q25': eps_q[1],
    #            'ensemble_q75': eps_q[2],
    #            'ensemble_q90': eps_q[3],
    #            'ensemble_mean':eps_q[4]}
    library = {'radar': rad_c,
               'deterministic_run': det_c,
               }
   
    
    return library



def calc_metric(event_dictionary, metric, threshold = 4):
    
    print('calculating {}..........'.format(metric))
    results = []
 
    
    env = event_dictionary
    
    for t in list(env.keys())[1:]:
        
            if metric == "ROC_auc":
                results.append(float(ROC(env['radar'], env[t]).roc_auc()))
                # a = ROC(env['radar'][0], env['radar'][idx])
                # a.roc_auc()
                # a.plot_roc()
            elif metric == "success_ratio":
                results.append(CONT(env['radar'], env[t], threshold).sr()*100)
            elif metric == "CSI":
                results.append(CONT(env['radar'], env[t], threshold).csi()*100)
            elif metric == "FAR":
                results.append(CONT(env['radar'], env[t], threshold).far()*100)
            elif metric == "PSS":
                results.append(CONT(env['radar'], env[t], threshold).pss()*100)
            elif metric == "cont.hits":
                results.append(CONT(env['radar'], env[t], threshold).hits())
            elif metric == "cont.misses":
                results.append(CONT(env['radar'], env[t], threshold).misses())
            elif metric == "cont.fal":
                results.append(CONT(env['radar'], env[t], threshold).false_alarms())
            elif metric == 'table':
                print('{} Contingency table \nthreshold={}'.format(t,threshold))
                CONT(env['radar'], env[t], threshold).print_cont_tab()
    if metric == 'cont':
        return res_cont
    else:
        return results







def evaluate_for_catchment(cats):

    variables = [
               'main_run'
               ]
    # TODO change
    max_lead = 24
    
    # initialize the result dictionaries
    roc_auc_dic = {key: [None]*int(max_lead/3)  for key in variables}
    cont_sr_dic = {key: [None]*int(max_lead/3)  for key in variables}
    con_csi_dic = {key: [None]*int(max_lead/3)  for key in variables}
    con_far_dic = {key: [None]*int(max_lead/3)  for key in variables}
    con_pss_dic = {key: [None]*int(max_lead/3)  for key in variables}
    con_hit_dic = {key: [None]*int(max_lead/3)  for key in variables}
    con_mis_dic = {key: [None]*int(max_lead/3)  for key in variables}
    con_fal_dic = {key: [None]*int(max_lead/3)  for key in variables}
    
    
    
    
    c = 0 # leadtime index counter
    
    for lead in range(3,max_lead+1,3):
        print('Catchment: {}\ncalculating for lead time: {}hrs'.format(cats,lead))
        # loading the files
        files = load_event_file(lead)
        # creating the instances for a specific catchment
        library = gen_for_catchment(cats,files)
        
        # calculate metrics
        roc_auc = calc_metric(library, "ROC_auc")
        cont_sr = calc_metric(library, "success_ratio")
        con_csi = calc_metric(library, "CSI")
        con_far = calc_metric(library, "FAR")
        con_pss = calc_metric(library, "PSS")
        con_hit = calc_metric(library, "cont.hits")
        con_mis = calc_metric(library, "cont.misses")
        con_fal = calc_metric(library, "cont.fal")
        
        
        # fill the dictionaries with values
        for var in variables:
            roc_auc_dic[var][c]= roc_auc[variables.index(var)]
            cont_sr_dic[var][c]= cont_sr[variables.index(var)]
            con_csi_dic[var][c]= con_csi[variables.index(var)]
            con_far_dic[var][c]= con_far[variables.index(var)]
            con_pss_dic[var][c]= con_pss[variables.index(var)]
            con_hit_dic[var][c]= con_hit[variables.index(var)]
            con_mis_dic[var][c]= con_mis[variables.index(var)]
            con_fal_dic[var][c]= con_fal[variables.index(var)]
        
        
        c+=1
    
    return {'Area under ROC curve':roc_auc_dic, 'Success ratio': cont_sr_dic,
            'CSI':con_csi_dic , 'False alarm ratio':con_far_dic,
            'PSS': con_pss_dic, 'Hits':con_hit_dic,
            'Misses': con_mis_dic, 'False Alarms': con_fal_dic}, library













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
        # plt.plot(leads, results[res]['ensemble_mean'],label='Ensemble mean' )
        
        # plt.fill_between(leads, results[res]['ensemble_q10'], 
        #                  results[res]['ensemble_q90'],
        #                  alpha=0.4, color='mediumturquoise',label='Quantile 10-90%')
        
        # plt.fill_between(leads, results[res]['ensemble_q25'],
        #                  results[res]['ensemble_q75'],
        #                  alpha=0.6, color='darkcyan',label='Quantile 25-75%')

        ax.legend()
        plt.xticks(leads)
        ax.set_title("Skill score of ICOND2 forecast performance \n{} catchment".format(catchment_name))
        ax.set_xlabel('Leadtime (hrs)')
        ax.set_ylabel(res)
        
        print('exporting results to file')
        plt.savefig("//vs-grp07.zih.tu-dresden.de/howa/work/students/Mohamed_Elghorab/results/{}/{}".format(catchment_name,res))        
    
    return plt.show()
    


###############################################################################

tic = time.time()
mugliz_results, lib_mu = evaluate_for_catchment("mugliz")
print((time.time()-tic)/3600) # operation time in hours
# plot_curves(mugliz_results, "Mügliz")

# tic = time.time()
# mugliz_results, lib_mu = evaluate_for_catchment("mandau")
# print((time.time()-tic)/3600) # operation time in hours
# plot_curves(mugliz_results, "Mandau")



# tic = time.time()
# mugliz_results, lib_mu = evaluate_for_catchment("weisse")
# print((time.time()-tic)/3600) # operation time in hours
# plot_curves(mugliz_results, "Weiße Elster")



list(lib_mu.keys())

# a = lib_mu['radar']
# b = lib_mu['deterministic_run']

# CONT(a, b, 15).misses()




# unique_timesteps1 = a.average.time.values

# # Get the unique timesteps in data_array2
# unique_timesteps2 = b.average.time.values

# # Find the different timesteps
# different_timesteps = set(unique_timesteps1) ^ set(unique_timesteps2)

# # Print the different timesteps
# print("Different timesteps:")
# for timestep in different_timesteps:
#     print(timestep in unique_timesteps2)













