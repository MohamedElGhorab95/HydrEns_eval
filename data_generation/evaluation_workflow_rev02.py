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


def load_event_file(leadtime):
    source_file_path = "data_generation/nc_paths.txt"
    with open(source_file_path, 'r') as f:
        # Read all the lines in the file
        lines = f.readlines()
        
        netcdf = lines[0] 
        
        
        # TODO fix file names
        radar_file = netcdf+ str(leadtime) + '/' + 'radRW_ICO.nc'
        deter_file = netcdf+ str(leadtime) + '/' + 'icond2_ev.nc' 
        ensem_file = netcdf+ str(leadtime) + '/' + 'icond2eps.nc'
        
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
    
    # creating instances 
    rad = Observation(locations[0]).gen_observation_field()
    det = Deterministic_run(locations[1]).gen_deterministic_field()
    eps = Ensemble_run(locations[2])
    
    
    # cropping for selected catchment
    rad_c = rad.extract_by_shp(load_catchment(catchment))
    det_c = det.extract_by_shp(load_catchment(catchment))
    eps_c = eps.eps_extract_by_shp(load_catchment(catchment))
    
    # average areal precipitation
    rad_c = rad_c.avg_areal_prec()
    det_c = det_c.avg_areal_prec()
    eps_c = eps_c.avg_areal_prec()
    
    # generating quantiles
    qs = [10,25,75,90,'mean']
    eps_q = [eps_c.gen_quantiles(i) for i in qs]
    
    
    library = {'radar': rad_c,
               'deterministic_run': det_c,
               'ensemble_q10': eps_q[0],
               'ensemble_q25': eps_q[1],
               'ensemble_q75': eps_q[2],
               'ensemble_q90': eps_q[3],
               'ensemble_mean':eps_q[4]}
    
    
    return library



def calc_metric(event_dictionary, metric, threshold = 10):
    
    results = []
    res_cont = [[],[],[]]
    
    for t in list(event_dictionary.keys())[1:]:
        
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
            elif metric == "cont":
                res_cont[0].append(CONT(env['radar'], env[t], threshold).hits())
                res_cont[1].append(CONT(env['radar'], env[t], threshold).misses())
                res_cont[2].append(CONT(env['radar'], env[t], threshold).false_alarms())
            elif metric == 'table':
                print('{} Contingency table \nthreshold={}'.format(t,threshold))
                CONT(env['radar'], env[t], threshold).print_cont_tab()
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
    
    
    
    # initialize the result dictionaries
    roc_auc_dic = {key: [None]*8  for key in variables}
    cont_sr_dic = {key: [None]*8  for key in variables}
    con_csi_dic = {key: [None]*8  for key in variables}
    con_far_dic = {key: [None]*8  for key in variables}
    con_pss_dic = {key: [None]*8  for key in variables}
    con_cont_dic = {key: [None]*8  for key in variables}
    
    
    
    
    c = 0 # leadtime index counter
    
    for lead in range(3,25,3):
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
        con_cont = calc_metric(library, "cont")
        
        
        # fill the dictionaries with values
        for var in variables:
            roc_auc_dic[var][c]= roc_auc[variables.index(var)]
            cont_sr_dic[var][c]= cont_sr[variables.index(var)]
            con_csi_dic[var][c]= con_csi[variables.index(var)]
            con_far_dic[var][c]= con_far[variables.index(var)]
            con_pss_dic[var][c]= con_pss[variables.index(var)]
            con_cont_dic[var][c]= con_cont[variables.index(var)]
        
        
        c+=1
    
    return {'Area under ROC curve':roc_auc_dic, 'Success ratio': cont_sr_dic,
            'CSI':con_csi_dic , 'False alarm ratio':con_far_dic,
            'PSS': con_pss_dic, 'contingency table':con_cont_dic}













def plot_curve(results, catchment_name):
    # PLOTTING
    
    
    leads = list(range(3,25,3))
    
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
        ax.set_title("Skill score of performance | {} catchment".format(catchment_name))
        ax.set_xlabel('Leadtime (hrs)')
        ax.set_ylabel(res)
    
    return plt.show()
    

def plot_TS(event_dictionary):
    
    leads = list(event_dictionary.keys())
    for t in leads:

        sns.set_theme(style="whitegrid")
        fig, ax = plt.subplots(dpi=750)
        
        event_dictionary[t][0].plot(label='Radar')
        event_dictionary[t][1].plot(label='ICOND2')
        event_dictionary[t][2].plot(label='ICOND2EPS|q=95')
        event_dictionary[t][3].plot(label='ICOND2EPS|mean')
        # for i in range(len(thresholds)):
        #     ax.hlines(thresholds[i], xmin=Lead_3[0].time[0], xmax=Lead_3[0].time[-1], colors='red', linestyles='dashed')
        
        
        ax.legend()
        ax.set_title("MAP | Müglitz catchment \nLeadtime = {}hrs".format(t))
        ax.set_ylabel('Rainfall (mm/3hrs)')
        plt.show()
        







