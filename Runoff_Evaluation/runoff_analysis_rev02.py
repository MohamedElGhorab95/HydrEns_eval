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




def calc_metric(event_dictionary, metric, threshold):
    print("calculating {}..........".format(metric),flush= True)
    
    results = []
 
    env = event_dictionary
    
    names = list(env.keys())
    names.remove('gauge')
    
    for t in names:
        
            
            if metric == "cont.hits":
                results.append(CONT(env['gauge'], env[t], threshold).hits())
            elif metric == "cont.misses":
                results.append(CONT(env['gauge'], env[t], threshold).misses())
            elif metric == "cont.fal":
                results.append(CONT(env['gauge'], env[t], threshold).false_alarms())
            elif metric == 'accuracy':
                results.append(CONT(env['gauge'], env[t], threshold).acc()*100)
            elif metric == 'roc':
                r = ROC(env['gauge'], env[t])
                results.append(r.roc_auc(threshold))
   
    
    return results



def evaluate_for_cat(catchment, th):
    
    threshold = 0.25*th
    
    if catchment == 'Oelsnitz':
        obs = R_Observation('C:/Users/User/Downloads/Oelsnitz.csv')
    elif catchment == 'Adorf':
        obs = R_Observation('C:/Users/User/Downloads/Adorf.csv')
    elif catchment == 'Bad Elster':
        obs = R_Observation('C:/Users/User/Downloads/Bad_Elster.csv')
    
    variables = [
               'ensemble_q10',
               'ensemble_q25',
               'ensemble_q50',
               'ensemble_q75',
               'ensemble_q90',
               'ensemble_mean']
    
    
    # initialize the result dictionaries
    con_hit_dic = {key: [None]*int(24/3)  for key in variables}
    con_mis_dic = {key: [None]*int(24/3)  for key in variables}
    con_fal_dic = {key: [None]*int(24/3)  for key in variables}
    con_acc_dic = {key: [None]*int(24/3)  for key in variables}
    roc_auc_dic = {key: [None]*int(24/3)  for key in variables}
    
    
    
    auc = []
    crps = []
    rmse = []
    
    c = 0
    
    for lead in np.arange(3,25,3):
        
        if catchment == 'Oelsnitz':
            quans = xr.open_dataset('//vs-grp07.zih.tu-dresden.de/howa/work/students/Mohamed_Elghorab/netCDFs/Runoff/{}hrs_leadtime__5661371_Oelsnitz_WeisseElster_data.nc'.format(lead))
            full = xr.open_dataset('//vs-grp07.zih.tu-dresden.de/howa/work/students/Mohamed_Elghorab/netCDFs/Runoff/{}hrs_leadtime__5661371_Oelsnitz_WeisseElster_data_EPS.nc'.format(lead))
        elif catchment == 'Adorf':
            quans = xr.open_dataset('//vs-grp07.zih.tu-dresden.de/howa/work/students/Mohamed_Elghorab/netCDFs/Runoff/{}hrs_leadtime__5661311_Adorf_WeisseElster_data.nc'.format(lead))
            full = xr.open_dataset('//vs-grp07.zih.tu-dresden.de/howa/work/students/Mohamed_Elghorab/netCDFs/Runoff/{}hrs_leadtime__5661311_Adorf_WeisseElster_data_EPS.nc'.format(lead))
        elif catchment == 'Bad Elster':
            quans = xr.open_dataset('//vs-grp07.zih.tu-dresden.de/howa/work/students/Mohamed_Elghorab/netCDFs/Runoff/{}hrs_leadtime__56611313_BadElster_WeisseElster_data.nc'.format(lead))
            full = xr.open_dataset('//vs-grp07.zih.tu-dresden.de/howa/work/students/Mohamed_Elghorab/netCDFs/Runoff/{}hrs_leadtime__56611313_BadElster_WeisseElster_data_EPS.nc'.format(lead))
            
            
        library = {'gauge': obs,
                   'ensemble_q10': quans.Q_q10,
                   'ensemble_q25': quans.Q_q25,
                   'ensemble_q50': quans.Q_q50,
                   'ensemble_q75': quans.Q_q75,
                   'ensemble_q90': quans.Q_q90,
                   'ensemble_mean':quans.Q_mean}
    
        
        con_hit = calc_metric(library, "cont.hits", threshold)
        con_mis = calc_metric(library, "cont.misses", threshold)
        con_fal = calc_metric(library, "cont.fal", threshold)
        con_acc = calc_metric(library, "accuracy", threshold)
        con_auc = calc_metric(library, "roc", threshold)
        auc.append(ROC_totens(obs, full, threshold).area)
        crps.append(Full_Ens(obs.fr, full).crps())
        rmse.append(Full_Ens(obs.fr, full).rmse())
        # fill the dictionaries with values
        for var in variables:
            
            con_hit_dic[var][c]= con_hit[variables.index(var)]
            con_mis_dic[var][c]= con_mis[variables.index(var)]
            con_fal_dic[var][c]= con_fal[variables.index(var)]
            con_acc_dic[var][c]= con_acc[variables.index(var)]
            roc_auc_dic[var][c]= con_auc[variables.index(var)]
        
        
        c+=1
    
    return {'Hits':con_hit_dic,
            'Misses': con_mis_dic, 'False Alarms': con_fal_dic, 
            'Accuracy':con_acc_dic, 'Area under ROC curve':roc_auc_dic, 'Area under ROC curve | Ensemble':auc, 'CRPS':crps, 'RMSE':rmse}




def plot_curves(results, catchment_name):
    # PLOTTING
    # TODO  change
    max_lead = 24
    leads = list(range(3,max_lead+1,3))
    
    for res in list(results.keys()):
        # create the figure
        sns.set_theme(style="whitegrid")
        fig, ax = plt.subplots(dpi=750)
        
        if res in ['CRPS','Area under ROC curve | Ensemble'] :
            plt.plot(leads, results[res])
            ax.legend()
            plt.xticks(leads)
            ax.set_xlim(3,max_lead)
            ax.set_title("Performance of runoff forecast  \n{} catchment".format(catchment_name))
            ax.set_xlabel('Leadtime (hrs)')
            ax.set_ylabel(res)
        else:
            # plotting the results
            plt.plot(leads, results[res]['ensemble_q10'],ls='--',label='Q10%' )
            plt.plot(leads, results[res]['ensemble_q25'],ls='-.',label='Q25%' )
            # plt.plot(leads, results[res]['ensemble_q50'],label='Ensemble median' )
            plt.plot(leads, results[res]['ensemble_q75'],ls='-.',label='Q75%' )
            plt.plot(leads, results[res]['ensemble_q90'],ls='--',label='Q90%' )
            plt.plot(leads, results[res]['ensemble_mean'],label='Ensemble mean' )
            
            
    
            ax.legend()
            plt.xticks(leads)
            ax.set_xlim(3,max_lead)
            
            
            if res == 'Accuracy'            : ax.set_ylim(75,100)
            if res in ['Hits',  'Misses',  'False Alarms']  : ax.yaxis.set_major_locator(MaxNLocator(integer=True))
            
            
            ax.set_title("Performance of runoff forecast  \n{} catchment".format(catchment_name))
            ax.set_xlabel('Leadtime (hrs)')
            ax.set_ylabel(res)
            
        print('exporting results to file')
        # plt.savefig("//vs-grp07.zih.tu-dresden.de/howa/work/students/Mohamed_Elghorab/results/{}/{}".format(catchment_name,res))        
    
    return plt.show()



cat_library = {'Oelsnitz':25.6 , 'Adorf':13.3, 'Bad Elster': 3.72}

for c in cat_library.keys():
    
    res = evaluate_for_cat(c,cat_library[c])
    plot_curves(res, c)

    

