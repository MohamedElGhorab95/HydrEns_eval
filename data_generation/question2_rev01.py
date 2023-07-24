# -*- coding: utf-8 -*-
"""
Created on Sun Jul  9 23:05:00 2023

@author: M Elghorab
"""

import time
import matplotlib.pyplot as plt
import seaborn as sns
from HydrEns_eval.engine.fr_entities_tools import *
from HydrEns_eval.engine.fr_Conttools import *
from HydrEns_eval.engine.fr_ROCtools import *
from tqdm import tqdm


def load_event_file(leadtime):
    source_file_path = "data_generation/nc_paths.txt"
    with open(source_file_path, 'r') as f:
        # Read all the lines in the file
        lines = f.readlines()
        
        netcdf = lines[1] 
        
        
     
        radar_file = netcdf+ 'fertig' + '/' + 'radRW_ICO.nc'
        deter_file = netcdf+ str(leadtime) + '/' + '{}hour_icond2.nc'.format(leadtime)
        ensem_file = netcdf+ str(leadtime) + '/' + '{}hour_icond2eps.nc'.format(leadtime)
        
        radar_det_cos = netcdf + 'fertig' + '/' + 'radRW_cosmod2.nc'
        radar_ens_cos = netcdf + 'fertig' + '/' + 'radRW_cosmod2eps.nc'
        det_cos = netcdf+ str(leadtime) + '/' + '{}hour_cosmod2.nc'.format(leadtime)
        ens_cos = netcdf+ str(leadtime) + '/' + '{}hour_cosmod2eps.nc'.format(leadtime)
        
    files = [radar_file, deter_file, ensem_file, radar_det_cos, radar_ens_cos, det_cos, ens_cos]
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
    det = Deterministic_run(locations[1]).gen_deterministic_field()
    fore = Ensemble_run(locations[2])
    
    # instantiating the cosmo fields
    rad_det = Observation(locations[3]).gen_observation_field().aggr_temporal(3)
    rad_ens = Observation(locations[4]).gen_observation_field().aggr_temporal(3)
    det_cos = Deterministic_run(locations[5]).gen_deterministic_field()
    for_cos = Ensemble_run(locations[6])
   
    print('cropping to: {}\ncalculaing areal averages'.format(catchment),flush= True)
   
    # cropping for selected catchment
    rad_c = rad.extract_by_shp(load_catchment(catchment)).avg_areal_prec()
    det_c = det.extract_by_shp(load_catchment(catchment)).avg_areal_prec()
    fore_c = fore.eps_extract_by_shp(load_catchment(catchment)).avg_areal_prec()
    
    rad_det_c = rad_det.extract_by_shp(load_catchment(catchment)).avg_areal_prec()
    rad_ens_c = rad_ens.extract_by_shp(load_catchment(catchment)).avg_areal_prec()
    det_cos_c = det_cos.extract_by_shp(load_catchment(catchment)).avg_areal_prec()
    for_cos_c = for_cos.eps_extract_by_shp(load_catchment(catchment)).avg_areal_prec()
    
    # garbage collection
    del rad, det ,fore, rad_det, rad_ens, det_cos, for_cos   
    
    library_icon = {'radar': rad_c,
               'ICOND2': det_c,
               'ICOND2EPS': fore_c,
               }
    
    library_cosmo = {'radar': rad_det_c, 'COSMOD2':det_cos_c}
    
    
    library_cosmo_eps = {'radar': rad_ens_c, 'COSMOD2EPS': for_cos_c}
    
    
    
    
    return library_icon, library_cosmo, library_cosmo_eps


# def calc_metric(observation, forecast , metric):
#     print("calculating performance metric..........",flush= True)
    
    
#     if metric == 'Frequency bias':
#         return CONT(observation, forecast, 3).fbias(50)
#     elif metric == 'Critical success index':
#         return CONT(observation, forecast, 3).csi(50)
#     elif metric == "Pierce's skill score":
#         return CONT(observation, forecast, 3).pss(50)
#     elif metric == 'Area under ROC curve':
#         return (float(ROC(observation, forecast).roc_auc(quantile=50)))
      




def evaluate_for_catchment(cats):

   
    # TODO change
    max_lead = 24
    threshold = 3
    variables = ['ICOND2', 'ICOND2EPS', 'COSMOD2', 'COSMOD2EPS']
    
    roc_auc_dic = {key: [None]*int(max_lead/3)  for key in variables}   
    con_acc_dic = {key: [None]*int(max_lead/3)  for key in variables}  
    con_pss_dic = {key: [None]*int(max_lead/3)  for key in variables}
    con_fbs_dic = {key: [None]*int(max_lead/3)  for key in variables}
    
    c = 0 # leadtime index counter
    
    for lead in range(3,max_lead+1,3):
        print('Catchment: {}\ncalculating for lead time: {}hrs'.format(cats,lead),flush= True)
        
        # loading the files
        files = load_event_file(lead)
        # creating the instances for a specific catchment
        lib_ic, lib_cos, lib_coseps = gen_for_catchment(cats,files)
        
        
        
        # Create a progress bar
        progress_bar = tqdm(total=4, desc='Calculating evaluation metrics.....')
        
        for x in variables[:2]:
            roc_auc_dic[x][c] = ROC(lib_ic['radar'], lib_ic[x]).roc_auc(quantile=50)
   
            con_acc_dic[x][c] = CONT(lib_ic['radar'], lib_ic[x],threshold).acc(50)*100
      
            con_pss_dic[x][c] = CONT(lib_ic['radar'], lib_ic[x],threshold).pss(50)*100
      
            con_fbs_dic[x][c] = CONT(lib_ic['radar'], lib_ic[x],threshold).fbias(50)
            progress_bar.update(1)
        
        roc_auc_dic[variables[2]][c] = ROC(lib_cos['radar'], lib_cos[variables[2]]).roc_auc(quantile=50)
   
        con_acc_dic[variables[2]][c] = CONT(lib_cos['radar'], lib_cos[variables[2]],threshold).acc(50)*100
  
        con_pss_dic[variables[2]][c] = CONT(lib_cos['radar'], lib_cos[variables[2]],threshold).pss(50)*100
  
        con_fbs_dic[variables[2]][c] = CONT(lib_cos['radar'], lib_cos[variables[2]],threshold).fbias(50)
        progress_bar.update(1)    
        
        roc_auc_dic[variables[-1]][c] = ROC(lib_coseps['radar'], lib_coseps[variables[-1]]).roc_auc(quantile=50)
   
        con_acc_dic[variables[-1]][c] = CONT(lib_coseps['radar'], lib_coseps[variables[-1]],threshold).acc(50)*100
  
        con_pss_dic[variables[-1]][c] = CONT(lib_coseps['radar'], lib_coseps[variables[-1]],threshold).pss(50)*100
  
        con_fbs_dic[variables[-1]][c] = CONT(lib_coseps['radar'], lib_coseps[variables[-1]],threshold).fbias(50)
        progress_bar.update(1)    
            
            
            
        # Close the progress bar
        progress_bar.close()
            
        del lib_ic, lib_cos, lib_coseps
            
       
        c+=1
    
    return {'Area under ROC curve':roc_auc_dic, 'Accuracy':con_acc_dic , 
            'PSS': con_pss_dic, 'Frequency bias': con_fbs_dic}



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
        plt.plot(leads, results[res]['ICOND2'],label='Main run | ICOND2' )
        plt.plot(leads, results[res]['ICOND2EPS'],label='Ensemble median | ICOND2' )
        
        plt.plot(leads, results[res]['COSMOD2'],'--',label='Main run | COSMOD2' )
        plt.plot(leads, results[res]['COSMOD2EPS'],'--',label='Ensemble median | COSMOD2' )
        
        

        ax.legend()
        plt.xticks(leads)
        ax.set_xlim(3,max_lead)
        ax.set_title("Skill score of different forecast products performance \n{} catchment".format(catchment_name))
        ax.set_xlabel('Leadtime (hrs)')
        ax.set_ylabel(res)
        
        print('exporting results to file')
        plt.savefig("//vs-grp07.zih.tu-dresden.de/howa/work/students/Mohamed_Elghorab/results/product analysis/{}/{}".format(catchment_name,res))        
    
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