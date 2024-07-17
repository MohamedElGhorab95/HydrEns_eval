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
import warnings
warnings.filterwarnings("ignore")
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
    ensem_file = netcdf + '/' + str(leadtime) + '/' + '{}hour_icond2eps.nc'.format(leadtime)
    
    
    rad_cos = netcdf + '/' +  'radRW_COS.nc'
    # det_cos = netcdf+ str(leadtime) + '/' + '{}hour_cosmod2.nc'.format(leadtime)
    ens_cos = netcdf+ '/' + str(leadtime) + '/' + '{}hour_cosmod2eps.nc'.format(leadtime)
        
    # files = [radar_file, deter_file, ensem_file, rad_cos, det_cos, ens_cos]
    files = [radar_file, ensem_file, rad_cos, ens_cos]
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
    fore = Ensemble_run(locations[1])
    
    # instantiating the cosmo fields
    rad_cos = Observation(locations[2]).gen_observation_field()
    for_cos = Ensemble_run(locations[3])
   
    print('cropping to: {}\ncalculaing areal averages'.format(catchment),flush= True)
   
    # cropping for selected catchment
    rad_c = rad.extract_by_shp(load_catchment(catchment)).avg_areal_prec().aggr_temporal(3)
   
    fore_c = fore.eps_extract_by_shp(load_catchment(catchment)).avg_areal_prec()
    
    rad_cos_c= rad_cos.extract_by_shp(load_catchment(catchment)).avg_areal_prec().aggr_temporal(3) 

    for_cos_c = for_cos.eps_extract_by_shp(load_catchment(catchment)).avg_areal_prec()
    
    # garbage collection
    del rad ,fore, rad_cos, for_cos   
    
    library_icon = {'radar': rad_c,
               'ICOND2EPS': fore_c,
               }
    
    library_cosmo = {'radar': rad_cos_c, 
                    'COSMOD2EPS': for_cos_c }
    

    
    del rad_c,fore_c,rad_cos_c,for_cos_c
    
    return library_icon, library_cosmo      




def evaluate_for_catchment(cats):

    # TODO change
    max_lead = 24
    
    variables = [ 'ICOND2EPS q10%','ICOND2EPS q25%','ICOND2EPS q50%','ICOND2EPS q75%','ICOND2EPS q90%',
                'COSMOD2EPS q10%','COSMOD2EPS q25%','COSMOD2EPS q50%','COSMOD2EPS q75%','COSMOD2EPS q90%',]
    
    roc_auc_dic = {key: [None]*int(max_lead/3)  for key in variables}   
    
    
    c = 0 # leadtime index counter
    
    for lead in range(3,max_lead+1,3):
        print('Catchment: {}\ncalculating for lead time: {}hrs'.format(cats,lead),flush= True)
        
        # loading the files
        files = load_event_file(lead)
        # creating the instances for a specific catchment
        lib_ic, lib_cos = gen_for_catchment(cats,files)
        
        
        
        # Create a progress bar
        progress_bar = tqdm(total=10, desc='Calculating evaluation metrics.....')
        
        icon = ROC(lib_ic['radar'], lib_ic['ICOND2EPS'])
        cosm = ROC(lib_cos['radar'], lib_cos['COSMOD2EPS'])
 
    
        roc_auc_dic[variables[0]][c] = icon.roc_auc(10)
        progress_bar.update(1)   
        roc_auc_dic[variables[1]][c] = icon.roc_auc(25)
        progress_bar.update(1)   
        roc_auc_dic[variables[2]][c] = icon.roc_auc(50)
        progress_bar.update(1)   
        roc_auc_dic[variables[3]][c] = icon.roc_auc(75)
        progress_bar.update(1)    
        roc_auc_dic[variables[4]][c] = icon.roc_auc(90)
        progress_bar.update(1)    
        roc_auc_dic[variables[5]][c] = cosm.roc_auc(10)
        progress_bar.update(1)    
        roc_auc_dic[variables[6]][c] = cosm.roc_auc(25)
        progress_bar.update(1)        
        roc_auc_dic[variables[7]][c] = cosm.roc_auc(50)
        progress_bar.update(1)   
        roc_auc_dic[variables[8]][c] = cosm.roc_auc(75)
        progress_bar.update(1)    
        roc_auc_dic[variables[9]][c] = cosm.roc_auc(90)
        progress_bar.update(1)    
  
        
        # Close the progress bar
        progress_bar.close()
            
        del lib_ic, lib_cos
            
       
        c+=1
    
    return {'Area under ROC curve':roc_auc_dic}



def plot_curves_ICON(results, catchment_name):
    # PLOTTING
    # TODO  change
    max_lead = 24
    leads = list(range(3,max_lead+1,3))
    
    for res in list(results.keys()):
        # create the figure
        sns.set_theme(style="whitegrid")
        fig, ax = plt.subplots(dpi=750)
        
        # plotting the results
        ax.fill_between(leads, results[res]['ICOND2EPS q90%'], results[res]['ICOND2EPS q10%'], color= 'skyblue' , label='q90-10%')
        plt.plot(leads, results[res]['ICOND2EPS q50%'],'-', color = 'blue',label='q50%' )
        ax.fill_between(leads, results[res]['ICOND2EPS q75%'], results[res]['ICOND2EPS q25%'], color= 'steelblue', label='q75-25%' )
        
        
        

        ax.legend()
        plt.xticks(leads)
        ax.set_xlim(3,max_lead)
        ax.set_title("Performance of ICOND2EPS  \n{} catchment".format(catchment_name))
        ax.set_xlabel('Leadtime (hrs)')
        ax.set_ylabel(res)
        ax.set_ylim(0.74, 0.91)
        
        print('exporting results to file')
        # plt.savefig("//vs-grp07.zih.tu-dresden.de/howa/work/students/Mohamed_Elghorab/results/product analysis/{}/{}".format(catchment_name,res), bbox_inches = 'tight')        
    
    return plt.show()



def plot_curves_COSMO(results, catchment_name):
    # PLOTTING
    # TODO  change
    max_lead = 24
    leads = list(range(3,max_lead+1,3))
    
    for res in list(results.keys()):
        # create the figure
        sns.set_theme(style="whitegrid")
        fig, ax = plt.subplots(dpi=750)
        
        # plotting the results
        
        ax.fill_between(leads, results[res]['COSMOD2EPS q90%'], results[res]['COSMOD2EPS q10%'], color= 'bisque' , label='q90-10%')
        plt.plot(leads, results[res]['COSMOD2EPS q50%'],'-', color = 'chocolate',label='q50%' )
        ax.fill_between(leads, results[res]['COSMOD2EPS q75%'], results[res]['COSMOD2EPS q25%'], color= 'orange', label='q75-25%' )
        
        

        ax.legend()
        plt.xticks(leads)
        ax.set_xlim(3,max_lead)
        ax.set_title("Performance of COSMOD2EPS  \n{} catchment".format(catchment_name))
        ax.set_xlabel('Leadtime (hrs)')
        ax.set_ylabel(res)
        ax.set_ylim(0.74, 0.91)
        print('exporting results to file')
        # plt.savefig("//vs-grp07.zih.tu-dresden.de/howa/work/students/Mohamed_Elghorab/results/product analysis/{}/{}".format(catchment_name,res), bbox_inches = 'tight')        
    
    return plt.show()

###############################################################################

def evaluate_single(cname):
    
    results = evaluate_for_catchment(cname)
    
    plot_curves_ICON(results, cname)
    plot_curves_COSMO(results, cname)

###############################################################################





    
################################################################################################

if __name__ == '__main__':
################################################################################################


    tic = time.time()
    
    cats = ["Dohna", "Oelsnitz", 'Zittau', 'Adorf', 'Geising', 'Grossschoenau',
            'Bad Elster', 'Lauenstein', 'Niederoderwitz','Seifhennersdorf']

   
    for cat in cats:
        evaluate_single(cat)
        
    print((time.time()-tic)/3600) # operation time in hours
    
    
  
