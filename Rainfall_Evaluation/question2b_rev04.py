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
    source_file_path = "in/nc_paths.txt"
    with open(source_file_path, 'r') as f:
        # Read all the lines in the file
        lines = f.readlines()
        
        netcdf = lines[0] 
        
        
     
        radar_file = netcdf+ 'fertig' + '/' + 'radRW_ICO.nc'
        deter_file = netcdf+ str(leadtime) + '/' + '{}hour_icond2.nc'.format(leadtime)
        ensem_file = netcdf+ str(leadtime) + '/' + '{}hour_icond2eps.nc'.format(leadtime)
        
        
        rad_cos = netcdf + 'fertig' + '/' + 'radRW_COS.nc'
        det_cos = netcdf+ str(leadtime) + '/' + '{}hour_cosmod2.nc'.format(leadtime)
        ens_cos = netcdf+ str(leadtime) + '/' + '{}hour_cosmod2eps.nc'.format(leadtime)
        
    files = [radar_file, deter_file, ensem_file, rad_cos, det_cos, ens_cos]
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
    fore = Ensemble_run(locations[2])
    
    # instantiating the cosmo fields
    rad_cos = Observation(locations[3]).gen_observation_field()
    
    det_cos = Deterministic_run(locations[4]).gen_deterministic_field()
    for_cos = Ensemble_run(locations[5])
   
    print('cropping to: {}\ncalculaing areal averages'.format(catchment),flush= True)
   
    # cropping for selected catchment
    rad_c = rad.extract_by_shp(load_catchment(catchment)).avg_areal_prec().aggr_temporal(3)
    det_c = det.extract_by_shp(load_catchment(catchment)).avg_areal_prec()
    fore_c = fore.eps_extract_by_shp(load_catchment(catchment)).avg_areal_prec()
    
    rad_cos_c= rad_cos.extract_by_shp(load_catchment(catchment)).avg_areal_prec().aggr_temporal(3) 
    det_cos_c = det_cos.extract_by_shp(load_catchment(catchment)).avg_areal_prec()
    for_cos_c = for_cos.eps_extract_by_shp(load_catchment(catchment)).avg_areal_prec()
    
    # garbage collection
    del rad, det ,fore, rad_cos, det_cos, for_cos   
    
    library_icon = {'radar': rad_c,
               'ICOND2': det_c,
               'ICOND2EPS': fore_c,
               }
    
    library_cosmo = {'radar': rad_cos_c, 
                     'COSMOD2':det_cos_c,
                    'COSMOD2EPS': for_cos_c }
    

    
    
    
    return library_icon, library_cosmo      




def evaluate_for_catchment(cats):

   
    # TODO change
    max_lead = 24
    thr_ICO = 2.5
    thr_COS = 2
    variables = ['ICOND2', 'ICOND2EPS', 'COSMOD2', 'COSMOD2EPS']
    
    roc_auc_dic = {key: [None]*int(max_lead/3)  for key in variables}   
    
    
    c = 0 # leadtime index counter
    
    for lead in range(3,max_lead+1,3):
        print('Catchment: {}\ncalculating for lead time: {}hrs'.format(cats,lead),flush= True)
        
        # loading the files
        files = load_event_file(lead)
        # creating the instances for a specific catchment
        lib_ic, lib_cos = gen_for_catchment(cats,files)
        
        
        
        # Create a progress bar
        progress_bar = tqdm(total=4, desc='Calculating evaluation metrics.....')
        
        # for x in variables[:2]:
        #     roc_auc_dic[x][c] = ROC_totens(lib_ic['radar'], lib_ic[x],threshold).area
   
            
        #     progress_bar.update(1)
        
        roc_auc_dic[variables[0]][c] = ROC(lib_ic['radar'], lib_ic[variables[0]]).roc_auc(thr_ICO)
        progress_bar.update(1)   
        roc_auc_dic[variables[1]][c] = ROC_totens(lib_ic['radar'], lib_ic[variables[1]],thr_ICO).area
        progress_bar.update(1)   
        roc_auc_dic[variables[2]][c] = ROC(lib_cos['radar'], lib_cos[variables[2]]).roc_auc(thr_COS)
        progress_bar.update(1)   
        roc_auc_dic[variables[-1]][c] = ROC_totens(lib_cos['radar'], lib_cos[variables[-1]],thr_COS).area
        progress_bar.update(1)    
            
            
            
        # Close the progress bar
        progress_bar.close()
            
        del lib_ic, lib_cos
            
       
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
        plt.plot(leads, results[res]['ICOND2'],'-', color = 'blue',label='Main run | ICOND2' )
        plt.plot(leads, results[res]['ICOND2EPS'],'--',color = 'blue',label='Ensemble | ICOND2' )
        
        plt.plot(leads, results[res]['COSMOD2'],'-',color = 'orange',label='Main run | COSMOD2' )
        plt.plot(leads, results[res]['COSMOD2EPS'],'--',color = 'orange',label='Ensemble | COSMOD2' )
        
        

        ax.legend()
        plt.xticks(leads)
        ax.set_xlim(3,max_lead)
        ax.set_title("Performance of different forecast products  \n{} catchment".format(catchment_name))
        ax.set_xlabel('Leadtime (hrs)')
        ax.set_ylabel(res)
        
        print('exporting results to file')
        # plt.savefig("//vs-grp07.zih.tu-dresden.de/howa/work/students/Mohamed_Elghorab/results/product analysis/{}/{}".format(catchment_name,res), bbox_inches = 'tight')        
    
    return plt.show()


###############################################################################

def evaluate_single(cname):
    tic = time.time()
    results = evaluate_for_catchment(cname)
    print((time.time()-tic)/3600) # operation time in hours
    plot_curves(results, cname)


# cats = ['Dohna', "Oelsnitz", 'Zittau', 'Adorf', 'Geising', 'Grossschoenau',
#         'Bad Elster', 'Lauenstein', 'Niederoderwitz','Seifhennersdorf']

cats = [ 'Oelsnitz']

for cat in cats:
    evaluate_single(cat)