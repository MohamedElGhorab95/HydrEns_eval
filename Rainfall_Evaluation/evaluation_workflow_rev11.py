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
import time
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
        
    files = [radar_file, deter_file, ensem_file]
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
   
   
    
    library = {'radar': rad_c,
               'deterministic_run': det_c,
               'ensemble_run': eps_c,
               }
   
    del rad_c, det_c, eps_c
    
    return library










def evaluate_for_catchment(cats):

    variables = {
               'main_run':['deterministic_run',None],
               'ensemble_q10':['ensemble_run',10],
               'ensemble_q25':['ensemble_run',25],
               'ensemble_q75':['ensemble_run',75],
               'ensemble_q90':['ensemble_run',90],
               'ensemble_mean':['ensemble_run','mean']} 
    # TODO change
    max_lead = 24
    threshold = 2.5
    
    # initialize the result dictionaries
    roc_auc_dic = {key: [None]*int(max_lead/3)  for key in variables}   
    roc_tot = []
    
    
    
    
    c = 0 # leadtime index counter
    
    for lead in range(3,max_lead+1,3):
        print('Catchment: {}\ncalculating for lead time: {}hrs'.format(cats,lead),flush= True)
        
        # loading the files
        files = load_event_file(lead)
        # creating the instances for a specific catchment
        library = gen_for_catchment(cats,files)
        
        
        
        # Create a progress bar
        progress_bar = tqdm(total=6, desc='Calculating evaluation metrics.................')
        roc_tot.append(ROC_totens(library['radar'], library['ensemble_run'], 2.5).area)
        for x in variables.keys():
            roc_auc_dic[x][c] = ROC(library['radar'], library[variables[x][0]],variables[x][1]).roc_auc(threshold)
   
            
            progress_bar.update(1)
            
            
        # Close the progress bar
        progress_bar.close()
            
        del library
        
        
        
        c+=1
    roc_auc_dic['Full Ensemble'] = roc_tot
    
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
        plt.plot(leads, results[res]['Full Ensemble'],label='Full Ensemble' )
        plt.plot(leads, results[res]['ensemble_q10'],ls='--',label='Q10%' )
        plt.plot(leads, results[res]['ensemble_q90'],ls='--',label='Q90%' )
        plt.plot(leads, results[res]['ensemble_q25'],ls='-.',label='Q25%' )
        plt.plot(leads, results[res]['ensemble_q75'],ls='-.',label='Q75%' )
        
        ax.legend()
        plt.xticks(leads)
        ax.set_xlim(3,max_lead)
        
        if res   == 'Area under ROC curve': ax.set_ylim(0.5,1.0)
        
        
        ax.set_title("Performance of ICOND2 forecast  \n{} catchment".format(catchment_name))
        ax.set_xlabel('Leadtime (hrs)')
        ax.set_ylabel(res)
        
        print('exporting results to file')
        plt.savefig("//vs-grp07.zih.tu-dresden.de/howa/work/students/Mohamed_Elghorab/results/{}/{}".format(catchment_name,res), bbox_inches = 'tight')        
    
    return plt.show()
    


###############################################################################

def evaluate_single(cname):
    
    results = evaluate_for_catchment(cname)
    plot_curves(results, cname)




tic = time.time()
cats = ['Dohna', "Oelsnitz", 'Zittau', 'Adorf', 'Geising', 'Grossschoenau',
        'Bad Elster', 'Lauenstein', 'Niederoderwitz','Seifhennersdorf']

for cat in cats:
    evaluate_single(cat)
    
print((time.time()-tic)/3600) # operation time in hours    
   