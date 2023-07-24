# -*- coding: utf-8 -*-
"""
Created on Sun Jul  2 10:09:02 2023

@author: M Elghorab
"""

import time
import matplotlib.pyplot as plt
import seaborn as sns
from HydrEns_eval.engine.fr_entities_tools import *
from HydrEns_eval.engine.fr_Conttools import *
from HydrEns_eval.engine.fr_ROCtools import *

################################################################################
################################################################################
################################################################################




def load_event_file(leadtime):
    source_file_path = "data_generation/nc_paths.txt"
    with open(source_file_path, 'r') as f:
        # Read all the lines in the file
        lines = f.readlines()
        # TODO change  to 0 to read from server
        netcdf = lines[1] 
        
        
     
        radar_file = netcdf+ 'fertig' + '/' + 'radRW_ICO.nc'
        
        ensem_file = netcdf+ str(leadtime) + '/' + '{}hour_icond2eps.nc'.format(leadtime)
        return [radar_file, ensem_file]
    
   


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

    fore = Ensemble_run(locations[1])
    
    
   
    print('cropping to: {}'.format(catchment),flush= True)
   
    # cropping for selected catchment
    rad_c = rad.extract_by_shp(load_catchment(catchment))
    
    fore_c = fore.eps_extract_by_shp(load_catchment(catchment))
    # garbage collection
    del rad,  fore

   
    # average areal precipitation
    print('calculaing areal averages',flush= True)
   
    rad_c = rad_c.avg_areal_prec()
  
    fore_c = fore_c.avg_areal_prec()
     
    
    return rad_c, fore_c



def calc_metric(observation, forecast , metric):
    print("calculating performance metric..........",flush= True)
    
    
    if metric == 'Frequency bias':
        return CONT(observation, forecast, 3).fbias(50)
    elif metric == 'Accuracy':
        return CONT(observation, forecast, 3).acc(50)
    elif metric == "Pierce's skill score":
        return CONT(observation, forecast, 3).pss(50)
    elif metric == 'Area under ROC curve':
        return (float(ROC(observation, forecast).roc_auc(quantile=50)))
        





def evaluate_for_catchment(cats, met):

   
    # TODO change
    max_lead = 24
   
    results = []
    
    for lead in range(3,max_lead+1,3):
        print('Catchment: {}\ncalculating for lead time: {}hrs'.format(cats,lead),flush= True)
        
        # loading the files
        files = load_event_file(lead)
        # creating the instances for a specific catchment
        rad, fore = gen_for_catchment(cats,files)
        
        # TODO calculate metrics | change the metric if required
        results.append(calc_metric(rad, fore, met)) 
       
    return results


def create_for_regions(metric):
    
    region1 = {"Oelsnitz":328, 'Bad Elster':48 ,'Adorf':170 }
    
    region2 = {'Dohna':200, 'Lauenstein':76, 'Geising':26 }
    
    region3 = {'Niederoderwitz':29, 'Zittau':279, 'Seifhennersdorf':75 , 'Grossschoenau': 162}
    
    
    data1   = {key: [None,region1[key]]  for key in region1.keys()}
    for a in data1.keys():
            data1[a][0]=evaluate_for_catchment(a, metric)   
    data1 = {data1[key][1]: [data1[key][0],key] for key in data1.keys()}
    
    data2   = {key: [None,region2[key]]  for key in region2.keys()}
    for a in data2.keys():
            data2[a][0]=evaluate_for_catchment(a, metric)
    data2 = {data2[key][1]: [data2[key][0],key] for key in data2.keys()}
    
    data3   = {key: [None,region3[key]]  for key in region3.keys()}
    for a in data3.keys():
            data3[a][0]=evaluate_for_catchment(a, metric)
    data3 = {data3[key][1]: [data3[key][0],key] for key in data3.keys()}
    
    return [data1, data2, data3]
    
        
def plot_spatial(list_of_dics, metric):
    
    max_lead = 24
    leads = list(range(3,max_lead+1,3))
    
    region = 0
    regions = ['Weiße Elster region', 'Müglitz region', 'Eastern Saxony region']
    for l in list_of_dics:
        sns.set_theme(style="whitegrid")
        fig, ax = plt.subplots(dpi=1000)
        
        l = dict(sorted(l.items()))
        
        for u in l.keys():
            
           
            
            plt.plot(leads, l[u],label='{}km\u00b2 |{}'.format(u,l[u][1]) )
        
        # TODO change metric name in y axis and folder    
        ax.legend()
        plt.xticks(leads)
        ax.set_title("Performance of ICOND2EPS forecast \n{}".format(regions[region]))
        ax.set_xlabel('Leadtime (hrs)')
        ax.set_ylabel('{}\nEnsemble median'.format(metric))
        plt.xlim(3, max_lead)
        plt.savefig("//vs-grp07.zih.tu-dresden.de/howa/work/students/Mohamed_Elghorab/results/spatial extent analysis/{}/{}".format(metric, regions[region]), bbox_inches = 'tight')        
        region +=1
        

    return plt.show()



mets = ['Frequency bias' , 'Accuracy', "Pierce's skill score", 'Area under ROC curve'  ]

for m in mets:

    tic = time.time()
    df = create_for_regions(m)
    print("Operation time : {} hours".format((time.time()-tic)/3600), flush=True)   
    plot_spatial(df,m)
