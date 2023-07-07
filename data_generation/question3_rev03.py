# -*- coding: utf-8 -*-
"""
Created on Sun Jul  2 10:09:02 2023

@author: M Elghorab
"""

import numpy as np
import time
import matplotlib.pyplot as plt
import seaborn as sns
# leads = np.arange(3, 19, 3)

# areas = {"Weiße Elster":327.61,  'Niederoderwitz':29.08,
#           'Mügliz':242.38, 'Mandau':279.29,'Lauenstein':75.5,
#           'Bad Elster':47.7,'Adorf':170.36 }


# area = [327.61,29.08,242.38,279.29,75.5,47.7,170.36]
# area.sort()



# valzo = dict(zip(area, [[random.random()*50 for i in range(6)] for x in range(8)]))

# df = pd.DataFrame.from_dict(valzo, orient='index')

# df.columns = leads




################################################################################
################################################################################
################################################################################


from HydrEns_eval.engine.fr_entities_tools import *
from HydrEns_eval.engine.fr_Conttools import *
from HydrEns_eval.engine.fr_ROCtools import *






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
        if name == 'Mandau':
            return lines[2].strip()
        if name == 'Mügliz':
            return lines[4].strip()
        if name == 'Weiße Elster':
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
 
   
    
    
    # generating quantiles
    print('generating quantiles',flush= True)
    fore_c = fore_c.gen_quantiles(50)
    
    
    return rad_c, fore_c



def calc_metric(observation, forecast):
    print("calculating ROC..........",flush= True)
    
    
 
    return (float(ROC(observation, forecast).roc_auc()))
        # a = ROC(observation[0], observation[idx])
        # a.roc_auc()
        # a.plot_roc()
    





def evaluate_for_catchment(cats):

   
    # TODO change
    max_lead = 21
   
    results = []
    
    for lead in range(3,max_lead+1,3):
        print('Catchment: {}\ncalculating for lead time: {}hrs'.format(cats,lead),flush= True)
        
        # loading the files
        files = load_event_file(lead)
        # creating the instances for a specific catchment
        rad, fore = gen_for_catchment(cats,files)
        
        # TODO calculate metrics | change the metric if required
        results.append(calc_metric(rad, fore)) 
       
    return results


def create_for_regions():
    
    region1 = {"Weiße Elster":328, 'Bad Elster':48 ,'Adorf':170 }
    
    region2 = {'Mügliz':242, 'Lauenstein':76 }
    
    region3 = {'Niederoderwitz':29, 'Mandau':279, 'Seifhennersdorf':75}
    
    
    data1   = {key: [None,region1[key]]  for key in region1.keys()}
    for a in data1.keys():
            data1[a][0]=evaluate_for_catchment(a)   
    data1 = {data1[key][1]: data1[key][0] for key in data1.keys()}
    
    data2   = {key: [None,region2[key]]  for key in region2.keys()}
    for a in data2.keys():
            data2[a][0]=evaluate_for_catchment(a)
    data2 = {data2[key][1]: data2[key][0] for key in data2.keys()}
    
    data3   = {key: [None,region3[key]]  for key in region3.keys()}
    for a in data3.keys():
            data3[a][0]=evaluate_for_catchment(a)
    data3 = {data3[key][1]: data3[key][0] for key in data3.keys()}
    
    return [data1, data2, data3]
    
        
def plot_spatial(list_of_dics):
    
    max_lead = 21
    leads = list(range(3,max_lead+1,3))
    
    region = 0
    regions = ['Western Saxony', 'Central Saxony', 'Eastern Saxony']
    for l in list_of_dics:
        sns.set_theme(style="whitegrid")
        fig, ax = plt.subplots(dpi=750)
        
        l = dict(sorted(l.items()))
        
        for u in l.keys():
            
           
            
            plt.plot(leads, l[u],label='{}km\u00b2'.format(u) )
            
        ax.legend()
        plt.xticks(leads)
        ax.set_title("Skill score of ICOND2EPS forecast performance\n{}".format(regions[region]))
        ax.set_xlabel('Leadtime (hrs)')
        ax.set_ylabel('Area under ROC curve\nEnsemble median')
        plt.xlim(3, max_lead)
        region +=1
    return plt.show()





tic = time.time()
df = create_for_regions()
print("Operation time : {} hours".format((time.time()-tic)/3600), flush=True)   
plot_spatial(df)
