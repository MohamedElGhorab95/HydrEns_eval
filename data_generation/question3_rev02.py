# -*- coding: utf-8 -*-
"""
Created on Sun Jul  2 10:09:02 2023

@author: M Elghorab
"""
import pandas as pd
import numpy as np
import time
import matplotlib.pyplot as plt

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






def load_event_file(leadtime, prod):
    source_file_path = "data_generation/nc_paths.txt"
    with open(source_file_path, 'r') as f:
        # Read all the lines in the file
        lines = f.readlines()
        # TODO change  to 0 to read from server
        netcdf = lines[1] 
        
        
     
        radar_file = netcdf+ 'fertig' + '/' + 'radRW_ICO.nc'
        if prod == "deter":
            deter_file = netcdf+ str(leadtime) + '/' + '{}hour_icond2.nc'.format(leadtime)
            return [radar_file, deter_file]
        else:
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



def gen_for_catchment(catchment, locations, product, q):
    
    print('generating instances..........',flush= True)
    # print(f"generating instances..........")
    # progress bar graphics
   
    
    # creating instances 
    rad = Observation(locations[0]).gen_observation_field().aggr_temporal(3)
    if product == 'deter':
        fore = Deterministic_run(locations[1]).gen_deterministic_field()
    else:
        fore = Ensemble_run(locations[1])
    
    
   
    print('cropping to: {}'.format(catchment),flush= True)
   
    # cropping for selected catchment
    rad_c = rad.extract_by_shp(load_catchment(catchment))
    if product == 'deter':
        fore_c = fore.extract_by_shp(load_catchment(catchment))
    else:
        fore_c = fore.eps_extract_by_shp(load_catchment(catchment))
    # garbage collection
    del rad,  fore

   
    # average areal precipitation
    print('calculaing areal averages',flush= True)
   
    rad_c = rad_c.avg_areal_prec()
  
    fore_c = fore_c.avg_areal_prec()
 
   
    
    
    # generating quantiles
    if product != 'deter':
        print('generating quantiles',flush= True)
        fore_c = fore_c.gen_quantiles(q)
    
    
    return rad_c, fore_c



def calc_metric(observation, forecast, metric, threshold = 4):
    print("calculating {}..........".format(metric),flush= True)
    
    
    if metric == "ROC_auc":
        return (float(ROC(observation, forecast).roc_auc()))
        # a = ROC(observation[0], observation[idx])
        # a.roc_auc()
        # a.plot_roc()
    elif metric == "success_ratio":
        return (CONT(observation, forecast, threshold).sr()*100)
    elif metric == "CSI":
        return (CONT(observation, forecast, threshold).csi()*100)
    elif metric == "FAR":
        return (CONT(observation, forecast, threshold).far()*100)
    elif metric == "PSS":
        return (CONT(observation, forecast, threshold).pss()*100)
    elif metric == "cont.hits":
        return (CONT(observation, forecast, threshold).hits())
    elif metric == "cont.misses":
        return (CONT(observation, forecast, threshold).misses())
    elif metric == "cont.fal":
        return (CONT(observation, forecast, threshold).false_alarms())
    elif metric == 'f1':
        return (CONT(observation, forecast, threshold).f1())
    elif metric == 'f2':
        return (CONT(observation, forecast, threshold).f2())





def evaluate_for_catchment(cats,product,q):

   
    # TODO change
    max_lead = 21
   
    results = []
    
    for lead in range(3,max_lead+1,3):
        print('Catchment: {}\ncalculating for lead time: {}hrs'.format(cats,lead),flush= True)
        
        # loading the files
        files = load_event_file(lead,product)
        # creating the instances for a specific catchment
        rad, fore = gen_for_catchment(cats,files,product,q)
        
        # TODO calculate metrics | change the metric if required
        results.append(calc_metric(rad, fore, "CSI")) 
       
    return results


def create_dataframe(product,q=90):
    
    areas = {"Weiße Elster":327.61,  'Niederoderwitz':29.08,
             'Mügliz':242.38, 'Mandau':279.29,'Lauenstein':75.5,
             'Bad Elster':47.7,'Adorf':170.36 }
    
    data = {key: [None,areas[key]]  for key in areas.keys()}
    # data = {key: [key,areas[key]]  for key in areas.keys()}
    
    for a in areas.keys():
   
            data[a][0]=evaluate_for_catchment(a, product,q)   
    
    
    data2 = {data[key][1]: data[key][0] for key in data.keys()}
    # TODO  change this
    leads = np.arange(3, 21+1, 3)
    df = pd.DataFrame.from_dict(data2, orient='index')
    df.columns = leads
    df_sorted = df.sort_index()
    return df_sorted



def plot_contour(dataframe):
    # Convert the DataFrame to a 2D numpy array
    x = dataframe
    data = x.values

    # Create meshgrid for x and y values
    x, y = np.meshgrid(dataframe.columns, dataframe.index)
    
    fig, ax = plt.subplots(dpi=750)
    # Create the contour plot
    
    plt.contourf(x, y, data)
    colorbar = plt.colorbar()
    colorbar.set_label('Critical success index\nEnsemble 90% Quantile')
    plt.xticks(dataframe.columns)
    ax.set_xlabel('Leadtime hrs')
    ax.set_ylabel('Catchment area km2')
    ax.set_title('Skill score of ICOND2EPS forecast performance\nSpatial extent relationship')
    ax.grid(True, linestyle='dashed', alpha=0.3)
    # Save the figure with the specified DPI
    plt.savefig('contour_plot_csi75.png')
    
    return plt.show()



tic = time.time()
df = create_dataframe('eps',90)  
print("Operation time : {} hours".format((time.time()-tic)/3600), flush=True)   
plot_contour(df)
