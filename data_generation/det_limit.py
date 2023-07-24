# -*- coding: utf-8 -*-
"""
Created on Thu Jul 13 18:02:46 2023

@author: M Elghorab
"""

import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import root_scalar
from scipy.interpolate import interp1d
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
    # eps = Ensemble_run(locations[2])
   
    print('cropping to: {}'.format(catchment),flush= True)
    # print(f"cropping to: {catchment}")
    # cropping for selected catchment
    rad_c = rad.extract_by_shp(load_catchment(catchment))
    det_c = det.extract_by_shp(load_catchment(catchment))
    # eps_c = eps.eps_extract_by_shp(load_catchment(catchment))
    # garbage collection
    # del rad, det, eps
   
    # average areal precipitation
    print('calculaing areal averages',flush= True)
    # print(f"calculaing areal averages")
    rad_c = rad_c.avg_areal_prec()
    det_c = det_c.avg_areal_prec()
    # eps_c = eps_c.avg_areal_prec()
   
   
    
    # library = {'radar': rad_c,
               # 'deterministic_run': det_c,
               # 'ensemble_run': eps_c,
               # }
               
    library = {'radar': rad_c,
               'deterministic_run': det_c,
               
               }
   
    # del rad_c, det_c, eps_c
    
    return library

def get_limit(x,y,z):
    # Interpolate the values of y and z
    interp_y = interp1d(x, y)
    interp_z = interp1d(x, z)
    
    # Create a function that calculates the difference between y and z at a given x-coordinate
    diff_func = lambda x_val: interp_y(x_val) - interp_z(x_val)
    try:
        intersection_result = root_scalar(diff_func, bracket=[x[0], x[-1]], method='brentq')
        return int(intersection_result.root)
    except:
        return 0
    
    
def evaluate_for_catchment(cats, threshold):

    hits   = []
    misfal = []
    max_lead = 9
    
    
    for lead in range(3,max_lead+1,3):
        print('Catchment: {}\ncalculating for lead time: {}hrs'.format(cats,lead),flush= True)
        
        # loading the files
        files = load_event_file(lead)
        # creating the instances for a specific catchment
        library = gen_for_catchment(cats,files)
        
        hits.append(CONT(library['radar'], library['deterministic_run'],threshold).hits(None))
  
        misfal.append(CONT(library['radar'], library['deterministic_run'],threshold).misses(None) + CONT(library['radar'], library['deterministic_run'],threshold).false_alarms(None))
       
      
    return np.array(hits), np.array(misfal)


limit = []
lt = np.arange(3,10,3)
for t in range(1,11):

    h, m = evaluate_for_catchment('Mandau', t)
    
    plt.plot(lt,h)
    plt.plot(lt,m)
    
    limit.append(get_limit(lt, h, m))

plt.plot(limit, range(1,11))
plt.xlim(0,np.max(limit))