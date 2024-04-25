# -*- coding: utf-8 -*-
"""
Created on Mon Mar 25 10:16:51 2024

@author: Administrator
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
        
        
        rad_cos = netcdf + 'fertig' + '/' + 'radRW_COS.nc'

        
    files = [ rad_cos ]
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
    
     
    
     
     # instantiating the cosmo fields
     rad_cos = Observation(locations[0]).gen_observation_field()
     

    
     print('cropping to: {}\ncalculaing areal averages'.format(catchment),flush= True)
    
     # cropping for selected catchment

     
     # rad_cos_c= rad_cos.extract_by_shp(load_catchment(catchment)).avg_areal_prec().aggr_temporal(3) 

     
     # garbage collection
     # del  rad_cos
  
     
     # library_cosmo = {'radar': rad_cos_c, 
     #                  'COSMOD2':det_cos_c,
     #                 'COSMOD2EPS': for_cos_c }
     

     
     
     
     return   rad_cos.extract_by_shp(load_catchment(catchment)).avg_areal_prec().aggr_temporal(3).average
 
    
 
   
# cats = ['Dohna', "Oelsnitz", 'Zittau', 'Adorf', 'Geising', 'Grossschoenau',
#         'Bad Elster', 'Lauenstein', 'Niederoderwitz','Seifhennersdorf']


cats = [ "Oelsnitz"]

data = []

for cat in cats:
       
    data.append(gen_for_catchment(cat,load_event_file(1)).values)


new = list(np.concatenate(data))





plt.plot([i for i in range(len(new))],new)

np.nanquantile(new, 0.987116)