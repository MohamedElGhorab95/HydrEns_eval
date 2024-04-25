# -*- coding: utf-8 -*-
"""
Created on Wed Jul 12 15:03:45 2023

@author: M Elghorab
"""
from engine.fr_entities_tools import *
import numpy as np
import xarray as xr
import xskillscore as xs
import seaborn as sns
from sklearn.metrics import mean_squared_error

class Full_Ens(object):
    
    def __init__(self,observation_object, forecast_object):
        
        if isinstance(forecast_object,xr.Dataset):
             self.forecast_array = forecast_object
        else:
            
            self.forecast_array = forecast_object.average
            
        if isinstance(observation_object,xr.DataArray):
             self.observation_array = observation_object
        else:
            if not isinstance(observation_object.average,int):
                self.observation_array = observation_object.average
            else:
                self.observation_array = observation_object.rtrn_arr()
        
        
        # selecting the same time window
        # selecting the same time window
        missing_times = set(self.observation_array.time.values) - set(self.forecast_array.time.values)
        
        obs_df = self.observation_array.to_dataframe()
        obs_df = obs_df[~obs_df.index.isin(missing_times)]

        # Convert DataFrame back to DataArray
        try:
            self.observation_array = obs_df.to_xarray()['rainfall radar observation | Radolan-RW']
        except:
            self.observation_array = obs_df.to_xarray()['runoff_discharge']
        
        # transform dataset into data array
        variable_names = list(self.forecast_array.data_vars)
        member_arrays = [self.forecast_array[var_name] for var_name in variable_names]
        self.forecast_array = xr.concat(member_arrays, dim='member')
        self.forecast_array = self.forecast_array.assign_coords(member=np.arange(20))

    def crps(self):
        
        return float(xs.crps_ensemble(self.observation_array, self.forecast_array, dim='time'))
    
    
    def brier(self, threshold):
        
        
        return float(xs.brier_score(self.observation_array > threshold,
               (self.forecast_array > threshold).mean('member'),
               dim="time"))
    
    def rmse(self):
       
      
       error = []
      
       for mem in range(20):
           error.append(float(xs.rmse(self.observation_array,self.forecast_array[mem] ) / np.mean(self.observation_array)) ) 
           
       return np.mean(error)  
       
    def disc_dia(self):
        
        dia = xs.discrimination(self.observation_array >0,(self.forecast_array>0).mean('member') )
        
        score = np.sum(np.abs(dia.values[0]-dia.values[1]))
        
        
        prob = dia.forecast_probability.values
        # create the figure
        sns.set_theme(style="whitegrid")
        fig, ax = plt.subplots(dpi = 750)
        ax.bar(prob-0.025, dia.values[0], 0.05, label='Observed')
        ax.bar(prob+0.025, dia.values[1], 0.05, label='Not observed')
        
        # Add some text for labels, title and custom x-axis tick labels, etc.
        ax.set_ylabel('Likelihood')
        ax.set_xlabel('Forecast probability')
        ax.set_xticks(prob)
        ax.legend(loc='upper left', ncols=2)
        ax.set_ylim(0, 0.9)
        
        plt.show()

        return np.around(score,3)
    
    
    
# rad = Observation("Z:/work/students/Mohamed_Elghorab/Mohamed Elghorab_MSc_Working_files/netCDFs/fertig/radRW_ICO.nc").gen_observation_field().extract_by_shp("C:/Users/user/Documents/GitHub/HydrEns_eval/shp/Mugliz/mugliz_cats.shp").aggr_temporal(3).avg_areal_prec()

# icond2eps = Ensemble_run("Z:/work/students/Mohamed_Elghorab/Mohamed Elghorab_MSc_Working_files//netCDFs/3/3hour_icond2eps.nc").eps_extract_by_shp("C:/Users/user/Documents/GitHub/HydrEns_eval/shp/Mugliz/mugliz_merged.shp").avg_areal_prec()    
    
# x = Full_Ens(rad, icond2eps)

# # crps = x.crps()
# # brier = x.brier(3)


# # er = x.rmse(10)
# a = x.disc_dia()

