# -*- coding: utf-8 -*-
"""
Created on Mon Mar 25 13:43:30 2024

@author: Administrator
"""

import xarray as xr
import xskillscore as xs
import numpy as np
from sklearn import metrics
from tsmoothie.smoother import *
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
from engine.fr_Conttools import *
from engine.fr_entities_tools import *
from sklearn.metrics import accuracy_score, confusion_matrix, roc_auc_score, roc_curve


class ROC(object):
    def __init__(self, observation_object, forecast_object,quantile = None):
        
        
        
        
        # self.q = quantile
        
        
        
        # if self.q != None:
            
        #     forecast_object = forecast_object.gen_quantiles(self.q).arr
        
        
        
        if isinstance(forecast_object, R_Forecast) == False:
            if isinstance(forecast_object,xr.DataArray):
                 self.forecast_array = forecast_object
            else:
                
                self.forecast_array = forecast_object.average
        else:
            
            self.forecast_array = forecast_object.fr
            
        if isinstance(observation_object, R_Observation) == False:    
            if isinstance(observation_object,xr.DataArray):
                 self.observation_array = observation_object
            else:
                if not isinstance(observation_object.average,int):
                    self.observation_array = observation_object.average
                else:
                    self.observation_array = observation_object.rtrn_arr()
        else:
            
            self.observation_array = observation_object.fr
       
        self.pod  = [] # probability of detection
        self.pofd = [] # probability of false detection
        
        
        
        
        # selecting the same time window
        
        missing_times = set(self.observation_array.time.values) - set(self.forecast_array.time.values)
        
        obs_df = self.observation_array.to_dataframe()
        obs_df = obs_df[~obs_df.index.isin(missing_times)]
        
        
        # Convert DataFrame back to DataArray
        try:
            self.observation_array = obs_df.to_xarray()['rainfall radar observation | Radolan-RW']
        except:
            self.observation_array = obs_df.to_xarray()['runoff_discharge']
            
            
    def roc_auc(self, quantile= None):
        '''
        This method returns the area under the ROC curve, the default setting is to 
        return the area for the whole extent of the data array, however,
        if calculation_extents = "cell_wise" the method returns a data array with the area
        for each pixel
   
        Parameters
        ----------
        calculation_extents : str, optional
            could be "cell_wise". The default is "global".
   
        Returns
        -------
        auc : data array
            area under the curve for the ROC.
   
        '''
        
        # get the upper and lower bounds of the rainfall fields
        upper = self.observation_array.max(skipna=True).values
        lower = self.observation_array.min(skipna=True).values
        bins = 19
        bin_width = (upper-lower)/bins
        
        thresholds = np.arange(lower, upper + bin_width, bin_width)
        

        # check if the data is determenistic or Ensemble

        try:       
            members = list(self.forecast_array.variables)[3:]
            
            # create the dataframes to store all curves 
            pod_df = pd.DataFrame(columns=members)
            pofd_df = pd.DataFrame(columns=members)
            
            
            
            
            for mem in members:
              
    
                roc = xs.roc(self.observation_array, self.forecast_array[mem], thresholds, dim=['time'],return_results='all_as_tuple')
                
                pod = roc[1]
                pofd = roc[0]
                
                pod_df[mem]  = pod
                pofd_df[mem] = pofd
                
               
            
            # generate quantiles 
            pod_df= pod_df.quantile(quantile/100,axis=1)
            pofd_df= pofd_df.quantile(quantile/100,axis=1)   
            self.ar = np.around(metrics.auc(pofd_df.sort_values(),pod_df.sort_values()),3)
        except:
            # deterministic run
            roc = xs.roc(self.observation_array, self.forecast_array, thresholds, dim=['time'],return_results='all_as_tuple')
            pod_df = roc[1]
            pofd_df = roc[0]
            self.ar = roc[2].values
 
        
        sns.set_theme(style="whitegrid")
        fig, ax = plt.subplots(dpi=750)
        
        
        ax.plot(pofd_df, pod_df)
        # Add a 45-degree line
        ax.plot([0, 1], [0, 1], transform=ax.transAxes, ls='--', color='k')
            
        ax.set_xlabel('PROBABILITY OF FALSE DETECTION')
        ax.set_ylabel('PROBABILITY OF DETECTION')
        
        ax.set_xlim(0,1)
        ax.set_ylim(0,1)
        
        plt.show()
    

        
       
        
        # self.ar = np.around(metrics.auc(pofd_df,pod_df),3)
        # self.ar = np.around(np.trapz(pofd_df,pod_df),3)
        # return auc
        return self.ar
        
    
    
    
# ========================================================================================================================================================


# ==============================================================================
###########################   Testing | Examples  #############################
# ==============================================================================


if __name__ == '__main__':
    from engine.fr_entities_tools import *
    
    
    rad = Observation("E:/Data2/netCDFs/fertig/radRW_ICO.nc").gen_observation_field().extract_by_shp("shp/Mugliz/mugliz_cats.shp").aggr_temporal(3).avg_areal_prec()
    # icond2 = Deterministic_run("E:/Data2/netCDFs/3/3hour_icond2.nc").gen_deterministic_field().extract_by_shp("shp/Mugliz/mugliz_cats.shp").avg_areal_prec()
    icond2eps = Ensemble_run("E:/Data2/netCDFs/3/3hour_icond2eps.nc").eps_extract_by_shp("shp/Mugliz/mugliz_cats.shp").avg_areal_prec()
    
    b = ROC(rad,icond2eps)
    b.roc_auc(quantile=10)
    # c = ROC(rad,icond2)
    # c.roc_auc()   
    
    
    fpr, tpr, thresholds = roc_curve(b.observation_array.values>0, b.forecast_array[members[0]].values)
