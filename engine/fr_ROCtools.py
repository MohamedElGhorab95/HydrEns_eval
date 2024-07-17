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
from sklearn.metrics import roc_curve


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
            
        
        
        self.calc_roc()
        #################################################################################################    
    
    
    def calc_roc(self):
        # calculate ROC    
        # get the upper and lower bounds of the rainfall fields
        upper = self.observation_array.max(skipna=True).values
        lower = self.observation_array.min(skipna=True).values
        bins = 1999
        bin_width = (upper-lower)/bins
        
        thresholds = np.arange(lower, upper , bin_width)[1:]
        
       
        # check if the data is determenistic or Ensemble
    
        try:       
            members = list(self.forecast_array.variables)[3:]
            
            # create the dataframes to store all curves 
            self.pod_df = pd.DataFrame(columns=members)
            self.pofd_df = pd.DataFrame(columns=members)
            self.area_df = []
            
            
            # loop over the members 
            for mem in members:
                
                pod_kl  = [0]
                pofd_kl = [0]
                
                # calculate pod and pofd for each threshold
                for thr in thresholds:
                    fpr, tpr, t = roc_curve(self.observation_array.values>thr, self.forecast_array[mem].values>thr)
                    pod_kl.append(tpr[1])
                    pofd_kl.append(fpr[1])
       
    
                pod_kl.append(1)
                pofd_kl.append(1)
                
                # Sort list A ascendingly and get the sorted indices
                sorted_indices = sorted(range(len(pofd_kl)), key=lambda i: pofd_kl[i])
                
                # Sort list A ascendingly using the sorted indices
                pofd_kl = [pofd_kl[i] for i in sorted_indices]
                
                # Rearrange list B based on the sorted indices of list A
                pod_kl = [pod_kl[i] for i in sorted_indices]

                
                # add all points to the members dataframe 
                self.pod_df[mem]  = pod_kl
                self.pofd_df[mem] = pofd_kl
                
              
                
                self.area_df.append(np.around(metrics.auc(self.pofd_df[mem],self.pod_df[mem]),3))
               
            del self.observation_array, self.forecast_array
          
        except:
            # deterministic run
            
            pod_kl  = [0]
            pofd_kl = [0]
            
            for thr in thresholds:
                fpr, tpr, t = roc_curve(self.observation_array.values>thr, self.forecast_array.values>thr)
                pod_kl.append(tpr[1])
                pofd_kl.append(fpr[1])
       
            
            del self.observation_array, self.forecast_array
            pod_kl.append(1)
            pofd_kl.append(1)
            
            # create the dataframes to store all curves 
            self.pod_df = pd.DataFrame(pod_kl)
            self.pofd_df = pd.DataFrame(pofd_kl)
            self.area_df = pd.DataFrame(np.around(metrics.auc(self.pofd_df,self.pod_df),3))
           
            
            
            
            
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
        
      
        # check if the data is determenistic or Ensemble

        try:       
            
            # generate quantiles 
            # self.pod_df_q= self.pod_df.quantile(quantile/100,axis=1)
            # self.pofd_df_q= self.pofd_df.quantile(quantile/100,axis=1)   
            # # sorting the data
            # self.pofd_df_q = self.pofd_df_q.sort_values()
            # self.pod_df_q = self.pod_df_q.reindex(self.pofd_df_q.index)
            # # computing area under the curve
            # self.ar = np.around(metrics.auc(self.pofd_df_q,self.pod_df_q),3)

            self.ar = np.quantile(self.area_df, quantile/100)
            
            
            
        except:
                 
        
            # sorting the data
            self.pofd_df = self.pofd_df.sort_values(by=0)
            self.pod_df = self.pod_df.reindex(self.pofd_df.index)
            # computing area under the curve
            self.ar = np.around(metrics.auc(self.pofd_df,self.pod_df),3)
     
        
     
        
        # sns.set_theme(style="whitegrid")
        # fig, ax = plt.subplots(dpi=750)
        
        
        
        # try:
        #     ax.plot(self.pofd_df_q, self.pod_df_q)
        # except:
        #     ax.plot(self.pofd_df, self.pod_df)
        # # Add a 45-degree line
        # ax.plot([0, 1], [0, 1], transform=ax.transAxes, ls='--', color='k')
            
        # ax.set_xlabel('PROBABILITY OF FALSE DETECTION')
        # ax.set_ylabel('PROBABILITY OF DETECTION')
        
        # ax.set_xlim(0,1)
        # ax.set_ylim(0,1)
        
        # plt.show()
        

        
       
        
      
        return self.ar
        
    
    
    
# ========================================================================================================================================================


# ==============================================================================
###########################   Testing | Examples  #############################
# ==============================================================================


if __name__ == '__main__':
    from engine.fr_entities_tools import *
    
    
    # rad = Observation("E:/Data2/netCDFs/fertig/radRW_COS.nc").gen_observation_field().extract_by_shp("shp/Mugliz/mugliz_cats.shp").aggr_temporal(3).avg_areal_prec()
    # icond2 = Deterministic_run("E:/Data2/netCDFs/3/3hour_icond2.nc").gen_deterministic_field().extract_by_shp("shp/Mugliz/mugliz_cats.shp").avg_areal_prec()
    # icond2eps = Ensemble_run("E:/Data2/netCDFs/15/15hour_icond2eps.nc").eps_extract_by_shp("shp/Mugliz/mugliz_cats.shp").avg_areal_prec()
    # cosmod2eps = Ensemble_run("E:/Data2/netCDFs/24/24hour_cosmod2eps.nc").eps_extract_by_shp("shp/Mugliz/mugliz_cats.shp").avg_areal_prec()
    
    # b = ROC(rad,icond2eps)    
    # b.roc_auc(quantile=10)
    # c = ROC(rad,icond2)
    # c.roc_auc()  
    d = ROC(rad,cosmod2eps)
    # d.roc_auc(quantile=10)
    
    