
# -*- coding: utf-8 -*-
"""
Created on Wed Apr 26 16:23:18 2023

@author: M Elghorab
"""
import xarray as xr
import xskillscore as xs
import numpy as np
from sklearn.metrics import auc
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
from engine.fr_Conttools import *
from engine.fr_entities_tools import *

class ROC(object):
    def __init__(self, observation_object, forecast_object,quantile = None):
        
        
        
        
        self.q = quantile
        
        
        
        if self.q != None:
            
            forecast_object = forecast_object.gen_quantiles(self.q).arr
        
        
        
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
        
        
        
    def roc_auc(self, threshold):
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
                
        tab = CONT(self.observation_array, self.forecast_array, threshold)
        
        self.pofd = [0,tab.pofd(),1]
        self.pod  = [0, tab.pod(),1]
       
        self.ar = np.around(auc(self.pofd,self.pod),3)
        # return auc
        return self.ar
    
    
    def plot_roc(self,lead):
        '''
        Returns
        -------
        a single ROC plot for the entire extents of the data array.

        '''
        
       
        sns.set_theme(style="whitegrid")
        
        fig, ax = plt.subplots(dpi=750)
        
        ax.plot(self.pofd,self.pod,color='r',label='area under curve = {}'.format(self.ar))
        # ax.scatter(self.pofd,self.pod,color='r', marker='*')
        
        # Add a 45-degree line
        ax.plot([0, 1], [0, 1], transform=ax.transAxes, ls='--', color='k')
            
        ax.set_xlabel('PROBABILITY OF FALSE DETECTION')
        ax.set_ylabel('PROBABILITY OF DETECTION')
        
        ax.set_xlim(0,1)
        ax.set_ylim(0,1)
        ax.legend()
        if type(self.q) == int:
            plt.title("ROC at lead time: {}hrs\nagainst ICOND2EPS {}% Quantile".format(lead,self.q))
        else:
            plt.title("ROC at lead time: {}hrs\nagainst ICOND2EPS ensemble mean".format(lead))

        # Show the plot
        plt.show()
    
        
    
    
    
class ROC_totens(object):
    def __init__(self, observation_object, forecast_object, threshold):
        
        
       
                    
        
        
        # if isinstance(forecast_object,xr.DataArray):
        #      self.forecast_array = forecast_object
        # else:
            
        #     self.forecast_array = forecast_object.average
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
            
            
        # if isinstance(observation_object,xr.DataArray):
        #      self.observation_array = observation_object
        # else:
        #     if not isinstance(observation_object.average,int):
        #         self.observation_array = observation_object.average
        #     else:
        #         self.observation_array = observation_object.rtrn_arr()
                
                
        
       
        pod = []
        pofd = []
       
        
        quans = np.arange(10,100,10)
        for q in quans:
            if isinstance(forecast_object, xr.Dataset):
                per = forecast_object.to_array(dim='new')
                per = per.quantile(q=q / 100, dim='new', skipna=True)
                forecast_object = forecast_object.assign(ensemble_val=per)
                self.forecast_array = forecast_object.ensemble_val
            else:
                self.forecast_array = forecast_object.gen_quantiles(q).arr
            tab = CONT(self.observation_array, self.forecast_array, threshold)
            pod.append(tab.pod())
            pofd.append(tab.pofd())
        
        self.pod = [0] + pod + [1]
        self.pofd = [0]+ pofd+ [1]
        self.area = np.around(auc(self.pofd,self.pod),3)
    
    def plot_roc(self,lead):
        '''
        Returns
        -------
        a single ROC plot for the entire extents of the data array.

        '''
        
       
        sns.set_theme(style="whitegrid")
        
        fig, ax = plt.subplots(dpi=750)
        
        ax.plot(self.pofd,self.pod,color='r',label='auc= {} | {}hrs'.format(self.area,lead))
        # ax.scatter(self.pofd,self.pod,color='r', marker='*')
        
        # Add a 45-degree line
        ax.plot([0, 1], [0, 1], transform=ax.transAxes, ls='--', color='k')
            
        ax.set_xlabel('PROBABILITY OF FALSE DETECTION')
        ax.set_ylabel('PROBABILITY OF DETECTION')
        
        ax.set_xlim(0,1)
        ax.set_ylim(0,1)
        ax.legend()
        plt.title("ROC for ICOND2 Ensemble")

        # Show the plot
        plt.show()


# ========================================================================================================================================================


# ==============================================================================
###########################   Testing | Examples  #############################
# ==============================================================================


if __name__ == '__main__':
    from engine.fr_entities_tools import *
    
    # rad = Observation("C:/Project/radRW_juli21.nc", 24).gen_observation_field()
    # eps = Ensemble_run("C:/Project/icond2eps_3_juli21.nc", 1, 24).gen_ensemble_field(95)
    
    # roc = ROC(eps, rad)
    # area = roc.roc_auc()
    # roc.plot_roc()
    
    # mugliz = "D:/Erasmus_FRM/05.Masterarbeit/03.Bearbeitung/01.Code/Workspace/shp/Mugliz/mugliz_cats"
    # rad = Observation("C:/Project/radRW_juli21.nc", 24).limit_to_shp(mugliz+".shp")
    # rad.arr.isel(time=5).plot()
    
    # eps = Ensemble_run("C:/Project/icond2eps_3_juli21.nc", 1, 24).eps_accelerate_by_shp(mugliz+".shp").gen_ensemble_field(95)
    # eps.arr.isel(time=5).plot()
    # roc = ROC(rad,eps)
    # area = roc.roc_auc()
    # roc.plot_roc()
    
    # rad = Observation("//vs-grp07.zih.tu-dresden.de/howa/work/students/Mohamed_Elghorab/netCDFs/fertig/radRW_ICO.nc").gen_observation_field().extract_by_shp("shp/Mugliz/mugliz_cats.shp").avg_areal_prec()
    # icond2 = Ensemble_run("//vs-grp07.zih.tu-dresden.de/howa/work/students/Mohamed_Elghorab/netCDFs//3/3hour_icond2eps.nc").eps_extract_by_shp("shp/Mugliz/mugliz_cats.shp")
    # .avg_areal_prec()
    # ROC(icond2, rad).roc_auc().plot_roc()
    # rad = Observation("//vs-grp07.zih.tu-dresden.de/howa/work/students/Mohamed_Elghorab/netCDFs/fertig/radRW_ICO.nc").gen_observation_field().extract_by_shp("shp/Mugliz/mugliz_cats.shp").avg_areal_prec()
    # icond2 = Deterministic_run("//vs-grp07.zih.tu-dresden.de/howa/work/students/Mohamed_Elghorab/netCDFs//3/3hour_icond2.nc").gen_deterministic_field().extract_by_shp("shp/Mugliz/mugliz_cats.shp").avg_areal_prec()
    # icond2eps = Ensemble_run("//vs-grp07.zih.tu-dresden.de/howa/work/students/Mohamed_Elghorab/netCDFs//3/3hour_icond2eps.nc").eps_extract_by_shp("shp/Mugliz/mugliz_cats.shp").avg_areal_prec().gen_quantiles(95)
    
    # rad = Observation("C:/netCDFs/fertig/radRW_ICO.nc").gen_observation_field().extract_by_shp("shp/Mugliz/mugliz_cats.shp").aggr_temporal(3).avg_areal_prec()
    # icond2 = Deterministic_run("C:/netCDFs/3/3hour_icond2.nc").gen_deterministic_field().extract_by_shp("shp/Mugliz/mugliz_cats.shp").avg_areal_prec()
    # icond2eps = Ensemble_run("C:/netCDFs/3/3hour_icond2eps.nc").eps_extract_by_shp("shp/Mandau/egp6741491_Zittau_5_TEZG_neu.shp").avg_areal_prec()
    
    # b = ROC(rad,icond2eps,90)
    # b.roc_auc(quantile=90)
    # b.plot_roc(21)    
    # a = icond2eps.gen_quantiles(90)
    
    
    # b = ROC(rad,icond2eps,90)
    # b.roc_auc(2.5)
    # b.plot_roc(3)
    # pod = pod + [1]
    # pofd = pofd+ [1]        
    # plt.plot(pofd, pod)
    
    c = ROC_totens(rad, icond2eps, 2.5)
    c.plot_roc(3)
