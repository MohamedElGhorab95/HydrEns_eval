# -*- coding: utf-8 -*-
"""
Created on Wed Apr 26 16:23:18 2023

@author: M Elghorab
"""
import xarray as xr
import xskillscore as xs
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd

class ROC(object):
    def __init__(self, observation_object, forecast_object):
        
        if isinstance(forecast_object,xr.DataArray):
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
       
        self.pod  = () # probability of detection
        self.pofd = () # probability of false detection
        
        # selecting the same time window
        
        missing_times = set(self.observation_array.time.values) - set(self.forecast_array.time.values)
        
        obs_df = self.observation_array.to_dataframe()
        obs_df = obs_df[~obs_df.index.isin(missing_times)]

        # Convert DataFrame back to DataArray
        self.observation_array = obs_df.to_xarray()['rainfall radar observation | Radolan-RW']
        
        
        
        
    def roc_auc(self, calculation_extents="global",quantile = None):
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
        
       
        if isinstance(self.observation_array, xr.DataArray) and isinstance(self.forecast_array, xr.DataArray):
            # creating bins to categorize the data
            # get the upper and lower bounds of the rainfall fields
            upper = max(self.forecast_array.max(skipna=True),self.observation_array.max(skipna=True)).values
            lower = min(self.forecast_array.min(skipna=True),self.observation_array.min(skipna=True)).values
            bins = 199
            bin_width = (upper-lower)/bins
            
            category_edges = np.arange(lower, upper + bin_width, bin_width)
            
            if calculation_extents != "global":
                # auc per cell
               self.pofd, self.pod, auc = xs.roc(self.observation_array, self.forecast_array, bin_edges= category_edges,dim='time' ,return_results='all_as_tuple')
            else:
                # auc globally
                if "lat" in self.observation_array.coords:
                    self.pofd, self.pod, auc = xs.roc(self.observation_array, self.forecast_array, bin_edges= category_edges,dim=['time','lat','lon'] ,return_results='all_as_tuple')
                else:
                    self.pofd, self.pod, auc = xs.roc(self.observation_array, self.forecast_array, bin_edges= category_edges,dim=['time'] ,return_results='all_as_tuple')
        else:
            members = list(self.forecast_array.variables)[3:]
            pod = []
            pofd = []
            area = []
            for mem in members:
                r = ROC(self.observation_array,self.forecast_array[mem])
                area.append(float(r.roc_auc()))
                pod.append(r.pod.values)
                pofd.append(r.pofd.values)
                
            
            if quantile == "mean":
                auc = np.mean(area)
                # self.pod = [np.mean([arr[i] for arr in pod]) for i in range(len(pod[1]))]
                # self.pofd = [np.mean([arr[i] for arr in pofd]) for i in range(len(pofd[1]))]

            else:
                auc = np.quantile(area, (quantile/100))
            
                # self.pod = [np.quantile([arr[i] for arr in pod], q=quantile/100) for i in range(len(pod[1]))]
                # self.pofd = [np.quantile([arr[i] for arr in pofd], q=quantile/100) for i in range(len(pofd[1]))]

        self.ar = auc
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
        
        ax.plot(self.pofd,self.pod,color='r')
        # ax.scatter(self.pofd,self.pod,color='r', marker='*')
        
        # Add a 45-degree line
        ax.plot([0, 1], [0, 1], transform=ax.transAxes, ls='--', color='k')
            
        ax.set_xlabel('PROBABILITY OF FALSE DETECTION')
        ax.set_ylabel('PROBABILITY OF DETECTION')
        
        ax.set_xlim(0,1)
        ax.set_ylim(0,1)
        
        plt.title("ROC at lead time: {}hrs\nagainst ICOND2EPS 90% Quantile".format(lead))

        # Show the plot
        plt.show()
        



# ========================================================================================================================================================


# ==============================================================================
###########################   Testing | Examples  #############################
# ==============================================================================


if __name__ == '__main__':
    from fr_entities_tools import *
    
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
    
    rad = Observation("C:/netCDFs/fertig/radRW_ICO.nc").gen_observation_field().extract_by_shp("shp/Mugliz/mugliz_cats.shp").aggr_temporal(3).avg_areal_prec()
    # icond2 = Deterministic_run("C:/netCDFs/3/3hour_icond2.nc").gen_deterministic_field().extract_by_shp("shp/Mugliz/mugliz_cats.shp").avg_areal_prec()
    icond2eps = Ensemble_run("C:/netCDFs/21/21hour_icond2eps.nc").eps_extract_by_shp("shp/Mugliz/mugliz_cats.shp").avg_areal_prec()
    
    b = ROC(rad,icond2eps)
    b.roc_auc(quantile=90)
    b.plot_roc(21)    
    