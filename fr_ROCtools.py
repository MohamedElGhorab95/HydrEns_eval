# -*- coding: utf-8 -*-
"""
Created on Wed Apr 26 16:23:18 2023

@author: M Elghorab
"""
import xskillscore as xs
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns


class ROC(object):
    def __init__(self, observation_object, forecast_object):
        
        
        self.forecast_array = forecast_object.rtrn_arr()
        self.observation_array = observation_object.rtrn_arr()
        self.pod  = () # probability of detection
        self.pofd = () # probability of false detection
        
        # selecting the same time window
        start = max(self.observation_array.time[0], self.forecast_array.time[0])
        end = min(self.observation_array.time[-1], self.forecast_array.time[-1])
        # slicing both arrays to the same time frame
        self.observation_array = self.observation_array.sel(time=slice(start, end))
        self.forecast_array = self.forecast_array.sel(time=slice(start, end))
        
        
    def roc_auc(self, calculation_extents="global"):
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
        
       
        
        # creating bins to categorize the data
        # get the upper and lower bounds of the rainfall fields
        upper = max(self.forecast_array.max(skipna=True),self.observation_array.max(skipna=True)).values
        lower = min(self.forecast_array.min(skipna=True),self.observation_array.min(skipna=True)).values
        bins = 4
        bin_width = (upper-lower)/bins
        
        category_edges = np.arange(lower, upper + bin_width, bin_width)
        
        if calculation_extents != "global":
            # auc per cell
           self.pofd, self.pod, auc = xs.roc(self.observation_array, self.forecast_array, bin_edges= category_edges,dim='time' ,return_results='all_as_tuple')
        else:
            # auc globally
           self.pofd, self.pod, auc = xs.roc(self.observation_array, self.forecast_array, bin_edges= category_edges,dim=['time','lat','lon'] ,return_results='all_as_tuple')
        
        
        return auc
    
    
    def plot_roc(self):
        '''
        Returns
        -------
        a single ROC plot for the entire extents of the data array.

        '''
        
       
        sns.set_theme(style="whitegrid")
        
        fig, ax = plt.subplots()
        
        ax.plot(self.pofd,self.pod,color='r')
        
        # Add a 45-degree line
        ax.plot([0, 1], [0, 1], transform=ax.transAxes, ls='--', color='k')
            
        ax.set_xlabel('PROBABILITY OF FALSE DETECTION')
        ax.set_ylabel('PROBABILITY OF DETECTION')
        
        ax.set_xlim(0,1)
        ax.set_ylim(0,1)
        
        plt.title("ROC")

        # Show the plot
        plt.show()
        



# ========================================================================================================================================================


# ==============================================================================
###########################   Testing | Examples  #############################
# ==============================================================================


if __name__ == '__main__':
    from HydrEns_eval.fr_entities_tools import *
    
    rad = Observation("C:/Project/radRW_juli21.nc", 24).gen_observation_field()
    eps = Ensemble_run("C:/Project/icond2eps_3_juli21.nc", 1, 24).gen_ensemble_field(95)
    
    roc = ROC(eps, rad)
    area = roc.roc_auc()
    roc.plot_roc()
    
