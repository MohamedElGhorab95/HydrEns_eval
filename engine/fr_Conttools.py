import numpy as np
import xarray as xr
import xskillscore as xs
from sklearn.metrics import f1_score
from sklearn.metrics import fbeta_score
import inspect

class CONT(object):

    def __init__(self,  observation_object, forecast_object, forecast_threshold, observation_threshold=None):
        '''
        Parameters
        ----------
        forecast_object : forecast type data
            forecast data array.
        observation_object :observation type data
            observation data array.
        forecast_threshold : float or list
            rainfall warning threshold.
        observation_threshold : float or list, optional
            warning thresholds to be used only
            in the case if it is different than
            the forecast thresholds. The default is None.

        Returns
        -------
         Contingency type object
            Has all the methods after https://xskillscore.readthedocs.io/en/stable/api/xskillscore.Contingency.html .

        '''

         # checking whether the data is xarray type
         
        # if isinstance(forecast_object,xr.DataArray):
        #      self.forecast_array = forecast_object
        # else:

        # # checking if the inputs are objects with fields or not i.e. checking if they already have an array(field) attribure

        #     # if 'arr' in dir(forecast_object):
        #     if forecast_object.dssid ==2: # deterministic run   
        #         self.forecast_array = forecast_object.rtrn_arr()
    
        #     if 'gen_deterministic_field' in dir(forecast_object):
        #         if forecast_object.dssid == 0:
        #             self.forecast_array = forecast_object.gen_deterministic_field().rtrn_arr()
            
        #         if not isinstance(forecast_object.average, int):
        #             self.forecast_array = forecast_object.average
        
        if isinstance(forecast_object,xr.DataArray):
             self.forecast_array = forecast_object
        else:
            
            self.forecast_array = forecast_object.average
        






        # if isinstance(observation_object,xr.DataArray):
        #      self.observation_array = observation_object
        # else:
        #     if 'arr' in dir(observation_object):
        #         self.observation_array = observation_object.rtrn_arr()
        #     else:
        #         self.observation_array = observation_object.gen_observation_field().rtrn_arr()
        
        # if not isinstance(observation_object.average,int):
        #     self.observation_array = observation_object.average
        
        if isinstance(observation_object,xr.DataArray):
             self.observation_array = observation_object
        else:
            if not isinstance(observation_object.average,int):
                self.observation_array = observation_object.average
            else:
                self.observation_array = observation_object.rtrn_arr()
        
        
        self.observation_threshold = observation_threshold
        self.forecast_threshold = forecast_threshold

        # selecting the same time window
        # selecting the same time window
        missing_times = set(self.observation_array.time.values) - set(self.forecast_array.time.values)
        
        obs_df = self.observation_array.to_dataframe()
        obs_df = obs_df[~obs_df.index.isin(missing_times)]

        # Convert DataFrame back to DataArray
        self.observation_array = obs_df.to_xarray()['rainfall radar observation | Radolan-RW']
        
        
        # self.forecast_array = self.forecast_array.combine_first(xr.DataArray(0, dims=('time',), coords={'time': list(missing_times)}))
 
        # start = max(self.observation_array.time[0], self.forecast_array.time[0])
        # end = min(self.observation_array.time[-1], self.forecast_array.time[-1])
        # # slicing both arrays to the same time frame
        # self.observation_array = self.observation_array.sel(time=slice(start, end))
        # self.forecast_array = self.forecast_array.sel(time=slice(start, end))

        # define the warning threshold
        # optional in case two different thresholds exist for forecast and observation
        if self.observation_threshold == None:
            self.observation_threshold = self.forecast_threshold

        # create the first dichotomous (2-category booleean array)
        # as type changes from true/false to 1/0
        self.bol_fr = (self.forecast_array >= self.forecast_threshold).astype(np.uint8)
        self.bol_obs = (self.observation_array >= self.observation_threshold).astype(np.uint8)

        # create the contingfency table

        # define the category edges
        cat_eds = np.linspace(0, 2, 3)  # this creates two bins 0-1 for the False and 1-2 for the true values

        # checking whether the data is one or two dimensional before creating the table

        if len(self.observation_array.coords) > 1 and "lat" in self.observation_array.coords:  # 2Dimensional data
            self.tab = xs.Contingency( self.bol_obs,self.bol_fr, cat_eds, cat_eds, dim=['lat', 'lon', 'time'])

        elif isinstance(self.forecast_array, xr.Dataset):
            pass
            
        else:
            self.tab = xs.Contingency( self.bol_obs,self.bol_fr, cat_eds, cat_eds, dim=['time'])

    # =============================================================================================================================

    def print_cont_tab(self):
        '''


        Returns
        -------
        str
            prints the contingency table values.

        '''

        return print('Hits    {}  | False Alarms {}\nMisses {} | Correct Negatives {}'.format(int(self.tab.hits()),
                                                                                              int(self.tab.false_alarms()),
                                                                                              int(self.tab.misses()),
                                                                                              int(self.tab.correct_negatives())))
    
    
    def quant(self, quantile):
        members = list(self.forecast_array.variables)[3:]
    
        metric = []
        calling_function = inspect.currentframe().f_back.f_code.co_name
        
        for mem in members:
            r = CONT(self.observation_array, self.forecast_array[mem], self.forecast_threshold)
            
            if calling_function == 'hits':
                metric.append(r.hits())
            elif calling_function == 'misses':
                metric.append(r.misses())
            elif calling_function == 'false_alarms':
                metric.append(r.false_alarms())
            elif calling_function == 'correct_negatives':
                metric.append(r.correct_negatives())
            elif calling_function == 'sr':
                metric.append(r.sr())
            elif calling_function == 'far':
                metric.append(r.far())
            elif calling_function == 'pod':
                metric.append(r.pod())
            elif calling_function == 'pofd':
                metric.append(r.pofd())
            elif calling_function == 'fbias':
                metric.append(r.fbias())
            elif calling_function == 'csi':
                metric.append(r.csi())
            elif calling_function == 'pss':
                metric.append(r.pss())
            elif calling_function == 'f1':
                metric.append(r.f1())
            elif calling_function == 'f2':
                metric.append(r.f2())
            elif calling_function == 'acc':
                metric.append(r.acc())
            # Add more conditions for other calling functions if needed
            
        if quantile == "mean":
            metric = np.mean(metric)
        else:
            metric = np.quantile(metric, (quantile/100))
        
        return metric
        
   
   

    # =============================================================================================================================

    # basic contingency table metrics getters

    # def hits(self):
    #     return self.tab.hits()
    def hits(self,q=None):
        if isinstance(self.observation_array, xr.DataArray) and isinstance(self.forecast_array, xr.DataArray):
            return self.tab.hits()
        else:
            return int(self.quant(q))
        
    def misses(self,q=None):
        if isinstance(self.observation_array, xr.DataArray) and isinstance(self.forecast_array, xr.DataArray):
            return self.tab.misses()
        else:
            return int(self.quant(q))

    def false_alarms(self,q=None):
        if isinstance(self.observation_array, xr.DataArray) and isinstance(self.forecast_array, xr.DataArray):
            return self.tab.false_alarms()
        else:
            return int(self.quant(q))

    def correct_negatives(self,q=None):
        if isinstance(self.observation_array, xr.DataArray) and isinstance(self.forecast_array, xr.DataArray):
            return self.tab.correct_negatives()
        else:
            return int(self.quant(q))

    # =============================================================================================================================

    # contingency table performance metrics

    def acc(self,q=None):
        '''
        .. math:: Accuracy = \\frac{hits+correct negatives}{total number of instances}

        Returns
        -------
        float
            SUCCESS RATIO OF THE FORECAST
        The success ratio shows the fraction of the correctly identified “Yes” instances over the total number of forecasted “Yes” instances.

        '''
        if isinstance(self.observation_array, xr.DataArray) and isinstance(self.forecast_array, xr.DataArray):
            return float((self.hits() + self.correct_negatives()) / (self.hits() + self.misses() + self.false_alarms() + self.correct_negatives()))
        else:
            return self.quant(q)


    
    def sr(self,q=None):
        '''
        .. math:: SR = \\frac{hits}{hits+false~alarms}

        Returns
        -------
        float
            SUCCESS RATIO OF THE FORECAST
        The success ratio shows the fraction of the correctly identified “Yes” instances over the total number of forecasted “Yes” instances.

        '''
        if isinstance(self.observation_array, xr.DataArray) and isinstance(self.forecast_array, xr.DataArray):
            return float(self.tab.success_ratio())
        else:
            return self.quant(q)

    def far(self,q=None):
        '''
        .. math:: FAR = \\frac{false~alarms}{hits+false~alarms}

        Returns
        -------
        float
            FAR RATIO OF THE FORECAST
        The false alarm ratio indicates the ratio of the false alarms to the total number of forecasted events.

        '''
        if isinstance(self.observation_array, xr.DataArray) and isinstance(self.forecast_array, xr.DataArray):
            return float(self.tab.false_alarm_ratio())
        else:
            return self.quant(q)

    def pod(self,q=None):
        '''
        .. math:: POD = \\frac{hits}{hits+misses}

        Returns
        -------
        float
            PROBABILITY OF DETECTION OF THE FORECAST
        The probability of detection indicates the ratio of the correctly forecasted hits to total number of observed events.

        '''
        if isinstance(self.observation_array, xr.DataArray) and isinstance(self.forecast_array, xr.DataArray):
            return float(self.hits() / (self.hits() + self.misses()))
        else:
            return self.quant(q)

    def pofd(self,q=None):
        '''
        .. math:: POFD = \\frac{false~alarms}{false~alarms+correct~negatives}

        Returns
        -------
        float
            PROBABILITY OF FALSE DETECTION OF THE FORECAST
        The probability of false detection indicates the ratio of the false alarms to the total number of observed non-event instances.

        '''
        if isinstance(self.observation_array, xr.DataArray) and isinstance(self.forecast_array, xr.DataArray):
            return float(self.false_alarms() / (self.false_alarms() + self.correct_negatives()))
        else:
            return self.quant(q)
            
    def fbias(self,q=None):
        '''
        .. math:: FBIAS = \\frac{hits+false~alarms}{hits+misses}

        Returns
        -------
        float
            FREQUENCY BIAS OF THE FORECAST
        The frequency bias compares the number of forecasted events to the number of observed events. It indicated whether the forecast has a trend of over or under forecasting.

        '''
        if isinstance(self.observation_array, xr.DataArray) and isinstance(self.forecast_array, xr.DataArray):
            return float((self.hits() + self.false_alarms()) / (self.hits() + self.misses()))
        else:
            return self.quant(q)
    
    def csi(self,q=None):
        '''
        .. math:: CSI = \\frac{hits}{Total~sum~of~instances}

        Returns
        -------
        float
            CRITICAL SUCCESS INDEX OF THE FORECAST
        The critical success index explicitly indicates to what degree the forecast product can capture flooding events.

        '''
        if isinstance(self.observation_array, xr.DataArray) and isinstance(self.forecast_array, xr.DataArray):
            return float(self.hits() / (self.hits() + self.false_alarms() + self.misses()))
        else:
            return self.quant(q)
        
    def pss(self, q=None):
        '''
        .. math:: PSS = \\frac{hits}{hits+misses} - \\frac{false~alarms}{false~alarms+correct~negatives}

        Returns
        -------
        float
            PIERCE'S SKILL SCORE OF THE FORECAST
        Peirce's skill score shows how well the forecast discriminates between the positive and the negative event instances

        '''
        if isinstance(self.observation_array, xr.DataArray) and isinstance(self.forecast_array, xr.DataArray):
            return float(self.tab.peirce_score())
        else:
            return self.quant(q)
    
    def f1(self, q=None):
        '''
        .. math:: F1 = 2*\\frac{Precision*POD}{Precision+POD} 

        Returns
        -------
        float
            F1 SKILL SCORE OF THE FORECAST
        The F1 score is a measure of a model's accuracy in binary classification tasks, combining both precision and recall into a single metric. It is the harmonic mean of precision and recall, providing a balanced evaluation of a model's performance
        where
            .. math:: Precision = SR = \\frac{hits}{hits+false~alarms}
        '''
        if isinstance(self.observation_array, xr.DataArray) and isinstance(self.forecast_array, xr.DataArray):
            # # Convert the DataArrays to numpy arrays
            predicted_np = self.bol_fr.values
            actual_np = self.bol_obs.values
            
            # Calculate the F1 score
            f1 = f1_score(actual_np, predicted_np)
            
                
            return float(f1)
        else:
            return self.quant(q)
    
    def f2(self,q=None):
        '''
        .. math:: F2 = 5*\\frac{Precision*POD}{4*Precision+POD} 

        Returns
        -------
        float
            F2 SKILL SCORE OF THE FORECAST
         F2 score emphasizes the importance of correctly identifying positive instances, which is particularly useful in scenarios where false negatives are more critical than false positives.
        where
            .. math:: Precision = SR = \\frac{hits}{hits+false~alarms}
        '''
        if isinstance(self.observation_array, xr.DataArray) and isinstance(self.forecast_array, xr.DataArray):
            # Convert the DataArrays to numpy arrays
            predicted_np = self.bol_fr.values
            actual_np = self.bol_obs.values
            
            # Calculate the F1 score
            f2 = fbeta_score(actual_np, predicted_np, beta=2)
           
                
            return float(f2)
        else:
            return self.quant(q)

# ========================================================================================================================================================


# ==============================================================================
###########################   Testing | Examples  #############################
# ==============================================================================


if __name__ == '__main__':
    from engine.fr_entities_tools import *


    def test_CONT():
        '''
        Testing the Cont Class

        '''

        # # instantiating an Xarray type object for radolan
        rad = Observation("C:/Project/radRW_juli21.nc", 3)
        # # # instantiating objects for Icond2 for a horizon of 9 hours
        icond2 = Deterministic_run("C:/Project/icond2_3_juli21.nc", 1, 9)
        icond2eps = Ensemble_run("C:/Project/icond2eps_3_juli21.nc", 1, 9)

        # # instantiating a contingency table object from icond2 and rad, with a
        # # threshold of 10 mm/hr
        cont = CONT(icond2, rad, 10)
        # instantiating a contingency table object from icond2eps' 95th quantile
        # and rad, with a threshold of 10 mm/hr
        conteps = CONT(icond2eps, rad, 10, 95)

        # instantiating for MAP

        # # instantiating an Xarray type object for radolan
        radav = Observation("C:/Project/radRW_juli21.nc"
                            , 9).gen_observation_field().avg_areal_prec()

        # # # instantiating an  objects for Icond2 for a horizon of 9 hours
        # # and 95th quantile for the ensemble
        icond2av = Deterministic_run("C:/Project/icond2_3_juli21.nc", 2
                                     , 9).gen_deterministic_field().avg_areal_prec()
        icond2epsav = Ensemble_run("C:/Project/icond2eps_3_juli21.nc", 2
                                   , 9).gen_ensemble_field(95).avg_areal_prec()

        # # instantiating a contingency table object from icond2 and rad, with a
        # # threshold of 10 mm/hr
        contav = CONT(icond2av, radav, 10)
        # # instantiating a contingency table object from icond2eps' 95th quantile
        # # and rad, with a threshold of 10 mm/hr
        contepsav = CONT(icond2epsav, radav, 10)

        # # some evaluation metrics

        cont.sr()
        # cont.far()
        # cont.pod()
        # cont.pofd()
        cont.fbias()
        cont.csi()
        cont.pss()

        return cont.print_cont_tab(), conteps.print_cont_tab(), contav.print_cont_tab(), contepsav.print_cont_tab()


    # test_CONT()

    def test_CONT_sub():
        '''
        Testing the Cont Class for subareas

        '''

        # # instantiating an Xarray type object for radolan
        rad = Observation("C:/Project/radRW_juli21.nc", 3).extract_by_coords(50.02, 51.01, 11.66, 12.68)
        # # # instantiating objects for Icond2 for a horizon of 9 hours
        icond2 = Deterministic_run("C:/Project/icond2_3_juli21.nc", 1, 9).extract_by_coords(50.02, 51.01, 11.66, 12.68)
        icond2eps = Ensemble_run("C:/Project/icond2eps_3_juli21.nc", 1, 9).eps_extract_by_coords(50.02, 51.01, 11.66,
                                                                                                 12.68)

        # # instantiating a contingency table object from icond2 and rad, with a
        # # threshold of 10 mm/hr
        cont = CONT(icond2, rad, 10)
        # instantiating a contingency table object from icond2eps' 95th quantile
        # and rad, with a threshold of 10 mm/hr
        conteps = CONT(icond2eps, rad, 10, ensemble_averaging_method=95)

        # instantiating for MAP

        # # instantiating an Xarray type object for radolan
        radav = Observation("C:/Project/radRW_juli21.nc"
                            , 9).extract_by_coords(50.02, 51.01, 11.66, 12.68).gen_observation_field().avg_areal_prec()

        # # # instantiating an  objects for Icond2 for a horizon of 9 hours
        # # and 95th quantile for the ensemble
        icond2av = Deterministic_run("C:/Project/icond2_3_juli21.nc", 2
                                     , 9).extract_by_coords(50.02, 51.01, 11.66,
                                                            12.68).gen_deterministic_field().avg_areal_prec()
        icond2epsav = Ensemble_run("C:/Project/icond2eps_3_juli21.nc", 2
                                   , 9).eps_extract_by_coords(50.02, 51.01, 11.66, 12.68).gen_ensemble_field(
            95).avg_areal_prec()

        # # instantiating a contingency table object from icond2 and rad, with a
        # # threshold of 10 mm/hr
        contav = CONT(icond2av, radav, 10)
        # # instantiating a contingency table object from icond2eps' 95th quantile
        # # and rad, with a threshold of 10 mm/hr
        contepsav = CONT(icond2epsav, radav, 10)

        # # some evaluation metrics

        # cont.sr()
        # cont.far()
        # cont.pod()
        # cont.pofd()
        # cont.fbias()
        # cont.csi()
        # cont.pss()

        return cont.print_cont_tab(), conteps.print_cont_tab(), contav.print_cont_tab(), contepsav.print_cont_tab()


    # test_CONT_sub()

    def test_CONT_cat():
        '''
        Testing the Cont Class for catchments
        für Weiße Elster Einzugsgebiet

        '''

        # # instantiating an Xarray type object for radolan
        rad = Observation("C:/Project/radRW_juli21.nc", 3).extract_by_shp(
            "C:/Project/shp/WeiseElster/OelsnitzTEZG_DHDN.shp")
        # # # instantiating objects for Icond2 for a horizon of 9 hours
        icond2 = Deterministic_run("C:/Project/icond2_3_juli21.nc", 1, 9).extract_by_shp(
            "C:/Project/shp/WeiseElster/OelsnitzTEZG_DHDN.shp")
        icond2eps = Ensemble_run("C:/Project/icond2eps_3_juli21.nc", 1, 9).eps_accelerate_by_shp(
            "C:/Project/shp/WeiseElster/OelsnitzTEZG_DHDN").gen_ensemble_field(95).eps_extract_by_shp(
            "C:/Project/shp/WeiseElster/OelsnitzTEZG_DHDN.shp")
        # # instantiating a contingency table object from icond2 and rad, with a
        # # threshold of 10 mm/hr
        cont = CONT(icond2, rad, 10)
        # instantiating a contingency table object from icond2eps' 95th quantile
        # and rad, with a threshold of 10 mm/hr
        conteps = CONT(icond2eps, rad, 10)

        # instantiating for MAP

        # # instantiating an Xarray type object for radolan
        radav = Observation("C:/Project/radRW_juli21.nc"
                            , 9).extract_by_shp(
            "C:/Project/shp/WeiseElster/OelsnitzTEZG_DHDN.shp").gen_observation_field().avg_areal_prec()

        # # # instantiating an  objects for Icond2 for a horizon of 9 hours
        # # and 95th quantile for the ensemble
        icond2av = Deterministic_run("C:/Project/icond2_3_juli21.nc", 2
                                     , 9).extract_by_shp(
            "C:/Project/shp/WeiseElster/OelsnitzTEZG_DHDN.shp").gen_deterministic_field().avg_areal_prec()
        icond2epsav = Ensemble_run("C:/Project/icond2eps_3_juli21.nc", 2
                                   , 9).eps_accelerate_by_shp(
            "C:/Project/shp/WeiseElster/OelsnitzTEZG_DHDN").gen_ensemble_field(95).eps_extract_by_shp(
            "C:/Project/shp/WeiseElster/OelsnitzTEZG_DHDN.shp").gen_ensemble_field(95).avg_areal_prec()

        # # instantiating a contingency table object from icond2 and rad, with a
        # # threshold of 10 mm/hr
        contav = CONT(icond2av, radav, 10)
        # # instantiating a contingency table object from icond2eps' 95th quantile
        # # and rad, with a threshold of 10 mm/hr
        contepsav = CONT(icond2epsav, radav, 10)

        # # some evaluation metrics

        # cont.sr()
        # cont.far()
        # cont.pod()
        # cont.pofd()
        # cont.fbias()
        # cont.csi()
        cont.pss()

        # return cont.print_cont_tab(), conteps.print_cont_tab()
        return cont.print_cont_tab(), conteps.print_cont_tab(), contav.print_cont_tab(), contepsav.print_cont_tab()


    # test_CONT_cat()
    
    
    
    
    
    
    rad = Observation("C:/netCDFs/fertig/radRW_ICO.nc").gen_observation_field().extract_by_shp("shp/Mugliz/mugliz_cats.shp").aggr_temporal(3).avg_areal_prec()
    # icond2 = Deterministic_run("C:/netCDFs/3/3hour_icond2.nc").gen_deterministic_field().extract_by_shp("shp/Mugliz/mugliz_cats.shp").avg_areal_prec()
    icond2eps = Ensemble_run("C:/netCDFs/3/3hour_icond2eps.nc").eps_extract_by_shp("shp/Mugliz/mugliz_merged.shp").avg_areal_prec()
    # .gen_quantiles(90)
    
    x = CONT(rad, icond2eps, 3)
    # x.hits(q=70)
    # x.csi(95)
    # x.f1(95)
    x.fbias('mean')
    x.pss(25)
    # x.sr()
    # a.f1()
    # # print(a.f1())
    # print(a.f2())
    x.acc(50)    
    x.pod()
    x.pofd()
    
    
    variable_names = list(icond2eps.average.data_vars)
    member_arrays = [icond2eps.average[var_name] for var_name in variable_names]
    forecast_dataarray = xr.concat(member_arrays, dim='member')
    
    crps = xs.crps_ensemble(rad.average, forecast_dataarray, dim='time')
