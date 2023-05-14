import numpy as np
import xskillscore as xs


class CONT(object):

    def __init__(self, forecast_object, observation_object, forecast_threshold, observation_threshold=None):
        '''
        Parameters
        ----------
        forecast_array : xarray data array
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

     

        # checking if the inputs are objects with fields or not i.e. checking if they already have an array(field) attribure

        if 'arr' in dir(forecast_object):
            self.forecast_array = forecast_object.rtrn_arr()

        if 'gen_deterministic_field' in dir(forecast_object):
            if forecast_object.dssid == 0:
                self.forecast_array = forecast_object.gen_deterministic_field().rtrn_arr()

        if 'arr' in dir(observation_object):
            self.observation_array = observation_object.rtrn_arr()
        else:
            self.observation_array = observation_object.gen_observation_field().rtrn_arr()

        self.observation_threshold = observation_threshold
        self.forecast_threshold = forecast_threshold

        # selecting the same time window
        start = max(self.observation_array.time[0], self.forecast_array.time[0])
        end = min(self.observation_array.time[-1], self.forecast_array.time[-1])
        # slicing both arrays to the same time frame
        self.observation_array = self.observation_array.sel(time=slice(start, end))
        self.forecast_array = self.forecast_array.sel(time=slice(start, end))

        # define the warning threshold
        # optional in case two different thresholds exist for forecast and observation
        if self.observation_threshold == None:
            self.observation_threshold = self.forecast_threshold

        # create the first dichotomous (2-category booleean array)
        # as type changes from true/false to 1/0
        bol_fr = (self.forecast_array >= self.forecast_threshold).astype(np.uint8)
        bol_obs = (self.observation_array >= self.observation_threshold).astype(np.uint8)

        # create the contingfency table

        # define the category edges
        cat_eds = np.linspace(0, 2, 3)  # this creates two bins 0-1 for the False and 1-2 for the true values

        # checking whether the data is one or two dimensional before creating the table

        if len(self.observation_array.coords) > 1 and "lat" in self.observation_array.coords:  # 2Dimensional data
            self.tab = xs.Contingency( bol_obs,bol_fr, cat_eds, cat_eds, dim=['lat', 'lon', 'time'])

        else:
            self.tab = xs.Contingency( bol_obs,bol_fr, cat_eds, cat_eds, dim=['time'])

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

    # =============================================================================================================================

    # basic contingency table metrics getters

    def hits(self):
        return self.tab.hits()

    def misses(self):
        return self.tab.misses()

    def false_alarms(self):
        return self.tab.false_alarms()

    def correct_negatives(self):
        return self.tab.correct_negatives()

    # =============================================================================================================================

    # contingency table performance metrics

    def sr(self):
        '''
        .. math:: SR = \\frac{hits}{hits+false~alarms}

        Returns
        -------
        float
            SUCCESS RATIO OF THE FORECAST
        The success ratio shows the fraction of the correctly identified “Yes” instances over the total number of forecasted “Yes” instances.

        '''

        return float(self.tab.success_ratio())

    def far(self):
        '''
        .. math:: FAR = \\frac{false~alarms}{hits+false~alarms}

        Returns
        -------
        float
            FAR RATIO OF THE FORECAST
        The false alarm ratio indicates the ratio of the false alarms to the total number of forecasted events.

        '''
        return float(self.tab.false_alarm_ratio())

    def pod(self):
        '''
        .. math:: POD = \\frac{hits}{hits+misses}

        Returns
        -------
        float
            PROBABILITY OF DETECTION OF THE FORECAST
        The probability of detection indicates the ratio of the correctly forecasted hits to total number of observed events.

        '''
        return float(self.hits() / (self.hits() + self.misses()))

    def pofd(self):
        '''
        .. math:: POFD = \\frac{false~alarms}{false~alarms+correct~negatives}

        Returns
        -------
        float
            PROBABILITY OF FALSE DETECTION OF THE FORECAST
        The probability of false detection indicates the ratio of the false alarms to the total number of observed non-event instances.

        '''
        return float(self.false_alarms() / (self.false_alarms() + self.correct_negatives()))

    def fbias(self):
        '''
        .. math:: FBIAS = \\frac{hits+false~alarms}{hits+misses}

        Returns
        -------
        float
            FREQUENCY BIAS OF THE FORECAST
        The frequency bias compares the number of forecasted events to the number of observed events. It indicated whether the forecast has a trend of over or under forecasting.

        '''
        return float((self.hits() + self.false_alarms()) / (self.hits() + self.misses()))

    def csi(self):
        '''
        .. math:: CSI = \\frac{hits}{Total~sum~of~instances}

        Returns
        -------
        float
            CRITICAL SUCCESS INDEX OF THE FORECAST
        The critical success index explicitly indicates to what degree the forecast product can capture flooding events.

        '''
        return float(self.hits() / (self.hits() + self.correct_negatives() + self.false_alarms() + self.misses()))

    def pss(self):
        '''
        .. math:: PSS = \\frac{hits}{hits+misses} + \\frac{false~alarms}{false~alarms+correct~negatives}

        Returns
        -------
        float
            PIERCE'S SKILL SCORE OF THE FORECAST
        Peirce's skill score shows how well the forecast discriminates between the positive and the negative event instances

        '''

        return float(self.tab.peirce_score())


# ========================================================================================================================================================


# ==============================================================================
###########################   Testing | Examples  #############################
# ==============================================================================


if __name__ == '__main__':
    from fr_entities_tools import *


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

        # cont.sr()
        # cont.far()
        # cont.pod()
        # cont.pofd()
        # cont.fbias()
        # cont.csi()
        # cont.pss()

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
        # cont.pss()

        # return cont.print_cont_tab(), conteps.print_cont_tab()
        return cont.print_cont_tab(), conteps.print_cont_tab(), contav.print_cont_tab(), contepsav.print_cont_tab()


    test_CONT_cat()