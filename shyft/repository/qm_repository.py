from . import interfaces
from shyft import api

class QMRepositoryError(Exception):
    pass

class QMRepository(interfaces.GeoTsRepository):
    def __init__(self, bbox, start_time, qm_interp_param, repo_prior_idx, qm_cfg_params):
        self.qm_cfg_params = qm_cfg_params
        # [{'repo': repo_obj, 'period': None, 'w': (850,300), 'ens': False},
        #                       {'repo': repo_obj, 'period': api.deltahours(6*24), 'w': (500,100), 'ens': False},
        #                       {'repo': repo_obj, 'period': None, 'w': (6,), 'ens': True}]
        self.repo_prior_idx = repo_prior_idx
        self.start_time = start_time
        self.bbox = bbox

    def _read_fcst(self, geo_ts_names):
        # Read forecast from repository
        raw_fcst_lst = []
        for i in range(len(self.qm_cfg_params)):
            nb_of_fcst = len(self.qm_cfg_params[i]['w'])
            is_ens = self.qm_cfg_params[i]['ens']
            if is_ens:
                #raw_fcst = self.qm_cfg_params[i]['repo'].get_forecast_ensemble(self.start_time, self.bbox)  # we read only one ensemble for now until ec_concat for ensemble is ready
                raw_fcst = self.qm_cfg_params[i]['repo'].get_forecast_ensemble(geo_ts_names, geo_location_criteria = self.bbox)
                # when we have get_forecast_ensembles (read multiple ensemble sets at once) method available
                # in arome_concat and ec_concat OR arome raw and ec raw then we use that
            else:
                # if using either arome_concat or ec_concat
                raw_fcst = self.qm_cfg_params[i]['repo'].get_forecasts(geo_ts_names,
                               {'latest_available_forecasts': {'number of forecasts': nb_of_fcst,
                                                               'forecasts_older_than': self.start_time}},
                               geo_location_criteria={'bbox': self.bbox})
                # if using arome raw and ec raw then we either need to call get_forecast two times OR
                # we need to make a get_forecasts method which internally fetches from multiples files
            raw_fcst_lst.append(raw_fcst)

        return raw_fcst_lst

    def _read_prior(self):
        self.prior = self.raw_fcst[self.repo_prior_idx][-1]  # to get the last (latest) forecast of multiple fcsts

    def _prep_fcst(self, raw_fcst_lst):
        # Make the time_axis which all fcsts should be projected to
        raw_fcst_lst = self._read_fcst()
        prep_fcst_lst = []
        for i in range(len(self.qm_cfg_params)):
            period = self.qm_cfg_params[i]['period']
            # Check if clipping is needed according to 'period'
            if period is None:  # No clipping av ts
                #
                prep_fcst_lst.append(raw_fcst_lst[i])
            else:  # clip ts according to 'period'
                # This intermediate step is necessary since the current qm core function does not take in arguments
                # which specify which portion of each forecast to use.
                # get orig time_axis from ts in fcst -> ta_orig
                # make a time_axis which is limited to period -> ta_clip
                # call GeoTsVector.average(ta_clip)
                # append clipped fcst to raw_fcst_lst
                pass
                #prep_fcst_lst.append(prep_fcst)

    def _prep_prior(self):
        pass

    def get_forecast_ensembles(self):
        # should replace 'get_forecast_ensemble' in the long-run
        pass

    def get_forecast_ensemble(self, input_source_types, utc_period,
                              t_c, geo_location_criteria=None):
        """
        Parameters
        ----------
        input_source_types: list
            List of source types to retrieve (precipitation, temperature, ...)
        utc_period: api.UtcPeriod
            The utc time period that should (as a minimum) be covered.
        t_c: long
            Forecast specification; return newest forecast older than t_c.
        geo_location_criteria: object
            Some type (to be decided), extent (bbox + coord.ref).

        Returns
        -------
        ensemble: list of same type as get_timeseries
        Important notice: The returned forecast time-series should at least cover the
            requested period. It could return *more* data than in
            the requested period, but must return sufficient data so
            that the f(t) can be evaluated over the requested period.
        """
        pass