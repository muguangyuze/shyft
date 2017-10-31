from . import interfaces
from shyft import api
import numpy as np

class QMRepositoryError(Exception):
    pass

class QMRepository(interfaces.GeoTsRepository):

    def __init__(self, bbox, start_time, qm_interp_param, repo_prior_idx, qm_cfg_params, qm_resolution):
        self.qm_cfg_params = qm_cfg_params
        # [{'repo': repo_obj, 'period': None, 'w': (850,300), 'ens': False},
        #                       {'repo': repo_obj, 'period': api.deltahours(6*24), 'w': (500,100), 'ens': False},
        #                       {'repo': repo_obj, 'period': None, 'w': (6,), 'ens': True}]
        self.qm_resolution = qm_resolution
        # qm_resolution = {'short': (0, api.TimeAxis(start_time, api.deltahours(1), 60)),
        #                  'mid': (1, api.TimeAxis(start_time + api.deltahours(60), api.deltahours(3), 24)),
        #                  'long': (2, api.TimeAxis(start_time + api.deltahours(60 + 3 * 24), api.deltahours(6), 24))}
        self.repo_prior_idx = repo_prior_idx
        self.qm_interp_param = qm_interp_param
        self.start_time = start_time
        self.bbox = bbox

    def _read_fcst(self, input_source_types):
        # Read forecasts from repository for each forecast group
        # Hierachy of return is list of forecast grup, list of forecasts, list of forecast members where each member
        # is a source-type keyed dictionary of geo time seies
        raw_fcst_lst = []
        for i in range(len(self.qm_cfg_params)):
            nb_of_fcst = len(self.qm_cfg_params[i]['w'])
            is_ens = self.qm_cfg_params[i]['ens']
            if is_ens:
                # raw_fcst = self.qm_cfg_params[i]['repo'].get_forecast_ensemble(self.start_time, self.bbox)
                # we read only one ensemble for now until ec_concat for ensemble is ready
                # when we have get_forecast_ensembles (read multiple ensemble sets at once) method available
                # in arome_concat and ec_concat OR arome raw and ec raw then we use that
                raw_fcst = self.qm_cfg_params[i]['repo'].get_forecast_ensemble(input_source_types, None, self.start_time,
                                                                               geo_location_criteria=self.bbox)
            else:
                # if using either arome_concat or ec_concat
                # return as list to match output from get_forecast_ensemble
                raw_fcst = [self.qm_cfg_params[i]['repo'].get_forecasts(input_source_types,
                               {'latest_available_forecasts': {'number of forecasts': nb_of_fcst,
                                                               'forecasts_older_than': self.start_time}},
                               geo_location_criteria={'bbox': self.bbox})]

                # if using arome raw and ec raw then we either need to call get_forecast two times OR
                # we need to make a get_forecasts method which internally fetches from multiples files

            # Examine first input_source_type to determine number of geo_points and use this to slice into raw_fcst to
            # return one list entry per forecast member
            fc_start_t, nb_geo_pts = np.unique([src.ts.time(0) for src in raw_fcst[0][input_source_types[0]]],
                                           return_counts=True)
            nb_geo_pts = int(nb_geo_pts[0])
            slices = [slice(i * nb_geo_pts, (i + 1) * nb_geo_pts) for i in range(nb_of_fcst)]
            forecasts = [[{src_type: memb[src_type][slc] for src_type in input_source_types} for memb in raw_fcst] for slc
                         in slices]
            # forecasts = [{src_type: memb[src_type][slc] for src_type in input_source_types} for memb in raw_fcst for
            #              slc in slices]

            raw_fcst_lst.append(forecasts)
        return raw_fcst_lst

    def _get_geo_pts(self, raw_fcst):
        # extract geo_pts from first souce vector and return as GeoPointVector
        first_fcst_member = raw_fcst[0][0]
        first_src_vec = first_fcst_member[list(first_fcst_member.keys())[0]]
        start_times, nb_geo_pts = np.unique([src.ts.time(0) for src in first_src_vec],
                                            return_counts=True)
        nb_geo_pts = nb_geo_pts[0]
        gpv = api.GeoPointVector()
        for k in range(0, nb_geo_pts):
            gpv.append(first_src_vec[k].mid_point())
        return gpv

    def _read_prior(self):
        # to get the last (latest) forecast of multiple fcsts
        pass

    def _prep_prior(self):
        pass

    def _downscaling(self, input_source_types, forecast, target_grid, ta_fixed_dt):
        # Using idw for time being
        prep_fcst = {}
        for src_type in input_source_types:
            if src_type == 'precipitation':
                # just setting som idw_params for time being
                idw_params = api.IDWPrecipitationParameter()
                idw_params.max_distance = 6000
                idw_params.max_members = 4
                idw_params.gradient_by_equation = False
                prep_fcst[src_type] = api.idw_precipitation(forecast[src_type], target_grid, ta_fixed_dt, idw_params)
            else:
                # just setting som idw_params for time being
                idw_params = api.IDWTemperatureParameter()
                idw_params.max_distance = 6000
                idw_params.max_members = 4
                idw_params.gradient_by_equation = False
                prep_fcst[src_type] = api.idw_temperature(forecast[src_type], target_grid, ta_fixed_dt, idw_params)
        return prep_fcst

    def _prep_fcst(self, input_source_types, resolution_key):
        # Identifies space-time resolution and call downscaling routene

        qm_cfg_params = self.qm_cfg_params
        raw_fcst_lst = self._read_fcst(input_source_types)
        qm_resolution_idx, ta = self.qm_resolution[resolution_key]

        # Use prior resolution as target resolution if nothing is specified
        if qm_resolution_idx is None and repo_prior_idx is not None:
            qm_resolution_idx = repo_prior_idx
        if qm_resolution_idx is not None:
            target_grid = self._get_geo_pts(raw_fcst_lst[qm_resolution_idx])
        # TODO: if qm_resolution_idx is None and repo_prior_idx is None read prior and use as target grid

        prep_fcst_lst = []
        for i in range(len(qm_cfg_params)):
            raw_fcst = raw_fcst_lst[i]
            period = qm_cfg_params[i]['period']

            # Adjust time axis end if neccesary
            if period is not None:
                period_end = self.start_time + api.deltahours(24 * period)
                ta_end = ta.total_period().end
                new_n = (min(period_end, ta_end) - ta.time(0)) // ta.fixed_dt.delta_t
                ta = api.TimeAxis(ta.time(0), ta.fixed_dt.delta_t, new_n)

            prep_fcst = [[self._downscaling(input_source_types, f_m, target_grid, ta.fixed_dt) for f_m in fct] for fct
                         in raw_fcst]
            # prep_fcst = [self._downscaling(input_source_types, f_m, target_grid, ta.fixed_dt) for f_m in raw_fcst]

            # for src_type in input_source_types:
            #     if src_type == 'precipitation':
            #         # just setting som idw_params for time being
            #         idw_params = api.IDWPrecipitationParameter()
            #         idw_params.max_distance = 6000
            #         idw_params.max_members = 4
            #         idw_params.gradient_by_equation = False
            #         prep_fcst = [{src_type: api.idw_precipitation(f_m[src_type], target_grid, ta.fixed_dt, idw_params)}
            #                      for f_m in raw_fcst]
            #     else:
            #         # just setting som idw_params for time being
            #         idw_params = api.IDWTemperatureParameter()
            #         idw_params.max_distance = 6000
            #         idw_params.max_members = 4
            #         idw_params.gradient_by_equation = False
            #
            #     prep_fcst = [{src_type: api.idw_temperature(f_m[src_type], target_grid, ta.fixed_dt, idw_params)}
            #                      for f_m in raw_fcst]

            prep_fcst_lst.append(prep_fcst)

        return prep_fcst_lst, target_grid, ta

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

        if self.repo_prior_idx is None:
            prior = self._read_prior()

        prep_fcst_lst, target_grid, ta = self._prep_fcst(input_source_types, 'long')

        for geo_pt_idx, geo_pt in enumerate(target_grid):
            for src in input_source_types:
                forecast_sets = api.TsVectorSet()
                # prior_data = api.TsVector()
                weight_sets = api.DoubleVector()
                for i, fcst_group in enumerate(prep_fcst_lst) :
                   for j, forecast in enumerate(fcst_group):
                        weight_sets.append(self.qm_cfg_params[i]['w'][j])
                        scenarios = api.TsVector()
                        for member in forecast:
                            scenarios.append(member[src][geo_pt_idx].ts)
                        forecast_sets.append(scenarios)
                        if i == self.repo_prior_idx and j==0:
                            prior_data = scenarios
                            # TODO: read prior if repo_prior_idx is None

                api.quantile_map_forecast(forecast_sets, weight_sets, prior_data, ta, self.qm_interp_param[0],
                                           self.qm_interp_param[1], True)
                # TODO: check return type
