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
        self.source_type_map = {"relative_humidity": api.RelHumSource,
                                "temperature": api.TemperatureSource,
                                "precipitation": api.PrecipitationSource,
                                "radiation": api.RadiationSource,
                                "wind_speed": api.WindSpeedSource}

        self.source_vector_map = {"relative_humidity": api.RelHumSourceVector,
                                "temperature": api.TemperatureSourceVector,
                                "precipitation": api.PrecipitationSourceVector,
                                "radiation": api.RadiationSourceVector,
                                "wind_speed": api.WindSpeedSourceVector}

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
                idw_params.max_distance = 15000
                idw_params.max_members = 4
                idw_params.gradient_by_equation = False
                prep_fcst[src_type] = api.idw_precipitation(forecast[src_type], target_grid, ta_fixed_dt, idw_params)
            else:
                # just setting som idw_params for time being
                idw_params = api.IDWTemperatureParameter()
                idw_params.max_distance = 15000
                idw_params.max_members = 4
                idw_params.gradient_by_equation = False
                prep_fcst[src_type] = api.idw_temperature(forecast[src_type], target_grid, ta_fixed_dt, idw_params)
        return prep_fcst

    def _prep_fcst(self, raw_fcst_lst, input_source_types, resolution_key):
        # Identifies space-time resolution and calls downscaling routine

        qm_cfg_params = self.qm_cfg_params
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
                ta_to_idw = api.TimeAxis(ta.time(0), ta.fixed_dt.delta_t, new_n)
            else:
                ta_to_idw = ta

            prep_fcst = [[self._downscaling(input_source_types, f_m, target_grid, ta_to_idw.fixed_dt) for f_m in fct]
                         for fct in raw_fcst]

            prep_fcst_lst.append(prep_fcst)

        return prep_fcst_lst, target_grid

    def _call_qm(self, prep_fcst_lst, geo_points, ta, input_source_types):

        # TODO: Extend handling to cover all cases and send out warnings
        # Check ta against interpolation start and end times
        # Simple logic for time being, should be refined for the overlap cases
        interp_start = self.qm_interp_param[0]
        interp_end = self.qm_interp_param[1]
        ta_start = ta.time(0)
        ta_end = ta.time(ta.size()-1) # start of last time step
        if interp_start > ta_end:
            interp_start = api.no_utctime
            interp_end = api.no_utctime
        if interp_end > ta_end:
            interp_end = ta_end

        dict = {}
        for src in input_source_types:
            qm_scenarios = []
            for geo_pt_idx, geo_pt in enumerate(geo_points):
                forecast_sets = api.TsVectorSet()

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

                qm_scenarios.append(api.quantile_map_forecast(forecast_sets, weight_sets, prior_data, ta,
                                                         interp_start, interp_end, True))
            dict[src] = np.array(qm_scenarios)
        # TODO: write function to extract prior info like number of scenarios
        nb_scen = dict[input_source_types[0]].shape[1]
        results = []
        for i in range(0,nb_scen):
            source_dict = {}
            for src in input_source_types:
                ts_vct = dict[src][:, i]
                vct = self.source_vector_map[src]()
                [vct.append(self.source_type_map[src](geo_pt, ts)) for geo_pt, ts in zip(geo_points, ts_vct)]
                # Alternatives:
                # vct[:] = [self.source_type_map[src](geo_pt, ts) for geo_pt, ts in zip(geo_points, ts_vct)]
                # vct = self.source_vector_map[src]([self.source_type_map[src](geo_pt, ts) for geo_pt, ts in zip(geo_points, ts_vct)])
                source_dict[src] = vct
            results.append(source_dict)
        return results

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

        raw_fcst_lst = self._read_fcst(input_source_types)
        results = {}
        for key in list(self.qm_resolution.keys()):
            ta = self.qm_resolution[key][1]
            prep_fcst_lst, geo_points = self._prep_fcst(raw_fcst_lst, input_source_types, key)
            results[key] = self._call_qm(prep_fcst_lst, geo_points, ta, input_source_types)
        return results
