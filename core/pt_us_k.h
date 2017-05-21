#pragma once

#include "priestley_taylor.h"
#include "kirchner.h"
#include "universal_snow.h"
#include "actual_evapotranspiration.h"
#include "precipitation_correction.h"
#include "glacier_melt.h"
#include "unit_conversion.h"
#include "routing.h"
namespace shyft {
  namespace core {
    namespace pt_us_k {
        using namespace std;

        /** \brief Simple parameter struct for the PTUSK method stack
         *
         * This struct contains the parameters to the methods used in the PTUSK assembly.
         *
         * \tparam PTParameter PriestleyTaylor parameter type that implements the interface:
         *    - PTParameter.albedo const --> double, land albedo parameter in PriestleyTaylor.
         *    - PTParameter.alpha const --> double, alpha parameter in PriestleyTaylor.
         * \tparam GSState UniversalSnow parameter type that implements the parameter interface for UniversalSnow.
         * \tparam KState Kirchner parameter type that implements the parameter interface for Kirchner.
         * \sa UniversalSnowParameter \sa Kirchner \sa PTUSK  \sa PriestleyTaylor
         */
        struct parameter {
            typedef priestley_taylor::parameter pt_parameter_t;
            typedef universal_snow::parameter us_parameter_t;
            typedef actual_evapotranspiration::parameter ae_parameter_t;
            typedef kirchner::parameter kirchner_parameter_t;
            typedef precipitation_correction::parameter precipitation_correction_parameter_t;
            typedef glacier_melt::parameter glacier_melt_parameter_t;
            typedef routing::uhg_parameter routing_parameter_t;
            parameter(pt_parameter_t pt,
                      us_parameter_t us,
                      ae_parameter_t ae,
                      kirchner_parameter_t k,
                      precipitation_correction_parameter_t p_corr,
                      glacier_melt_parameter_t gm=glacier_melt_parameter_t(),
                      routing_parameter_t routing=routing_parameter_t())
             : pt(pt), us(us), ae(ae), kirchner(k), p_corr(p_corr) ,gm(gm),routing(routing){ /*Do nothing */ }
			parameter()=default;
			parameter(const parameter&)=default;
			parameter(parameter&&)=default;
			parameter& operator=(const parameter &c)=default;
			parameter& operator=(parameter&&c)=default;


            pt_parameter_t pt;
            us_parameter_t us;
            ae_parameter_t ae;
            kirchner_parameter_t  kirchner;
            precipitation_correction_parameter_t p_corr;
            glacier_melt_parameter_t gm;
            routing_parameter_t routing;
            ///<calibration support, needs vector interface to params, size is the total count
            size_t size() const { return 28; }
            ///<calibration support, need to set values from ordered vector
            void set(const vector<double>& p) {
                if (p.size() != size())
                    throw runtime_error("PTUSK Parameter Accessor: .set size missmatch");
                int i = 0;
                kirchner.c1 = p[i++];
                kirchner.c2 = p[i++];
                kirchner.c3 = p[i++];
                ae.ae_scale_factor = p[i++];
                us.tx = p[i++];
                us.wind_scale = p[i++];
                us.max_water = p[i++];
                us.wind_const = p[i++];
                us.fast_albedo_decay_rate = p[i++];
                us.slow_albedo_decay_rate = p[i++];
                us.surface_magnitude = p[i++];
                us.max_albedo = p[i++];
                us.min_albedo = p[i++];
                us.snowfall_reset_depth = p[i++];
                us.snow_cv = p[i++];
                us.glacier_albedo = p[i++];
                p_corr.scale_factor = p[i++];
                us.snow_cv_forest_factor=p[i++];
                us.snow_cv_altitude_factor=p[i++];
				pt.albedo = p[i++];
				pt.alpha = p[i++];
				us.initial_bare_ground_fraction = p[i++];
				us.winter_end_day_of_year = size_t(p[i++]);
				us.calculate_iso_pot_energy = p[i++] != 0.0 ? true : false;
                gm.dtf = p[i++];
                routing.velocity = p[i++];
                routing.alpha = p[i++];
                routing.beta  = p[i++];
            }

            ///< calibration support, get the value of i'th parameter
            double get(size_t i) const {
                switch (i) {
                    case  0:return kirchner.c1;
                    case  1:return kirchner.c2;
                    case  2:return kirchner.c3;
                    case  3:return ae.ae_scale_factor;
                    case  4:return us.tx;
                    case  5:return us.wind_scale;
                    case  6:return us.max_water;
                    case  7:return us.wind_const;
                    case  8:return us.fast_albedo_decay_rate;
                    case  9:return us.slow_albedo_decay_rate;
                    case 10:return us.surface_magnitude;
                    case 11:return us.max_albedo;
                    case 12:return us.min_albedo;
                    case 13:return us.snowfall_reset_depth;
                    case 14:return us.snow_cv;
                    case 15:return us.glacier_albedo;
                    case 16:return p_corr.scale_factor;
                    case 17:return us.snow_cv_forest_factor;
                    case 18:return us.snow_cv_altitude_factor;
					case 19:return pt.albedo;
					case 20:return pt.alpha;
					case 21:return us.initial_bare_ground_fraction;
					case 22:return (double)us.winter_end_day_of_year;
					case 23:return us.calculate_iso_pot_energy ? 1.0 : 0.0;
                    case 24:return gm.dtf;
                    case 25:return routing.velocity;
                    case 26:return routing.alpha;
                    case 27:return routing.beta;

                default:
                    throw runtime_error("PTUSK Parameter Accessor:.get(i) Out of range.");
                }
                return 0;
            }

            ///< calibration and python support, get the i'th parameter name
            string get_name(size_t i) const {
                static const char *names[] = {
                    "kirchner.c1",
                    "kirchner.c2",
                    "kirchner.c3",
                    "ae.ae_scale_factor",
                    "us.tx",
                    "us.wind_scale",
                    "us.max_water",
                    "us.wind_const",
                    "us.fast_albedo_decay_rate",
                    "us.slow_albedo_decay_rate",
                    "us.surface_magnitude",
                    "us.max_albedo",
                    "us.min_albedo",
                    "us.snowfall_reset_depth",
                    "us.snow_cv",
                    "us.glacier_albedo",
                    "p_corr.scale_factor",
                    "us.snow_cv_forest_factor",
                    "us.snow_cv_altitude_factor",
                    "pt.albedo",
                    "pt.alpha",
					"us.initial_bare_ground_fraction",
					"us.winter_end_day_of_year",
					"us.calculate_iso_pot_energy",
                    "gm.dtf",
                    "routing.velocity",
                    "routing.alpha",
                    "routing.beta"
                };
                if (i >= size())
                    throw runtime_error("PTUSK Parameter Accessor:.get_name(i) Out of range.");
                return names[i];
            }
        };

        /** \brief Simple state struct for the PTUSK method stack
         *
         * This struct contains the states of the methods used in the PTUSK assembly.
         *
         * \tparam GSState UniversalSnow state type that implements the state interface for UniversalSnow.
         * \tparam KState Kirchner state type that implements the state interface for Kirchner.
         * \sa UniversalSnowState \sa Kirchner \sa PTUSK \sa PriestleyTaylor \sa UniversalSnow
         */
        struct state {
            typedef universal_snow::state us_state_t;
            typedef kirchner::state kirchner_state_t;
            state() {}
            state(us_state_t us, kirchner_state_t k) : us(us), kirchner(k) {}
            us_state_t us;
            kirchner_state_t kirchner;
            bool operator==(const state& x) const {return us==x.us && kirchner==x.kirchner;}
            x_serialize_decl();
        };


        /** \brief Simple response struct for the PTUSK method stack
         *
         * This struct contains the responses of the methods used in the PTUSK assembly.
         */
        struct response {
            // Model responses
            typedef priestley_taylor::response  pt_response_t;
            typedef universal_snow::response us_response_t;
            typedef actual_evapotranspiration::response  ae_response_t;
            typedef kirchner::response kirchner_response_t;
            pt_response_t pt;
            us_response_t us;
            ae_response_t ae;
            kirchner_response_t kirchner;
            double gm_melt_m3s;
            // Stack response
            double total_discharge;
            double charge_m3s;
        };

        /** \brief Calculation Model using assembly of PriestleyTaylor, UniversalSnow and Kirchner
         *
         * This model first uses PriestleyTaylor for calculating the potential
         * evapotranspiration based on time series data for temperature,
         * radiation and relative humidity. Then it uses the UniversalSnow method
         * to calculate the snow/ice adjusted runoff using time series data for
         * precipitation and wind speed in addition to the time series used in
         * the PriestleyTaylor method. The PriestleyTaylor potential evaporation is
         * is used to calculate the actual evapotranspiration that is passed on to the
         * last step, Kirchner.
         * Kirchner is run with the gamma snow output
         * and actual evaporation response from the two methods above to
         * calculate the discharge.
         *
         * TODO: This method should construct an internal time stepping PTUSK struct.
         * This stepper should be used as an iterator in time integration loops. This
         * paves the way for inter-cell communication (e.g. snow transport) without
         * touching this simple interface.
         *
         * \tparam TS Time series type that implements:
         *    - TS::source_accessor_type --> Type of accessor used to retrieve time series data.
         *    - TS.accessor(const T& time_axis) const --> TS::source_accessor_type. Accessor object
         *      for this time series.
         *
         * \tparam TS::source_accessor_type Time series source accessor type that implements:
         *    - TS::source_accessor_type(const TS& time_series, const T& time_axis) --> Construct
         *      accessor for the given time series and time axis.
         *    - TS::source_accessor_type.value(size_t i) --> double, -value of the time series for
         *      period i in the time axis.
         * \tparam T Time axis type that implements:
         *    - T.size() const --> Number of periods in the time axis.
         *    - T(size_t i) const --> shyft::core::utcperiod, - Start and end as shyft::core::utctime
         *      of the i'th period.
         * \tparam S State type that implements:
         *    - S::us_state_type --> State type for the UniversalSnow method.
         *    - S::kirchner_state_type --> State type for the Kirchner method.
         *    - S.gs --> S::us_state_type, - State variables for the UniversalSnow method
         *    - S.kirchner --> S::kirchner_state_type, - State variables for the Kirchner method
         * \tparam R Response type that implements:
         *    - R::us_response_type --> Response type for the UniversalSnow routine.
         *    - R.gs --> R::us_response_type, -Response object passed to the UniversalSnow routine.
         * \tparam P Parameter type that implements:
         *    - P::pt_parameter_type --> Parameter type for the PriestleyTaylor method
         *    - P::us_parameter_type --> Parameter type for the UniversalSnow method.
         *    - P::ae_parameter_type --> Parameter type for the ActualEvapotranspiration method.
         *    - P::kirchner_parameter_type --> Parameter type for the Kirchner method.
         *    - P.pt --> P::pt_parameter_type --> Parameters for the PriestleyTaylor method.
         *    - P.gs --> P::us_parameter_type --> Parameters for thge UniversalSnow method.
         *    - P.ae --> P::ae_parameter_type --> Parameters for thge ActualEvapotranspiration method.
         *    - P.kirchner --> P::kirchner_parameter_type --> Parameters for the Kirchner method.
         * \tparam SC State collector type that implements:
         *    - SC.collect(utctime t, const S& state) --> Possibly save some states at time t.
         * \tparam RC Response collector type that implements:
         *    - RC.collect(utctime t, const R& response) --> Possibly save some responses at time t.
         */
        template<template <typename, typename> class A, class R, class T_TS, class P_TS, class WS_TS, class RH_TS, class RAD_TS, class T,
        class S, class GCD, class P, class SC, class RC>
        void run_pt_us_k(const GCD& geo_cell_data,
            const P& parameter,
            const T& time_axis, int start_step, int  n_steps,
            const T_TS& temp,
            const P_TS& prec,
            const WS_TS& wind_speed,
            const RH_TS& rel_hum,
            const RAD_TS& rad,
            S& state,
            SC& state_collector,
            RC& response_collector
            ) {
            // Access time series input data through accessors of template A (typically a direct accessor).
            using temp_accessor_t = A<T_TS, T>;
            using prec_accessor_t = A<P_TS, T>;
            using wind_speed_accessor_t = A<WS_TS, T>;
            using rel_hum_accessor_t = A<RH_TS, T>;
            using rad_accessor_t = A<RAD_TS, T>;

            auto temp_accessor = temp_accessor_t(temp, time_axis);
            auto prec_accessor = prec_accessor_t(prec, time_axis);
            auto wind_speed_accessor = wind_speed_accessor_t(wind_speed, time_axis);
            auto rel_hum_accessor = rel_hum_accessor_t(rel_hum, time_axis);
            auto rad_accessor = rad_accessor_t(rad, time_axis);

            // Initialize the method stack
            precipitation_correction::calculator p_corr(parameter.p_corr.scale_factor);
            priestley_taylor::calculator pt(parameter.pt.albedo, parameter.pt.alpha);
            universal_snow::calculator<typename P::us_parameter_t, typename S::us_state_t, typename R::us_response_t> us;
            kirchner::calculator<kirchner::trapezoidal_average, typename P::kirchner_parameter_t> kirchner(parameter.kirchner);
            //
            // Get the initial states
            R response;
            const double forest_fraction=geo_cell_data.land_type_fractions_info().forest();
            const double total_lake_fraction = geo_cell_data.land_type_fractions_info().lake() + geo_cell_data.land_type_fractions_info().reservoir(); // both give direct response for now
            const double glacier_fraction = geo_cell_data.land_type_fractions_info().glacier();
            const double kirchner_fraction = 1 - glacier_fraction;
            const double cell_area_m2 = geo_cell_data.area();
            const double glacier_area_m2 = geo_cell_data.area()*glacier_fraction;
            const double altitude= geo_cell_data.mid_point().z;
            // Step through times in axis
            size_t i_begin = n_steps > 0 ? start_step : 0;
            size_t i_end = n_steps > 0 ? start_step + n_steps : time_axis.size();
            for (size_t i = i_begin ; i < i_end ; ++i) {
                utcperiod period = time_axis.period(i);
                double temp = temp_accessor.value(i);
                double rad = rad_accessor.value(i);
                double rel_hum = rel_hum_accessor.value(i);
                double prec = p_corr.calc(prec_accessor.value(i));
                state_collector.collect(i, state);///< \note collect the state at the beginning of each period (the end state is saved anyway)

                us.step(state.us, response.us, period.start, period.timespan(), parameter.us,
                        temp, rad, prec, wind_speed_accessor.value(i), rel_hum,forest_fraction,altitude);
                response.gm_melt_m3s = glacier_melt::step(parameter.gm.dtf, temp, geo_cell_data.area()*response.us.sca, glacier_area_m2);
                response.pt.pot_evapotranspiration = pt.potential_evapotranspiration(temp, rad, rel_hum)*calendar::HOUR; //mm/s -> mm/h
                response.ae.ae = actual_evapotranspiration::calculate_step(state.kirchner.q, response.pt.pot_evapotranspiration,
                                  parameter.ae.ae_scale_factor, std::max(response.us.sca,glacier_fraction), // a evap only on non-snow/non-glac area
                                  period.timespan());
                kirchner.step(period.start, period.end, state.kirchner.q, response.kirchner.q_avg, response.us.outflow, response.ae.ae); // all units mm/h over 'same' area

                double bare_lake_fraction = total_lake_fraction*(1.0 - response.us.sca);// only direct response on bare (no snow-cover) lakes
                response.total_discharge =
                      std::max(0.0,prec - response.ae.ae)*bare_lake_fraction // when it rains, remove ae. from direct response
                    + shyft::m3s_to_mmh(response.gm_melt_m3s,cell_area_m2)
                    + response.kirchner.q_avg*(kirchner_fraction-bare_lake_fraction);//let kirchner respond to all except glacier and direct lake
                response.charge_m3s =
                    + shyft::mmh_to_m3s(prec, cell_area_m2)
                    - shyft::mmh_to_m3s(response.ae.ae, cell_area_m2)
                    + response.gm_melt_m3s
                    - shyft::mmh_to_m3s(response.total_discharge, cell_area_m2);
                // Possibly save the calculated values using the collector callbacks.
                response_collector.collect(i, response);///< \note collect the response valid for the i'th period (current state is now at the end of period)
                if(i+1==i_end)
                    state_collector.collect(i+1, state);///< \note last iteration,collect the  final state as well.
            }
            response_collector.set_end_response(response);
        }
    } // pt_us_k
  } // core
} // shyft
  //-- serialization support shyft
x_serialize_export_key(shyft::core::pt_us_k::state);