#include "boostpython_pch.h"
#include <boost/python/docstring_options.hpp>
#include "core/utctime_utilities.h"
#include "core/priestley_taylor.h"
#include "core/actual_evapotranspiration.h"
#include "core/precipitation_correction.h"
#include "core/universal_snow.h"
#include "core/kirchner.h"
#include "core/pt_us_k.h"
#include "api/api.h"
#include "api/pt_us_k.h"
#include "core/pt_us_k_cell_model.h"
#include "core/region_model.h"
#include "core/model_calibration.h"
#include "expose_statistics.h"
#include "expose.h"

static char const* version() {
   return "v1.0";
}

namespace expose {
    namespace pt_us_k {
        using namespace boost::python;
        using namespace shyft::core;
        using namespace shyft::core::pt_us_k;

        static void
        parameter_state_response() {

            class_<parameter,bases<>,std::shared_ptr<parameter>>("PTUSKParameter",
                              "Contains the parameters to the methods used in the PTUSK assembly\n"
                              "priestley_taylor,universal_snow,actual_evapotranspiration,precipitation_correction,kirchner\n"
                )
                .def(init<priestley_taylor::parameter,universal_snow::parameter,actual_evapotranspiration::parameter,kirchner::parameter,precipitation_correction::parameter, optional<glacier_melt::parameter,routing::uhg_parameter>>(args("pt","us","ae","k","p_corr","gm","routing"),"create object with specified parameters"))
                .def(init<const parameter&>(args("p"),"clone a parameter"))
                .def_readwrite("pt",&parameter::pt,"priestley_taylor parameter")
                .def_readwrite("us",&parameter::us,"gamma-snow parameter")
                .def_readwrite("gm", &parameter::gm, "glacier melt parameter")
				.def_readwrite("ae",&parameter::ae,"actual evapotranspiration parameter")
                .def_readwrite("kirchner",&parameter::kirchner,"kirchner parameter")
                .def_readwrite("p_corr",&parameter::p_corr,"precipitation correction parameter")
                .def_readwrite("routing",&parameter::routing,"routing cell-to-river catchment specific parameters")
                .def("size",&parameter::size,"returns total number of calibration parameters")
                .def("set",&parameter::set,args("p"),"set parameters from vector/list of float, ordered as by get_name(i)")
                .def("get",&parameter::get,args("i"),"return the value of the i'th parameter, name given by .get_name(i)")
                .def("get_name",&parameter::get_name,args("i"),"returns the i'th parameter name, see also .get()/.set() and .size()")
                ;

            typedef std::map<int,parameter> PTUSKParameterMap;
            class_<PTUSKParameterMap>("PTUSKParameterMap","dict (int,parameter)  where the int is the catchment_id")
                .def(map_indexing_suite<PTUSKParameterMap>())
            ;

            class_<state>("PTUSKState")
                .def(init<universal_snow::state,kirchner::state>(args("us","k"),"initializes state with gamma-snow us and kirchner k"))
                .def_readwrite("us",&state::us,"gamma-snow state")
                .def_readwrite("kirchner",&state::kirchner,"kirchner state")
                ;

            typedef std::vector<state> PTUSKStateVector;
            class_<PTUSKStateVector,bases<>,std::shared_ptr<PTUSKStateVector> >("PTUSKStateVector")
                .def(vector_indexing_suite<PTUSKStateVector>())
                ;


            class_<response>("PTUSKResponse","This struct contains the responses of the methods used in the PTUSK assembly")
                .def_readwrite("pt",&response::pt,"priestley_taylor response")
                .def_readwrite("us",&response::us,"gamma-snnow response")
                .def_readwrite("gm_melt_m3s", &response::gm_melt_m3s, "glacier melt response[m3s]")
                .def_readwrite("ae",&response::ae,"actual evapotranspiration response")
                .def_readwrite("kirchner",&response::kirchner,"kirchner response")
                .def_readwrite("total_discharge",&response::total_discharge,"total stack response")
                ;
        }

        static void
        collectors() {
            typedef shyft::core::pt_us_k::all_response_collector PTUSKAllCollector;
            class_<PTUSKAllCollector>("PTUSKAllCollector", "collect all cell response from a run")
                .def_readonly("destination_area",&PTUSKAllCollector::destination_area,"a copy of cell area [m2]")
                .def_readonly("avg_discharge",&PTUSKAllCollector::avg_discharge,"Kirchner Discharge given in [m^3/s] for the timestep")
                .def_readonly("snow_sca",&PTUSKAllCollector::snow_sca," gamma snow covered area fraction, sca.. 0..1 - at the end of timestep (state)")
                .def_readonly("snow_swe",&PTUSKAllCollector::snow_swe,"gamma snow swe, [mm] over the cell sca.. area, - at the end of timestep")
                .def_readonly("snow_outflow",&PTUSKAllCollector::snow_outflow," gamma snow output [m^3/s] for the timestep")
                .def_readonly("glacier_melt", &PTUSKAllCollector::glacier_melt, " glacier melt (outflow) [m3/s] for the timestep")
                .def_readonly("ae_output",&PTUSKAllCollector::ae_output,"actual evap mm/h")
                .def_readonly("pe_output",&PTUSKAllCollector::pe_output,"pot evap mm/h")
                .def_readonly("end_reponse",&PTUSKAllCollector::end_reponse,"end_response, at the end of collected")
            ;

            typedef shyft::core::pt_us_k::discharge_collector PTUSKDischargeCollector;
            class_<PTUSKDischargeCollector>("PTUSKDischargeCollector", "collect all cell response from a run")
                .def_readonly("cell_area",&PTUSKDischargeCollector::cell_area,"a copy of cell area [m2]")
                .def_readonly("avg_discharge",&PTUSKDischargeCollector::avg_discharge,"Kirchner Discharge given in [m^3/s] for the timestep")
                .def_readonly("snow_sca",&PTUSKDischargeCollector::snow_sca," gamma snow covered area fraction, sca.. 0..1 - at the end of timestep (state)")
                .def_readonly("snow_swe",&PTUSKDischargeCollector::snow_swe,"gamma snow swe, [mm] over the cell sca.. area, - at the end of timestep")
                .def_readonly("end_reponse",&PTUSKDischargeCollector::end_response,"end_response, at the end of collected")
                .def_readwrite("collect_snow",&PTUSKDischargeCollector::collect_snow,"controls collection of snow routine")
                ;
            typedef shyft::core::pt_us_k::null_collector PTUSKNullCollector;
            class_<PTUSKNullCollector>("PTUSKNullCollector","collector that does not collect anything, useful during calibration to minimize memory&maximize speed")
                ;

            typedef shyft::core::pt_us_k::state_collector PTUSKStateCollector;
            class_<PTUSKStateCollector>("PTUSKStateCollector","collects state, if collect_state flag is set to true")
                .def_readwrite("collect_state",&PTUSKStateCollector::collect_state,"if true, collect state, otherwise ignore (and the state of time-series are undefined/zero)")
                .def_readonly("kirchner_discharge",&PTUSKStateCollector::kirchner_discharge,"Kirchner state instant Discharge given in m^3/s")
                .def_readonly("us_albedo",&PTUSKStateCollector::us_albedo,"")
                .def_readonly("us_lwc",&PTUSKStateCollector::us_lwc,"")
                .def_readonly("us_surface_heat",&PTUSKStateCollector::us_surface_heat,"")
                .def_readonly("us_alpha",&PTUSKStateCollector::us_alpha,"")
                .def_readonly("us_sdc_melt_mean",&PTUSKStateCollector::us_sdc_melt_mean,"")
                .def_readonly("us_acc_melt",&PTUSKStateCollector::us_acc_melt,"")
                .def_readonly("us_iso_pot_energy",&PTUSKStateCollector::us_iso_pot_energy,"")
                .def_readonly("us_temp_swe",&PTUSKStateCollector::us_temp_swe,"")
            ;

        }

        static void
        cells() {
              typedef shyft::core::cell<parameter, environment_t, state, state_collector, all_response_collector> PTUSKCellAll;
              typedef shyft::core::cell<parameter, environment_t, state, null_collector, discharge_collector> PTUSKCellOpt;
              expose::cell<PTUSKCellAll>("PTUSKCellAll","tbd: PTUSKCellAll doc");
              expose::cell<PTUSKCellOpt>("PTUSKCellOpt","tbd: PTUSKCellOpt doc");
              expose::statistics::universal_snow<PTUSKCellAll>("PTUSKCell");//it only gives meaning to expose the *All collect cell-type
              expose::statistics::actual_evapotranspiration<PTUSKCellAll>("PTUSKCell");
              expose::statistics::priestley_taylor<PTUSKCellAll>("PTUSKCell");
              expose::statistics::kirchner<PTUSKCellAll>("PTUSKCell");
              expose::cell_state_etc<PTUSKCellAll>("PTUSK");// just one expose of state

        }

        static void
        models() {
            typedef shyft::core::region_model<pt_us_k::cell_discharge_response_t, shyft::api::a_region_environment> PTUSKOptModel;
            typedef shyft::core::region_model<pt_us_k::cell_complete_response_t, shyft::api::a_region_environment> PTUSKModel;
            expose::model<PTUSKModel>("PTUSKModel","PTUSK");
            expose::model<PTUSKOptModel>("PTUSKOptModel","PTUSK");
            def_clone_to_similar_model<PTUSKModel, PTUSKOptModel>("create_opt_model_clone");
            def_clone_to_similar_model<PTUSKOptModel,PTUSKModel>("create_full_model_clone");
        }

        static void
        state_io() {
            expose::state_io<shyft::api::pt_us_k_state_io,shyft::core::pt_us_k::state>("PTUSKStateIo");
        }


        static void
        model_calibrator() {
            expose::model_calibrator<shyft::core::region_model<pt_us_k::cell_discharge_response_t,shyft::api::a_region_environment>>("PTUSKOptimizer");
        }
    }
}


BOOST_PYTHON_MODULE(_pt_us_k)
{

    boost::python::scope().attr("__doc__")="SHyFT python api for the pt_us_k model";
    boost::python::def("version", version);
	boost::python::docstring_options doc_options(true, true, false);// all except c++ signatures
    expose::pt_us_k::state_io();
    expose::pt_us_k::parameter_state_response();
    expose::pt_us_k::cells();
    expose::pt_us_k::models();
    expose::pt_us_k::collectors();
    expose::pt_us_k::model_calibrator();
}
