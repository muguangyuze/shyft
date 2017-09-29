#pragma once
#include "core/core_pch.h"
#include "core/pt_hps_k.h"

namespace shyft {
  namespace api {
    typedef shyft::core::pt_hps_k::state_t pt_hps_k_state_t;

    struct pt_hps_k_state_io {
        bool from_string(const std::string &str, pt_hps_k_state_t &s) const {
            return from_raw_string(str.c_str(), s);
        }

        bool from_raw_string(const char* str, pt_hps_k_state_t& s) const {
            if (str && *str) {
                if (sscanf(str, "pthpsk:%lf %lf %lf %lf",
                    &s.hps.sca,&s.hps.swe, &s.hps.surface_heat,
					&s.kirchner.q) == 4)
                    return true;
            }
            return false;
        }

        std::string to_string(const pt_hps_k_state_t& s) const {
            char r[500];
            sprintf(r, "pthpsk:%f %f %f %f\n",
				s.hps.sca,s.hps.swe, s.hps.surface_heat,
				s.kirchner.q);
            return r;
        }

        std::string to_string(const std::vector<pt_hps_k_state_t> &sv) const {
            std::string r; r.reserve(200*200*50);
            for (size_t i = 0; i<sv.size(); ++i) {
                r.append(to_string(sv[i]));
            }
            return r;
        }

        std::vector<pt_hps_k_state_t> vector_from_string(const std::string &s) const {
            std::vector<pt_hps_k_state_t> r;
            if (s.size() > 0) {
                r.reserve(200*200);
                const char *l = s.c_str();
                const char *h;
                pt_hps_k_state_t e;
                while (*l && (h = strstr(l, "pthpsk:"))) {
                    if (!from_raw_string(h, e))
                        break;
                    r.emplace_back(e);
                    l = h + 7;  // advance after pthpsk marker
                }
            }
            return r;
        }
    };
  }
}
