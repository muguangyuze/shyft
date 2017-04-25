#pragma once
#include "core/core_pch.h"
#include "core/pt_us_k.h"

namespace shyft {
  namespace api {
    typedef shyft::core::pt_us_k::state pt_us_k_state_t;

    struct pt_us_k_state_io {
        bool from_string(const std::string &str, pt_us_k_state_t &s) const {
            return from_raw_string(str.c_str(), s);
        }

        bool from_raw_string(const char* str, pt_us_k_state_t& s) const {
            if (str && *str) {
                if (sscanf(str, "ptusk:%lf %lf %lf %lf %lf %lf %lf %lf %lf",
                    &s.us.albedo, &s.us.alpha, &s.us.sdc_melt_mean,
                    &s.us.acc_melt, &s.us.iso_pot_energy, &s.us.temp_swe, &s.us.surface_heat, &s.us.lwc,
                    &s.kirchner.q) == 9)
                    return true;

                // support old 7 string state variable format
                if (sscanf(str, "ptusk:%lf %lf %lf %lf %lf %lf %lf",
                    &s.us.albedo, &s.us.alpha, &s.us.sdc_melt_mean,
                    &s.us.acc_melt, &s.us.iso_pot_energy, &s.us.temp_swe,
                    &s.kirchner.q) == 7)
                    return true;
            }
            return false;
        }

        std::string to_string(const pt_us_k_state_t& s) const {
            char r[500];
            sprintf(r, "ptusk:%f %f %f %f %f %f %f %f %f\n",
                s.us.albedo, s.us.alpha, s.us.sdc_melt_mean,
                s.us.acc_melt, s.us.iso_pot_energy, s.us.temp_swe, s.us.surface_heat, s.us.lwc,
                s.kirchner.q);
            return r;
        }

        std::string to_string(const std::vector<pt_us_k_state_t> &sv) const {
            std::string r; r.reserve(200*200*50);
            for (size_t i = 0; i<sv.size(); ++i) {
                r.append(to_string(sv[i]));
            }
            return r;
        }

        std::vector<pt_us_k_state_t> vector_from_string(const std::string &s) const {
            std::vector<pt_us_k_state_t> r;
            if (s.size() > 0) {
                r.reserve(200*200);
                const char *l = s.c_str();
                const char *h;
                pt_us_k_state_t e;
                while (*l && (h = strstr(l, "ptusk:"))) {
                    if (!from_raw_string(h, e))
                        break;
                    r.emplace_back(e);
                    l = h + 6;// advance after ptusk marker
                }
            }
            return r;
        }
    };
  }
}
