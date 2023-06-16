#ifndef vorticity_H
#define vorticity_H

#include "general_utils.hpp"
#include "structs.hpp"

void rossby_haurwitz_4(run_config& run_information, vector<double>& dynamics_state) {
    vector<double> curr_pos, latlon;
    double lat, lon;
    for (int i = 0; i < run_information.dynamics_initial_points; i++) {
        curr_pos = slice(dynamics_state, run_information.info_per_point * i, 1, 3); // gets the position of a particle
        latlon = lat_lon(curr_pos);
        lat = latlon[0];
        lon = latlon[1];
        dynamics_state[run_information.info_per_point * i + 3] = 2 * M_PI * sin(lat) / 7.0 + 30.0 * sin(lat) * pow(cos(lat), 4) * cos(4 * lon);
    }
}

void gauss_vortex(run_config& run_information, vector<double>& dynamics_state, vector<double>& dynamics_areas) {
    vector<double> curr_pos, latlon, p1;
    double lat, total_vor = 0.0;
    double center_lat = M_PI / 20.0;
    double center_lon = 0.0;
    for (int i = 0; i < run_information.dynamics_initial_points; i++) {
        curr_pos = slice(dynamics_state, run_information.info_per_point * i, 1, 3); // gets the position of a particle
        latlon = lat_lon(curr_pos);
        lat = latlon[0];
        p1 = sphere_to_cart(1.0, M_PI / 2.0 - center_lat, center_lon);
        vec_minus(p1, curr_pos);
        dynamics_state[run_information.info_per_point * i + 3] = 4 * M_PI * exp(-16 * pow(vec_norm(p1), 2));
    }
    for (int i = 0; i < run_information.dynamics_initial_points; i++) {
        total_vor += dynamics_state[run_information.info_per_point * i + 3] * dynamics_areas[i]; //
    }
    for (int i = 0; i < run_information.dynamics_initial_points; i++) {
        dynamics_state[run_information.info_per_point * i + 3] -= total_vor / (4.0 * M_PI);
    }
}

void ssw_initial(run_config& run_information, vector<double>& dynamics_state, vector<double>& dynamics_areas) {
    vector<double> curr_pos, latlon, p1;
}

#endif
