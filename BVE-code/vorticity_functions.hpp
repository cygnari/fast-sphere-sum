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
    vector<double> curr_pos, latlon;
    double theta0 = 15.0 * M_PI / 32.0, beta = 1.5;
    double lat, vor;
    for (int i = 0; i < run_information.dynamics_initial_points; i++) {
        curr_pos = slice(dynamics_state, run_information.info_per_point * i, 1, 3);
        latlon = lat_lon(curr_pos);
        lat = latlon[0];
        vor = 2 * cos(lat) * pow(beta, 2) * (cos(theta0) * sin(lat) - sin(theta0) * cos(lat)) + sin(lat);
        vor *= M_PI * exp(-2 * pow(beta, 2) * (1 - cos(theta0) * cos(lat) - sin(theta0 * sin(lat))));
        dynamics_state[run_information.info_per_point * i + 3] = vor;
    }
}

double ssw_force(vector<double>& curr_pos, double time, double omega) {
    vector<double> latlon;
    double Tp = 4, Tf = 15, lat, lon, Afunc, Bfunc, theta1 = M_PI / 3.0;
    latlon = lat_lon(curr_pos);
    lat = latlon[0];
    lon = latlon[1];
    if (time < Tp) {
        Afunc = 0.5 * (1 - cos(M_PI * time / Tp));
    } else if (time < Tf - Tp) {
        Afunc = 1;
    } else if (time < Tf) {
        Afunc = 0.5 * (1 - cos(M_PI * (time - Tf + Tp) / Tp + M_PI));
    } else {
        Afunc = 0;
    }
    if (lat > 0) {
        Bfunc = pow(tan(theta1), 2) / pow(tan(lat), 2) * exp(1 - pow(tan(theta1), 2) / pow(tan(lat), 2));
    }
    return 0.6 * omega * Afunc * Bfunc * cos(lon);
}

double vor_force_func(run_config& run_information, vector<double>& curr_pos, double time, double omega) {
    if (run_information.vor_forcing == "ssw") {
        return ssw_force(curr_pos, time, omega);
    } else {
        return 0;
    }
}

#endif
