#ifndef vorticity_H
#define vorticity_H

#include "general_utils.hpp"
#include "structs.hpp"

void rossby_haurwitz(run_config& run_information, vector<double>& dynamics_state, double omega) {
    vector<double> curr_pos, latlon;
    double lat, lon;
    int deg = run_information.init_cond_param1;
    double w = (run_information.init_cond_param2 * (1 + deg) * (2 + deg) + 2 * omega) / (deg * (3 + deg));
    for (int i = 0; i < run_information.dynamics_initial_points; i++) {
        curr_pos = slice(dynamics_state, run_information.info_per_point * i, 1, 3); // gets the position of a particle
        latlon = lat_lon(curr_pos);
        lat = latlon[0];
        lon = latlon[1];
        dynamics_state[run_information.info_per_point * i + 3] = 2 * w * sin(lat) + (pow(deg, 2) + 3 * deg + 2) * sin(lat) * pow(cos(lat), deg) * cos(deg * lon);
    }
}

void gauss_vortex(run_config& run_information, vector<double>& dynamics_state) {
    vector<double> curr_pos, latlon, p1;
    double lat, dist;
    double center_lat = run_information.init_cond_param2 * M_PI;
    double center_lon = 0.0;
    for (int i = 0; i < run_information.dynamics_initial_points; i++) {
        curr_pos = slice(dynamics_state, run_information.info_per_point * i, 1, 3); // gets the position of a particle
        latlon = lat_lon(curr_pos);
        lat = latlon[0];
        p1 = sphere_to_cart(run_information.radius, M_PI / 2.0 - center_lat, center_lon);
        vec_minus(p1, curr_pos);
        dist = vec_norm(p1);
        dynamics_state[run_information.info_per_point * i + 3] = 4 * M_PI * exp(-pow(run_information.init_cond_param1, 2) * pow(dist, 2));
    }
}

void rankine_vortex(run_config& run_information, vector<double>& dynamics_state) {
    vector<double> curr_pos, latlon, p1;
    double lat, dist, val;
    double center_lat = run_information.init_cond_param2 * M_PI;
    double center_lon = 0.0;
    double radius = run_information.init_cond_param1 * 0.01 * run_information.radius;
    for (int i = 0; i < run_information.dynamics_initial_points; i++) {
        curr_pos = slice(dynamics_state, run_information.info_per_point * i, 1, 3);
        latlon = lat_lon(curr_pos);
        lat = latlon[0];
        p1 = sphere_to_cart(run_information.radius, M_PI / 2.0 - center_lat, center_lon);
        vec_minus(p1, curr_pos);
        dist = vec_norm(p1);
        if (dist < radius) {
            val = dist / pow(radius, 2);
        } else {
            val = 1.0 / dist;
        }
        dynamics_state[run_information.info_per_point * i + 3] = 4 * M_PI * val;
    }
}

void polar_vortex(run_config& run_information, vector<double>& dynamics_state) {
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

double ssw_force(vector<double>& curr_pos, double time, double omega, int wavenumber, double time_dur) {
    vector<double> latlon;
    double Tp = 4, Tf = Tp + time_dur, lat, lon, Afunc, Bfunc, theta1 = M_PI / 3.0;
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
    return 0.6 * omega * Afunc * Bfunc * cos(wavenumber * lon) * wavenumber;
}

double ssw_blend(vector<double>& curr_pos, double time, double omega, double time_dur) {
    double Tp = 4, Tf = Tp + time_dur;
    if (time < Tp) {
        return ssw_force(curr_pos, time, omega, 1, time_dur);
    } else if (time < Tf - Tp) {
        double ratio = (time - Tp) / (Tf - Tp);
        return (1 - ratio) * ssw_force(curr_pos, time, omega, 1, time_dur) + ratio * ssw_force(curr_pos, time, omega, 2, time_dur);
    } else if (time < Tf) {
        return ssw_force(curr_pos, time, omega, 2, time_dur);
    } else {
        return 0;
    }
}

double vor_force_func(run_config& run_information, vector<double>& curr_pos, double time, double omega) {
    if (run_information.vor_forcing == "ssw") {
        return ssw_force(curr_pos, time, omega, run_information.forcing_param1, run_information.forcing_param2);
    } else if (run_information.vor_forcing == "ssw_blend") {
        return ssw_blend(curr_pos, time, omega, run_information.forcing_param2);
    } else {
        return 0;
    }
}

#endif
