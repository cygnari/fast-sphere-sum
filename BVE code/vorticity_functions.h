#ifndef vor_H
#define vor_H

#include <iostream>
#include <cmath>
#include <array>
#include <vector>
#include <Accelerate/Accelerate.h>
#include <cassert>
#include <algorithm>
#include "helpers.h"

void rossby_haurwitz_4(vector<double>& curr_state, int point_count) {
    vector<double> curr_pos, latlon;
    double lat, lon;
    for (int i = 0; i < point_count; i++) {
        curr_pos = slice(curr_state, 5 * i, 1, 3); // gets the position of a particle
        latlon = lat_lon(curr_pos);
        lat = latlon[0];
        lon = latlon[1];
        curr_state[5 * i + 4] = lat;
        curr_state[5 * i + 3] = 2 * M_PI * sin(lat) / 7.0 + 30.0 * sin(lat) * pow(cos(lat), 4) * cos(4 * lon);
    }
}

void gauss_vortex(vector<double>& curr_state, vector<double>& area, int point_count) {
    vector<double> curr_pos, latlon, p1;
    double lat, lon, total_vor = 0.0;
    double center_lat = M_PI / 20.0;
    double center_lon = 0.0;
    for (int i = 0; i < point_count; i++) {
        curr_pos = slice(curr_state, 5 * i, 1, 3); // gets the position of a particle
        latlon = lat_lon(curr_pos);
        lat = latlon[0];
        lon = latlon[1];
        p1 = sphere_to_cart(1.0, M_PI / 2.0 - center_lat, center_lon);
        vec_minus(p1, curr_pos);
        curr_state[5 * i + 4] = lat;
        // curr_state[5 * i + 3] = 4 * M_PI * exp(-16 * pow(great_circ_dist(curr_pos, p1, 1.0), 2));
        curr_state[5 * i + 3] = 4 * M_PI * exp(-16 * pow(vec_norm(p1), 2));
    }
    for (int i = 0; i < point_count; i++) {
        total_vor += curr_state[5 * i + 3] * area[i]; //
    }
    for (int i = 0; i < point_count; i++) {
        curr_state[5 * i + 3] -= total_vor / (4.0 * M_PI);
    }
}

void vorticity_init(vector<double>& curr_state, vector<double>& area, int point_count, char* distribution) {
    if (strcmp(distribution, "rh4") == 0) {
        rossby_haurwitz_4(curr_state, point_count);
    } else if (strcmp(distribution, "gv") == 0) {
        gauss_vortex(curr_state, area, point_count);
    }
}


#endif
