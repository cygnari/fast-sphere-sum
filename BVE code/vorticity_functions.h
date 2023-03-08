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
    double lat, lon, vor;
    for (int i = 0; i < point_count; i++) {
        curr_pos = slice(curr_state, 5 * i, 1, 3); // gets the position of a particle
        latlon = lat_lon(curr_pos);
        lat = latlon[0];
        lon = latlon[1];
        curr_state[5 * i + 4] = lat;
        curr_state[5 * i + 3] = 2 * M_PI * sin(lat) / 7.0 + 30.0 * sin(lat) * pow(cos(lat), 4) * cos(4 * lon);
    }
}

void vorticity_init(vector<double>& curr_state, int point_count, char* distribution) {
    if (strcmp(distribution, "rh4") == 0) {
        rossby_haurwitz_4(curr_state, point_count);
    }
}


#endif
