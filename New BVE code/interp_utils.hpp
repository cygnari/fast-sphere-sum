#ifndef interp_H
#define interp_H

#include <cmath>
#include <array>
#include <vector>
#include <Accelerate/Accelerate.h>
#include <cassert>
#include <algorithm>

#include "general_utils.hpp"

void tri_interp(int iv1, int iv2, int iv3, vector<double>& v1, vector<double>& v2, vector<double>& v3, vector<double>& curr_state, vector<double>& target_points, vector<double>& curr_target, int i, double omega) {
    vector<double> bary;
    bary = norm_barycoords(v1, v2, v3, curr_target);
    // interpolate relative vorticity directly
    target_points[5 * i + 3] = bary[0] * curr_state[5 * iv1 + 3] + bary[1] * curr_state[5 * iv2 + 3] + bary[2] * curr_state[5 * iv3 + 3];
    // interpolate passive tracer
    target_points[5 * i + 4] = bary[0] * curr_state[5 * iv1 + 4] + bary[1] * curr_state[5 * iv2 + 4] + bary[2] * curr_state[5 * iv3 + 4];
}

void project_points(vector<double>& curr_state, int point_count, double omega, double radius) {
    // project points to surface of sphere
    vector<double> projected;
    double delta_z;
    for (int i = 0; i < point_count; i++) {
        projected = slice(curr_state, 5 * i, 1, 3);
        project_to_sphere(projected, radius);
        delta_z = projected[2] - curr_state[5 * i + 2];
        for (int j = 0;j < 3; j++) curr_state[5 * i + j] = projected[j];
        curr_state[5 * i + 3] += -2 * omega * delta_z;
    }
}

#endif
