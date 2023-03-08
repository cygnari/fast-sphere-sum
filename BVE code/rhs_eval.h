#ifndef rhs_H
#define rhs_H

#include <iostream>
#include <cmath>
#include <array>
#include <vector>
#include <Accelerate/Accelerate.h>
#include <cassert>
#include <algorithm>
#include "helpers.h"
#include "fast_rhs.h"

void rhs_direct_sum(vector<double>& modify, vector<double>& curr_state, vector<double>& area, double omega, int points) { // direct summation for all of RHS
    vector<double> pos_change, particle_i, particle_j, contribution;
    for (int i = 0; i < points; i++) {
        pos_change = {0, 0, 0};
        particle_i = slice(curr_state, 5 * i, 1, 3);
        for (int j = 0; j < points; j++) {
            if (i != j) {
                particle_j = slice(curr_state, 5 * j, 1, 3);
                contribution = BVE_gfunc(particle_i, particle_j);
                scalar_mult(contribution, curr_state[5 * j + 3] * area[j]);
                vec_add(pos_change, contribution);
            }
        }
        scalar_mult(pos_change, -1.0 / (4.0 * M_PI));
        for (int j = 0; j < 3; j++) modify[5 * i + j] = pos_change[j];
        modify[5 * i + 3] = -2 * omega * pos_change[2];
    }
}

void rhs_func(vector<double>& modify, vector<double>& curr_state, vector<double>& area, double omega, int points, char type) {
    if (type == 'd') { // d for direct
        rhs_direct_sum(modify, curr_state, area, omega, points);
    } else if (type == 'f') { // f for fast
        // rhs_fast_sum()
    }

}

#endif
