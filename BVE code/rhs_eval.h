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
#include "struct_list.h"
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

void rhs_fast_sum(vector<double>& modify, vector<double>& curr_state, vector<double>& area, vector<interaction_pair>& interactions, icos_struct icos1, interp_struct interp1, double omega, int points) {
    for (int i = 0; i < interactions.size(); i++) {
        // cout << i << endl;
        // cout << interactions[i].type << endl;
        if (interactions[i].type == "pp") pp(modify, curr_state, area, icos1, interactions[i]);
        else if (interactions[i].type == "cp") cp(modify, curr_state, area, icos1, interp1, interactions[i]);
        // else if (interactions[i].type == "pc") pc(modify, curr_state, area, icos1, interp1, interactions[i]);
        else if (interactions[i].type == "pc") pp(modify, curr_state, area, icos1, interactions[i]);
        // else cc(modify, curr_state, area, icos1, interp1, interactions[i]);
        else cp(modify, curr_state, area, icos1, interp1, interactions[i]);
    }
    scalar_mult(modify, -1.0 / (4.0 * M_PI));
    for (int i = 0; i < points; i++) modify[5 * i + 3] = -2 * omega * modify[5 * i + 2];
}

void rhs_func(vector<double>& modify, vector<double>& curr_state, vector<double>& area, double omega, run_config config1) {
    fill(modify.begin(), modify.end(), 0);
    if (config1.use_fast) {
        rhs_fast_sum(modify, curr_state, area, config1.icos.interactions, config1.icos, config1.interp, omega, config1.point_count);
    } else { // f for fast
        rhs_direct_sum(modify, curr_state, area, omega, config1.point_count);
    }

}

#endif
