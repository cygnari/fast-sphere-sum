#ifndef fast_H
#define fast_H

#include <iostream>
#include <cmath>
#include <vector>
#include <queue>
#include <chrono>
#include <Accelerate/Accelerate.h>

#include "helpers.h"
#include "icos_utils.h"
#include "interp_utils.h"
// #include "rhs_eval.h"

struct interaction_pair {
    int lev_target; //icosahedron level of target/source
    int lev_source;
    int curr_target; // index of target/source
    int curr_source;
    int count_target;
    int count_source;
    string type; // pp, pc, cp, or cc
};

void tree_traverse(vector<interaction_pair>& interactions, icos_struct& icos1, double theta) {
    int curr_source, curr_target, lev_target, lev_source;
    int particle_count_target, particle_count_source;
}

void pp(vector<double>& modify, vector<double>& curr_state, icos_struct icos1, interaction_pair interact, double omega, vector<double>& area) {
// void direct_sum(vector<double>& modify, vector<double>& curr_state, vector<vector<vector<int>>>& tri_points, int lev_target, int lev_source, int curr_target, int curr_source, double omega, vector<double>& area) { // direct sum for P-P interaction
    // int particle_count_target, particle_count_source; // do direct summation for a pair of triangles
    int target_i, source_j;
    // particle_count_target = tri_points[lev_target][curr_target].size();
    // particle_count_source = tri_points[lev_source][curr_source].size();
    vector<double> particle_i, particle_j, contribution;
    for (int i = 0; i < interact.count_target; i++) {
        target_i = icos1.tri_points[interact.lev_target][interact.curr_target][i];
        vector<double> pos_change {0, 0, 0};
        particle_i = slice(curr_state, 5 * target_i, 1, 3);
        for (int j = 0; j < interact.count_source; j++) {
            source_j = icos1.tri_points[interact.lev_source][interact.curr_source][j];
            if (target_i != source_j) {
                particle_j = slice(curr_state, 5 * source_j, 1, 3);
                contribution = BVE_gfunc(particle_i, particle_j);
                scalar_mult(contribution, curr_state[5 * source_j + 3] * area[source_j]);
                vec_add(pos_change, contribution);
            }
        }
        for (int j = 0; j < 3; j++) {
            modify[5 * target_i + j] += pos_change[j];
        }
    }
}

void pc(vector<double>& modify, vector<double>& curr_state, icos_struct icos1, interaction_pair interact, double omega) {
    // stuff
}

#endif
