#ifndef rhs_H
#define rhs_H

#include <iostream>
#include <cmath>
#include <array>
#include <vector>
#include <Accelerate/Accelerate.h>
#include <cassert>
#include <algorithm>

#include "general_utils.hpp"
#include "fast_sum_utils.hpp"
#include "structs.hpp"

void rhs_direct_sum(run_config& run_information, vector<double>& modify, vector<double>& dynamics_state, vector<double>& dynamics_areas) { // direct summation for all of RHS
    vector<double> pos_change, particle_i, particle_j, contribution;
    for (int i = 0; i < run_information.dynamics_curr_point_count; i++) {
        pos_change = {0, 0, 0};
        particle_i = slice(dynamics_state, run_information.info_per_point * i, 1, 3);
        for (int j = 0; j < run_information.dynamics_curr_point_count; j++) {
            if (i != j) {
                particle_j = slice(dynamics_state, run_information.info_per_point * j, 1, 3);
                contribution = BVE_gfunc(particle_i, particle_j);
                scalar_mult(contribution, dynamics_state[run_information.info_per_point * j + 3] * dynamics_areas[j]);
                vec_add(pos_change, contribution);
            }
        }
        // scalar_mult(pos_change, -1.0 / (4.0 * M_PI));
        // for (int j = 0; j < 3; j++) modify[5 * i + j] = pos_change[j];
        vector_copy(modify, pos_change, run_information.info_per_point * i, 3);
        // modify[run_information.info_per_point * i + 3] = -2 * omega * pos_change[2];
    }
}

void rhs_fast_sum(run_config& run_information, vector<double>& modify, vector<double>& curr_state, vector<double>& area,
        vector<interaction_pair>& interactions, vector<vector<vector<int>>>& fast_sum_tree_tri_points, vector<vector<vector<int>>>& fast_sum_icos_tri_verts,
        vector<vector<double>>& fast_sum_icos_verts) {
    for (int i = 0; i < interactions.size(); i++) {
        if (interactions[i].type == "pp") pp(run_information, modify, curr_state, area, interactions[i], fast_sum_tree_tri_points);
        else if (interactions[i].type == "cp") cp(run_information, modify, curr_state, area, interactions[i], fast_sum_tree_tri_points, fast_sum_icos_tri_verts, fast_sum_icos_verts);
        // else if (interactions[i].type == "pc") pc(run_information, modify, curr_state, area, interactions[i], fast_sum_tree_tri_points, fast_sum_icos_tri_verts, fast_sum_icos_verts);
        else if (interactions[i].type == "pc") pp(run_information, modify, curr_state, area, interactions[i], fast_sum_tree_tri_points);
        // else if (interactions[i].type == "cc") cc(run_information, modify, curr_state, area, interactions[i], fast_sum_tree_tri_points, fast_sum_icos_tri_verts, fast_sum_icos_verts);
        else if (interactions[i].type == "cc") cp(run_information, modify, curr_state, area, interactions[i], fast_sum_tree_tri_points, fast_sum_icos_tri_verts, fast_sum_icos_verts);
        // else cc(run_information, modify, curr_state, area, interactions[i], fast_sum_tree_tri_points, fast_sum_icos_tri_verts, fast_sum_icos_verts);
        // else cp(run_information, modify, curr_state, area, interactions[i], fast_sum_tree_tri_points, fast_sum_icos_tri_verts, fast_sum_icos_verts);
    }
    // scalar_mult(modify, -1.0 / (4.0 * M_PI));
}

void rhs_func(run_config& run_information, vector<double>& modify, vector<double>& dynamics_state, vector<double>& dynamics_areas,
        double omega, vector<interaction_pair>& interactions, vector<vector<vector<int>>>& fast_sum_tree_tri_points, vector<vector<vector<int>>>& fast_sum_icos_tri_verts,
        vector<vector<double>>& fast_sum_icos_verts) {
    fill(modify.begin(), modify.end(), 0);
    if (run_information.use_fast) {
        rhs_fast_sum(run_information, modify, dynamics_state, dynamics_areas, interactions, fast_sum_tree_tri_points, fast_sum_icos_tri_verts, fast_sum_icos_verts);
    } else { // f for fast
        rhs_direct_sum(run_information, modify, dynamics_state, dynamics_areas);
    }
    // rhs_direct_sum(run_information, modify, dynamics_state, dynamics_areas, omega);
    scalar_mult(modify, -1.0 / (4.0 * M_PI));
    for (int i = 0; i < run_information.dynamics_curr_point_count; i++) modify[run_information.info_per_point * i + 3] = -2 * omega * modify[run_information.info_per_point * i + 2];
}

void project_points(run_config& run_information, vector<double>& dynamics_state, double omega) {
    vector<double> projected;
    double delta_z;
    for (int i = 0; i < run_information.dynamics_curr_point_count; i++) {
        projected = slice(dynamics_state, run_information.info_per_point * i, 1, 3);
        project_to_sphere(projected, run_information.radius);
        delta_z = projected[2] - dynamics_state[run_information.info_per_point * i + 2];
        for (int j = 0;j < 3; j++) dynamics_state[run_information.info_per_point * i + j] = projected[j];
        // vector_copy(dynamics_state, projected, run_information.info_per_point * i, 3);
        dynamics_state[run_information.info_per_point * i + 3] += -2 * omega * delta_z;
    }
}

#endif
