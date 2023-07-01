#ifndef rhs_H
#define rhs_H

#include "general_utils.hpp"
#include "fast_sum_utils.hpp"
#include "structs.hpp"
#include "vorticity_functions.hpp"
#include "green_funcs.hpp"

void rhs_direct_sum(run_config& run_information, vector<double>& modify, vector<double>& dynamics_state, vector<double>& dynamics_areas,
        double time, double omega, int nval, int point_offset, vector<double> (*func_point) (vector<double>&, vector<double>&)) { // direct summation for all of RHS
    vector<double> pos_change, particle_i, particle_j, contribution;
    double vor;
    for (int i = run_information.particle_lb; i < run_information.particle_ub; i++) {
        pos_change = {0, 0, 0};
        particle_i = slice(dynamics_state, run_information.info_per_point * i, 1, 3);
        for (int j = 0; j < run_information.dynamics_curr_point_count; j++) {
            if (i != j) {
                particle_j = slice(dynamics_state, run_information.info_per_point * j, 1, 3);
                contribution = func_point(particle_i, particle_j);
                vor = dynamics_state[run_information.info_per_point * j + 3];
                vor -= vor_force_func(run_information, particle_j, time, omega);
                scalar_mult(contribution, vor * dynamics_areas[j]);
                vec_add(pos_change, contribution);
            }
        }
        vector_copy(modify, pos_change, point_offset * i, nval);
    }
}

void rhs_fast_sum(run_config& run_information, vector<double>& modify, vector<double>& curr_state, vector<double>& area,
        vector<interaction_pair>& interactions, vector<vector<vector<int>>>& fast_sum_tree_tri_points, vector<vector<vector<int>>>& fast_sum_icos_tri_verts,
        vector<vector<double>>& fast_sum_icos_verts, double time, double omega, int nval, int point_offset, vector<double> (*func_point) (vector<double>&, vector<double>&)) {
    for (int i = 0; i < interactions.size(); i++) {
        if (i % run_information.mpi_P == run_information.mpi_ID) { // evenly split up interactions
            if (interactions[i].type == "pp") pp(run_information, modify, curr_state, area, interactions[i], fast_sum_tree_tri_points, time, omega, nval, point_offset, func_point);
            else if (interactions[i].type == "cp") cp(run_information, modify, curr_state, area, interactions[i], fast_sum_tree_tri_points, fast_sum_icos_tri_verts, fast_sum_icos_verts, time, omega, nval, point_offset, func_point);
            // else if (interactions[i].type == "pc") pc(run_information, modify, curr_state, area, interactions[i], fast_sum_tree_tri_points, fast_sum_icos_tri_verts, fast_sum_icos_verts, time, omega, nval, point_offset, func_point);
            else if (interactions[i].type == "pc") pp(run_information, modify, curr_state, area, interactions[i], fast_sum_tree_tri_points, time, omega, nval, point_offset, func_point);
            // else if (interactions[i].type == "cc") cc(run_information, modify, curr_state, area, interactions[i], fast_sum_tree_tri_points, fast_sum_icos_tri_verts, fast_sum_icos_verts, time, omega, nval, point_offset, func_point);
            else if (interactions[i].type == "cc") cp(run_information, modify, curr_state, area, interactions[i], fast_sum_tree_tri_points, fast_sum_icos_tri_verts, fast_sum_icos_verts, time, omega, nval, point_offset, func_point);
        }
    }
}

void convolve(run_config& run_information, vector<double>& modify, vector<double>& dynamics_state, vector<double>& dynamics_areas,
        double omega, vector<interaction_pair>& interactions, vector<vector<vector<int>>>& fast_sum_tree_tri_points, vector<vector<vector<int>>>& fast_sum_icos_tri_verts,
        vector<vector<double>>& fast_sum_icos_verts, double time, int nval, int point_offset, vector<double> (*func_point) (vector<double>&, vector<double>&)) {
    fill(modify.begin(), modify.end(), 0);
    if (run_information.use_fast) {
        rhs_fast_sum(run_information, modify, dynamics_state, dynamics_areas, interactions, fast_sum_tree_tri_points, fast_sum_icos_tri_verts, fast_sum_icos_verts, time, omega, nval, point_offset, func_point);
    } else {
        rhs_direct_sum(run_information, modify, dynamics_state, dynamics_areas, time, omega, nval, point_offset, func_point);
    }
}

void rhs_func(run_config& run_information, vector<double>& modify, vector<double>& dynamics_state, vector<double>& dynamics_areas,
        double omega, vector<interaction_pair>& interactions, vector<vector<vector<int>>>& fast_sum_tree_tri_points, vector<vector<vector<int>>>& fast_sum_icos_tri_verts,
        vector<vector<double>>& fast_sum_icos_verts, double time, int nval, int point_offset, vector<double> (*func_point) (vector<double>&, vector<double>&)) {
    convolve(run_information, modify, dynamics_state, dynamics_areas, omega, interactions, fast_sum_tree_tri_points, fast_sum_icos_tri_verts, fast_sum_icos_verts, time, nval, point_offset, func_point);
    for (int i = 0; i < run_information.dynamics_curr_point_count; i++) modify[run_information.info_per_point * i + 3] = -2 * omega * modify[run_information.info_per_point * i + 2];
}

void project_points(run_config& run_information, vector<double>& dynamics_state, double omega) {
    vector<double> projected;
    for (int i = 0; i < run_information.dynamics_curr_point_count; i++) {
        projected = slice(dynamics_state, run_information.info_per_point * i, 1, 3);
        project_to_sphere(projected, run_information.radius);
        for (int j = 0; j < 3; j++) dynamics_state[run_information.info_per_point * i + j] = projected[j];
    }
}

#endif
