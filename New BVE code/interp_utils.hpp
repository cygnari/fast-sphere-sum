#ifndef interp_H
#define interp_H

#include <cmath>
#include <array>
#include <vector>
#include <Accelerate/Accelerate.h>
#include <cassert>
#include <algorithm>
#include <iostream>

#include "general_utils.hpp"
#include "structs.hpp"

// void tri_interp(int iv1, int iv2, int iv3, vector<double>& v1, vector<double>& v2, vector<double>& v3, vector<double>& curr_state, vector<double>& target_points, vector<double>& curr_target, int i, double omega) {
//     vector<double> bary;
//     bary = normalized_barycoords(v1, v2, v3, curr_target);
//     // interpolate relative vorticity directly
//     target_points[5 * i + 3] = bary[0] * curr_state[5 * iv1 + 3] + bary[1] * curr_state[5 * iv2 + 3] + bary[2] * curr_state[5 * iv3 + 3];
//     // interpolate passive tracer
//     target_points[5 * i + 4] = bary[0] * curr_state[5 * iv1 + 4] + bary[1] * curr_state[5 * iv2 + 4] + bary[2] * curr_state[5 * iv3 + 4];
// }

void __attribute__((optnone)) fekete_init(vector<vector<double>>& points, int degree)  { // initializes fekete matrix, local
// void __attribute__((optimize(0))) fekete_init(vector<vector<double>>& points, int degree)  { // initializes fekete matrix, on GL
    double delta_x = 1.0 / degree;
    int index;
    double a, b, c;
    for (int i = 0; i < degree + 1; i++) {
        a = 1 - i * delta_x;
        for (int j = 0; j < i + 1; j++) {
            index = i * (i + 1) / 2 + j;
            c = j * delta_x;
            b = 1 - a - b;
            a = 0.5 * (1 + sin(M_PI / 2 * (2 * a - 1)));
            b = 0.5 * (1 + sin(M_PI / 2 * (2 * b - 1)));
            c = 0.5 * (1 + sin(M_PI / 2 * (2 * c - 1)));
            points[index][0] = a / (a + b + c);
            points[index][1] = b / (a + b + c);
            points[index][2] = c / (a + b + c);
        }
    }
}

void interp_mat_init(vector<double>& mat, vector<vector<double>>& points, int degree, int point_count) { // sets up matrix to interpolate with fekete points
    int index, place;
    double a, b;
    for (int i = 0; i < degree + 1; i++) {
        for (int j = 0; j < i + 1; j++) {
            index = i * (i + 1) / 2 + j;
            for (int k = 0; k < point_count; k++) {
                a = points[k][0];
                b = points[k][1];
                place = index * point_count + k;
                mat[place] = pow(a, i - j) * pow(b, j);
            }
        }
    }
}

double interp_eval(vector<double>& alphas, double s, double t, int degree) { // evaluate interpolation polynomial with coefficients alpha and barycentric point (s, t)
    double accum = 0;
    int index;
    for (int i = 0; i < degree + 1; i++) {
        for (int j = 0; j < i + 1; j++) {
            index = i * (i + 1) / 2 + j;
            accum += pow(s, i - j) * pow(t, j) * alphas[index = i * (i + 1) / 2 + j];
        }
    }
    return accum;
}

void biquadratic_interp(run_config& run_information, vector<double>& target_point, int target, int iv1, int iv2, int iv3, int iv4, int iv5, int iv6,
        vector<double>& target_state, vector<double>& dynamics_state) {

    vector<double> v1, v2, v3, v4, v5, v6, curr_alphas, bary_cords;
    vector<vector<double>> points(6, vector<double> (3, 0));
    vector<double> vorticity_values (6, 0);
    vector<double> tracer_values (6 * run_information.tracer_count, 0);
    vector<double> interp_matrix (36, 0);
    vector<int> ipiv (6, 0);
    char trans = 'N';
    int Num = 6, nrhs, info;

    v1 = slice(dynamics_state, run_information.info_per_point * iv1, 1, 3);
    v2 = slice(dynamics_state, run_information.info_per_point * iv2, 1, 3);
    v3 = slice(dynamics_state, run_information.info_per_point * iv3, 1, 3);
    v4 = slice(dynamics_state, run_information.info_per_point * iv4, 1, 3);
    v5 = slice(dynamics_state, run_information.info_per_point * iv5, 1, 3);
    v6 = slice(dynamics_state, run_information.info_per_point * iv6, 1, 3);

    points[0][0] = 1;
    points[1][1] = 1;
    points[2][2] = 1;
    points[3] = normalized_barycoords(v1, v2, v3, v4);
    points[4] = normalized_barycoords(v1, v2, v3, v5);
    points[5] = normalized_barycoords(v1, v2, v3, v6);

    vorticity_values[0] = dynamics_state[run_information.info_per_point * iv1 + 3];
    vorticity_values[1] = dynamics_state[run_information.info_per_point * iv2 + 3];
    vorticity_values[2] = dynamics_state[run_information.info_per_point * iv3 + 3];
    vorticity_values[3] = dynamics_state[run_information.info_per_point * iv4 + 3];
    vorticity_values[4] = dynamics_state[run_information.info_per_point * iv5 + 3];
    vorticity_values[5] = dynamics_state[run_information.info_per_point * iv6 + 3];

    for (int j = 0; j < run_information.tracer_count; j++) {
        tracer_values[6 * j] = dynamics_state[run_information.info_per_point * iv1 + 4 + j];
        tracer_values[6 * j + 1] = dynamics_state[run_information.info_per_point * iv2 + 4 + j];
        tracer_values[6 * j + 2] = dynamics_state[run_information.info_per_point * iv3 + 4 + j];
        tracer_values[6 * j + 3] = dynamics_state[run_information.info_per_point * iv4 + 4 + j];
        tracer_values[6 * j + 4] = dynamics_state[run_information.info_per_point * iv5 + 4 + j];
        tracer_values[6 * j + 5] = dynamics_state[run_information.info_per_point * iv6 + 4 + j];
    }

    interp_mat_init(interp_matrix, points, 2, 6);
    dgetrf_(&Num, &Num, &*interp_matrix.begin(), &Num, &*ipiv.begin(), &info);
    nrhs = 1;
    dgetrs_(&trans, &Num, &nrhs, &*interp_matrix.begin(), &Num, &*ipiv.begin(), &*vorticity_values.begin(), &Num, &info);
    // cout << "here 1 4" << endl;
    nrhs = run_information.tracer_count;
    dgetrs_(&trans, &Num, &nrhs, &*interp_matrix.begin(), &Num, &*ipiv.begin(), &*tracer_values.begin(), &Num, &info);

    bary_cords = barycoords(v1, v2, v3, target_point);
    target_state[run_information.info_per_point * target + 3] = interp_eval(vorticity_values, bary_cords[0], bary_cords[1], 2);
    for (int j = 0; j < run_information.tracer_count; j++) {
        curr_alphas = slice(tracer_values, 6 * j, 1, 6);
        target_state[run_information.info_per_point * target + 4 + j] = interp_eval(curr_alphas, bary_cords[0], bary_cords[1], 2);
    }
}

void remesh_points(run_config& run_information, vector<double>& target_points, vector<double>& dynamics_state,
        vector<vector<vector<int>>>& dynamics_triangles, vector<vector<bool>>& dynamics_triangles_is_leaf) {
    // remesh points back to regular point distribution
    vector<double> curr_target;
    // , v1, v2, v3, v4, v5, v6, bary_cords, vorticity_values (6, 0);
    // vector<vector<double>> points (6, vector<double> (3, 0));
    // vector<double> interp_matrix (36, 0);
    // vector<int> poss_tris;
    // vector<double> tracer_values (6 * (run_information.tracer_count), 0);
    // vector<double> curr_alphas;
    // double curr_vor;
    int iv1, iv2, iv3, iv4, iv5, iv6, curr_level, tri_loc, super_tri_loc; // , lb, ub;
    // bool found_leaf_tri, found_curr_level;
    // char trans = 'N';
    // int Num = 6, nrhs = 1, info;
    // vector<int> ipiv (6);
    for (int i = 0; i < run_information.dynamics_curr_point_count; i++) {

        curr_target = slice(target_points, run_information.info_per_point * i, 1, 3);

        tie(curr_level, tri_loc) = find_leaf_tri(curr_target, dynamics_state, dynamics_triangles, dynamics_triangles_is_leaf, run_information.info_per_point, run_information.dynamics_levels_max);
        super_tri_loc = floor(tri_loc / 4.0);

        iv1 = dynamics_triangles[curr_level-1][super_tri_loc][0];
        iv2 = dynamics_triangles[curr_level-1][super_tri_loc][1];
        iv3 = dynamics_triangles[curr_level-1][super_tri_loc][2];
        iv4 = dynamics_triangles[curr_level][4*super_tri_loc+3][0];
        iv5 = dynamics_triangles[curr_level][4*super_tri_loc+3][1];
        iv6 = dynamics_triangles[curr_level][4*super_tri_loc+3][2];

        biquadratic_interp(run_information, curr_target, i, iv1, iv2, iv3, iv4, iv5, iv6, target_points, dynamics_state);
    }
}

#endif
