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

void tri_interp(int iv1, int iv2, int iv3, vector<double>& v1, vector<double>& v2, vector<double>& v3, vector<double>& curr_state, vector<double>& target_points, vector<double>& curr_target, int i, double omega) {
    vector<double> bary;
    bary = normalized_barycoords(v1, v2, v3, curr_target);
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

void remesh_points(run_config& run_information, vector<double>& target_points, vector<double>& dynamics_state, vector<vector<vector<int>>>& dynamics_triangles, vector<vector<bool>>& dynamics_triangles_is_leaf) {
    // remesh points back to regular point distribution
    vector<double> curr_target, v1, v2, v3, v4, v5, v6, bary_cords, vorticity_values (6, 0);
    vector<vector<double>> points (6, vector<double> (3, 0));
    vector<double> interp_matrix (36, 0);
    vector<int> poss_tris;
    vector<double> tracer_values (6 * (run_information.info_per_point - 4), 0);
    vector<double> curr_alphas;
    // double curr_vor;
    int iv1, iv2, iv3, iv4, iv5, iv6, curr_level, tri_loc, super_tri_loc, lb, ub;
    bool found_leaf_tri, found_curr_level;
    char trans = 'N';
    int Num = 6, nrhs = 1, info;
    // int ipiv[12];
    vector<int> ipiv (6);
    for (int i = 0; i < run_information.dynamics_curr_point_count; i++) {
    // for (int i = 0; i < 1; i++) {
        // find vorticity at each target point
        found_leaf_tri = false;
        curr_target = slice(target_points, run_information.info_per_point * i, 1, 3);
        curr_level = 0;
        lb = 0;
        ub = 20;
        // cout << i << endl;
        // cout << curr_target[0] << " " << curr_target[1] << " " << curr_target[2] << endl;
        // found_curr_level = false;
        // while (not found_leaf_tri) {
        for (int level = 0; level < run_information.dynamics_levels_max; level++) {
        // for (int level = 0; level < 3; level++) {
            // cout << "level: " << level << endl;
            // cout << "lb: " << lb << " ub: " << ub << endl;
            found_curr_level = false;

            for (int j = lb; j < ub; j++) {
                // cout << j << endl;
                iv1 = dynamics_triangles[level][j][0];
                iv2 = dynamics_triangles[level][j][1];
                iv3 = dynamics_triangles[level][j][2];
                v1 = slice(dynamics_state, run_information.info_per_point * iv1, 1, 3);
                v2 = slice(dynamics_state, run_information.info_per_point * iv2, 1, 3);
                v3 = slice(dynamics_state, run_information.info_per_point * iv3, 1, 3);
                // cout << "v1: " << v1[0] << "," << v1[1] << "," << v1[2] << endl;
                // cout << "v2: " << v2[0] << "," << v2[1] << "," << v2[2] << endl;
                // cout << "v3: " << v3[0] << "," << v3[1] << "," << v3[2] << endl;
                bary_cords = barycoords(v1, v2, v3, curr_target);
                // cout << bary_cords[0] << " " << bary_cords[1] << " " << bary_cords[2] << endl;
                if (check_in_tri(v1, v2, v3, curr_target)) {
                    // cout << "here" << endl;
                    found_curr_level = true;
                    if (dynamics_triangles_is_leaf[level][j]) {
                        found_leaf_tri = true;
                        tri_loc = j;
                        break;
                    } else {
                        curr_level += 1;
                        lb = 4 * j;
                        ub = 4 * j + 4;
                        break;
                    }
                }
            }
            // cout << "here 4" << endl;
            if (found_leaf_tri) break;
            // cout << "here 5" << endl;
            // if (found_curr_level) break;
            if (not found_curr_level) {
                // cout << "here 3" << endl;
                for (int j = 0; j < 20 * pow(4, level); j++) {
                    // cout << j << endl;
                    iv1 = dynamics_triangles[level][j][0];
                    iv2 = dynamics_triangles[level][j][1];
                    iv3 = dynamics_triangles[level][j][2];
                    v1 = slice(dynamics_state, run_information.info_per_point * iv1, 1, 3);
                    v2 = slice(dynamics_state, run_information.info_per_point * iv2, 1, 3);
                    v3 = slice(dynamics_state, run_information.info_per_point * iv3, 1, 3);
                    // cout << "v1: " << v1[0] << "," << v1[1] << "," << v1[2] << endl;
                    // cout << "v2: " << v2[0] << "," << v2[1] << "," << v2[2] << endl;
                    // cout << "v3: " << v3[0] << "," << v3[1] << "," << v3[2] << endl;
                    bary_cords = barycoords(v1, v2, v3, curr_target);
                    // cout << bary_cords[0] << " " << bary_cords[1] << " " << bary_cords[2] << endl;
                    if (check_in_tri(v1, v2, v3, curr_target)) {
                        // cout << "here" << endl;
                        found_curr_level = true;
                        if (dynamics_triangles_is_leaf[level][j]) {
                            found_leaf_tri = true;
                            tri_loc = j;
                            break;
                        } else {
                            curr_level += 1;
                            lb = 4 * j;
                            ub = 4 * j + 4;
                            break;
                        }
                    }
                }
            }
            if (found_leaf_tri) break;
            // if ()
            // cout << i << " bad" << endl;
            // return;
        }
        super_tri_loc = floor(tri_loc / 4.0);
        // cout << tri_loc << " " << super_tri_loc << endl;
        iv1 = dynamics_triangles[curr_level-1][super_tri_loc][0];
        iv2 = dynamics_triangles[curr_level-1][super_tri_loc][1];
        iv3 = dynamics_triangles[curr_level-1][super_tri_loc][2];
        iv4 = dynamics_triangles[curr_level][4*super_tri_loc+3][0];
        iv5 = dynamics_triangles[curr_level][4*super_tri_loc+3][1];
        iv6 = dynamics_triangles[curr_level][4*super_tri_loc+3][2];

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
        for (int j = 0; j < run_information.info_per_point - 4; j++) {
            // for (int k = 0; k < 6; k++) {
            //     tracer_values[6 * j + k] = dynamics_state[]
            // }
            tracer_values[6 * j] = dynamics_state[run_information.info_per_point * iv1 + 4 + j];
            tracer_values[6 * j + 1] = dynamics_state[run_information.info_per_point * iv2 + 4 + j];
            tracer_values[6 * j + 2] = dynamics_state[run_information.info_per_point * iv3 + 4 + j];
            tracer_values[6 * j + 3] = dynamics_state[run_information.info_per_point * iv4 + 4 + j];
            tracer_values[6 * j + 4] = dynamics_state[run_information.info_per_point * iv5 + 4 + j];
            tracer_values[6 * j + 5] = dynamics_state[run_information.info_per_point * iv6 + 4 + j];
        }
        // cout << "here" << endl;
        interp_mat_init(interp_matrix, points, 2, 6);
        // cout << "here2" << endl;
        // for (int i = 0; i < interp_matrix.size(); i++) cout << interp_matrix[i] << " ";
        // for (int j = 0; j < 6; j++) {
        //     for (int k = 0; k < 6; k++) {
        //         cout << interp_matrix[j+6*k] << " ";
        //     }
        //     cout << endl;
        // }
        // cout << endl;
        // cout << interp_matrix.size() << " " << vorticity_values.size() << endl;
        // dgetrs_(&trans, &Num, &nrhs, &*interp_matrix.begin(), &Num, &*ipiv, &*vorticity_values.begin(), &Num, &info);
        dgetrf_(&Num, &Num, &*interp_matrix.begin(), &Num, &*ipiv.begin(), &info);
        dgetrs_(&trans, &Num, &nrhs, &*interp_matrix.begin(), &Num, &*ipiv.begin(), &*vorticity_values.begin(), &Num, &info);
        nrhs = run_information.info_per_point - 4;
        dgetrs_(&trans, &Num, &nrhs, &*interp_matrix.begin(), &Num, &*ipiv.begin(), &*tracer_values.begin(), &Num, &info);
        // dgesv_(&Num, &nrhs, &*interp_matrix.begin(), &Num, &*ipiv.begin(), &*vorticity_values.begin(), &Num, &info);
        bary_cords = barycoords(v1, v2, v3, curr_target);
        target_points[run_information.info_per_point * i + 3] = interp_eval(vorticity_values, bary_cords[0], bary_cords[1], 2);
        for (int j = 0; j < run_information.info_per_point - 4; j++) {
            curr_alphas = slice(tracer_values, 6 * j, 1, 6);
            target_points[run_information.info_per_point * i + 4 + j] = interp_eval(curr_alphas, bary_cords[0], bary_cords[1], 2);
        }
    }
}

#endif
