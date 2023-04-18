#ifndef init_H
#define init_H

#include <cmath>
#include <vector>
#include <iostream>
#include "general_utils.hpp"
#include "structs.hpp"
#include "vorticity_functions.hpp"

void dynamics_points_initialize(run_config& run_information, vector<double>& dynamics_state,
        vector<vector<vector<int>>>& dynamics_triangles,
        // vector<vector<vector<int>>>& dynamics_points_adj_triangles,
        // vector<vector<int>>& dynamics_parent_triangles, vector<vector<int>>& dynamics_child_triangles,
        vector<vector<int>>& dynamics_points_parents, vector<vector<bool>>& dynamics_triangles_is_leaf) {
    // creates all the dynamics points and corresponding triangles
    double phi = (1 + sqrt(5)) / 2;
    int iv1, iv2, iv3, iv12, iv23, iv31;
    vector<double> v1, v2, v3, v12, v23, v31;
    dynamics_state.resize(run_information.dynamics_initial_points * run_information.info_per_point, 0);
    // dynamics_triangles.resize(run_information.dynamics_initial_triangles, vector<int> (4, 0));
    dynamics_triangles.resize(run_information.dynamics_levels_max);
    dynamics_triangles[0] = vector<vector<int>> (20, vector<int> (4, 0));
    // dynamics_points_adj_triangles.resize(run_information.dynamics_initial_points, vector<vector<int>> (run_information.dynamics_levels_max));
    // dynamics_parent_triangles.resize(run_information.dynamics_levels_max);
    // dynamics_child_triangles.resize(run_information.dynamics_initial_triangles, vector<int> (4, -1));
    dynamics_points_parents.resize(run_information.dynamics_initial_points, vector<int> (2, -1));
    dynamics_triangles_is_leaf.resize(run_information.dynamics_levels_max);
    run_information.dynamics_curr_point_count = 12;
    run_information.dynamics_curr_tri_count = 20;
    vector_copy2(dynamics_state, project_to_sphere_2(vector<double> {0, 1, phi}, run_information.radius), 0, 3);
    vector_copy2(dynamics_state, project_to_sphere_2(vector<double> {0, -1, phi}, run_information.radius), run_information.info_per_point, 3);
    vector_copy2(dynamics_state, project_to_sphere_2(vector<double> {0, 1, -phi}, run_information.radius), 2 * run_information.info_per_point, 3);
    vector_copy2(dynamics_state, project_to_sphere_2(vector<double> {0, -1, -phi}, run_information.radius), 3 * run_information.info_per_point, 3);
    vector_copy2(dynamics_state, project_to_sphere_2(vector<double> {1, phi, 0}, run_information.radius), 4 * run_information.info_per_point, 3);
    vector_copy2(dynamics_state, project_to_sphere_2(vector<double> {1, -phi, 0}, run_information.radius), 5 * run_information.info_per_point, 3);
    vector_copy2(dynamics_state, project_to_sphere_2(vector<double> {-1, phi, 0}, run_information.radius), 6 * run_information.info_per_point, 3);
    vector_copy2(dynamics_state, project_to_sphere_2(vector<double> {-1, -phi, 0}, run_information.radius), 7 * run_information.info_per_point, 3);
    vector_copy2(dynamics_state, project_to_sphere_2(vector<double> {phi, 0, 1}, run_information.radius), 8 * run_information.info_per_point, 3);
    vector_copy2(dynamics_state, project_to_sphere_2(vector<double> {phi, 0, -1}, run_information.radius), 9 * run_information.info_per_point, 3);
    vector_copy2(dynamics_state, project_to_sphere_2(vector<double> {-phi, 0, 1}, run_information.radius), 10 * run_information.info_per_point, 3);
    vector_copy2(dynamics_state, project_to_sphere_2(vector<double> {-phi, 0, -1}, run_information.radius), 11 * run_information.info_per_point, 3);
    dynamics_triangles[0][0].insert(dynamics_triangles[0][0].begin(), {0, 1, 8, 0}); // 0, 1, 2 are indices of the three vertices
    dynamics_triangles[0][1].insert(dynamics_triangles[0][1].begin(), {0, 1, 10, 0}); // 20 starting faces
    dynamics_triangles[0][2].insert(dynamics_triangles[0][2].begin(), {0, 4, 6, 0});
    dynamics_triangles[0][3].insert(dynamics_triangles[0][3].begin(), {0, 4, 8, 0});
    dynamics_triangles[0][4].insert(dynamics_triangles[0][4].begin(), {0, 6, 10, 0});
    dynamics_triangles[0][5].insert(dynamics_triangles[0][5].begin(), {1, 5, 7, 0});
    dynamics_triangles[0][6].insert(dynamics_triangles[0][6].begin(), {1, 5, 8, 0});
    dynamics_triangles[0][7].insert(dynamics_triangles[0][7].begin(), {1, 7, 10, 0});
    dynamics_triangles[0][8].insert(dynamics_triangles[0][8].begin(), {2, 3, 9, 0});
    dynamics_triangles[0][9].insert(dynamics_triangles[0][9].begin(), {2, 3, 11, 0});
    dynamics_triangles[0][10].insert(dynamics_triangles[0][10].begin(), {2, 4, 6, 0});
    dynamics_triangles[0][11].insert(dynamics_triangles[0][11].begin(), {2, 4, 9, 0});
    dynamics_triangles[0][12].insert(dynamics_triangles[0][12].begin(), {2, 6, 11, 0});
    dynamics_triangles[0][13].insert(dynamics_triangles[0][13].begin(), {3, 5, 7, 0});
    dynamics_triangles[0][14].insert(dynamics_triangles[0][14].begin(), {3, 5, 9, 0});
    dynamics_triangles[0][15].insert(dynamics_triangles[0][15].begin(), {3, 7, 11, 0});
    dynamics_triangles[0][16].insert(dynamics_triangles[0][16].begin(), {4, 8, 9, 0});
    dynamics_triangles[0][17].insert(dynamics_triangles[0][17].begin(), {5, 8, 9, 0});
    dynamics_triangles[0][18].insert(dynamics_triangles[0][18].begin(), {6, 10, 11, 0});
    dynamics_triangles[0][19].insert(dynamics_triangles[0][19].begin(), {7, 10, 11, 0});
    // dynamics_points_adj_triangles[0][0] = {0, 1, 2, 3, 4};
    // dynamics_points_adj_triangles[1][0] = {0, 1, 5, 6, 7};
    // dynamics_points_adj_triangles[2][0] = {8, 9, 10, 11, 12};
    // dynamics_points_adj_triangles[3][0] = {8, 9, 13, 14, 15};
    // dynamics_points_adj_triangles[4][0] = {2, 3, 10, 11, 16};
    // dynamics_points_adj_triangles[5][0] = {5, 6, 13, 14, 17};
    // dynamics_points_adj_triangles[6][0] = {2, 4, 10, 12, 18};
    // dynamics_points_adj_triangles[7][0] = {5, 7, 13, 15, 19};
    // dynamics_points_adj_triangles[8][0] = {0, 3, 6, 16, 17};
    // dynamics_points_adj_triangles[9][0] = {8, 11, 14, 16, 17};
    // dynamics_points_adj_triangles[10][0] = {1, 4, 7, 18, 19};
    // dynamics_points_adj_triangles[11][0] = {9, 12, 15, 18, 19};
    for (int i = 0; i < run_information.dynamics_levels_min - 1; i++) {
        dynamics_triangles[i+1] = vector<vector<int>> (20 * pow(4, i+1), vector<int> (4, 0));
        for (int j = 0; j < 20 * pow(4, i); j++) {
            iv1 = dynamics_triangles[i][j][0];
            iv2 = dynamics_triangles[i][j][1];
            iv3 = dynamics_triangles[i][j][2];
            v1 = slice(dynamics_state, run_information.info_per_point * iv1, 1, 3);
            v2 = slice(dynamics_state, run_information.info_per_point * iv2, 1, 3);
            v3 = slice(dynamics_state, run_information.info_per_point * iv3, 1, 3);
            v12 = v1;
            v23 = v2;
            v31 = v3;
            vec_add(v12, v2);
            vec_add(v23, v3);
            vec_add(v31, v1);
            scalar_mult(v12, 0.5);
            scalar_mult(v23, 0.5);
            scalar_mult(v31, 0.5);
            project_to_sphere(v12, run_information.radius);
            project_to_sphere(v23, run_information.radius);
            project_to_sphere(v31, run_information.radius);
            iv12 = check_point_exist(dynamics_points_parents, run_information.dynamics_curr_point_count, min(iv1, iv2), max(iv1, iv2));
            iv23 = check_point_exist(dynamics_points_parents, run_information.dynamics_curr_point_count, min(iv2, iv3), max(iv2, iv3));
            iv31 = check_point_exist(dynamics_points_parents, run_information.dynamics_curr_point_count, min(iv3, iv1), max(iv3, iv1));
            if (iv12 == -1) {
                iv12 = run_information.dynamics_curr_point_count;
                run_information.dynamics_curr_point_count += 1;
                dynamics_points_parents[iv12] = {min(iv1, iv2), max(iv1, iv2)};
                vector_copy2(dynamics_state, v12, iv12 * run_information.info_per_point, 3);
            }
            if (iv23 == -1) {
                iv23 = run_information.dynamics_curr_point_count;
                run_information.dynamics_curr_point_count += 1;
                dynamics_points_parents[iv23] = {min(iv2, iv3), max(iv2, iv3)};
                vector_copy2(dynamics_state, v23, iv23 * run_information.info_per_point, 3);
            }
            if (iv31 == -1) {
                iv31 = run_information.dynamics_curr_point_count;
                run_information.dynamics_curr_point_count += 1;
                dynamics_points_parents[iv31] = {min(iv3, iv1), max(iv3, iv1)};
                vector_copy2(dynamics_state, v31, iv31 * run_information.info_per_point, 3);
            }
            dynamics_triangles[i+1][4*j].insert(dynamics_triangles[i+1][4*j].begin(), {iv1, iv12, iv31, i + 1});
            dynamics_triangles[i+1][4*j+1].insert(dynamics_triangles[i+1][4*j+1].begin(), {iv2, iv12, iv23, i + 1});
            dynamics_triangles[i+1][4*j+2].insert(dynamics_triangles[i+1][4*j+2].begin(), {iv3, iv23, iv31, i + 1});
            dynamics_triangles[i+1][4*j+3].insert(dynamics_triangles[i+1][4*j+3].begin(), {iv12, iv23, iv31, i + 1});
        }
    }
    // for (int i = 0; i < run_information.dynamics_initial_triangles; i++) {
    // dynamics_triangles_is_leaf[run_information.dynamics_levels_min].resize(run_information.dynamics_initial_triangles, true);
    // }
    for (int i = 0; i < run_information.dynamics_levels_max; i++) {
        if (i == run_information.dynamics_levels_min - 1) {
            dynamics_triangles_is_leaf[i].resize(20 * pow(4, i), true);
        } else {
            dynamics_triangles_is_leaf[i].resize(20 * pow(4, i), false);
        }

    }
}

void area_initialize(run_config& run_information, vector<double>& dynamics_state, vector<vector<vector<int>>>& dynamics_triangles, vector<double>& dynamics_areas) {
    // initialize areas, node patch area for each point
    int iv1, iv2, iv3;
    vector<double> v1, v2, v3;
    double tri_area;
    for (int i = 0; i < run_information.dynamics_initial_triangles; i++) {
        iv1 = dynamics_triangles[run_information.dynamics_levels_min-1][i][0];
        iv2 = dynamics_triangles[run_information.dynamics_levels_min-1][i][1];
        iv3 = dynamics_triangles[run_information.dynamics_levels_min-1][i][2];
        v1 = slice(dynamics_state, run_information.info_per_point * iv1, 1, 3);
        v2 = slice(dynamics_state, run_information.info_per_point * iv2, 1, 3);
        v3 = slice(dynamics_state, run_information.info_per_point * iv3, 1, 3);
        tri_area = sphere_tri_area(v1, v2, v3, run_information.radius);
        dynamics_areas[iv1] += tri_area / 3.0;
        dynamics_areas[iv2] += tri_area / 3.0;
        dynamics_areas[iv3] += tri_area / 3.0;
    }
}

void area_init(run_config& run_information, vector<double>& curr_state, vector<double>& area, vector<vector<int>>& triangles, int tri_count) {
    int iv1, iv2, iv3;
    double curr_area;
    vector<double> v1, v2, v3;
    for (int i = 0; i < tri_count; i++) {
        iv1 = triangles[i][0];
        iv2 = triangles[i][1];
        iv3 = triangles[i][2];
        v1 = slice(curr_state, 5 * iv1, 1, 3);
        v2 = slice(curr_state, 5 * iv2, 1, 3);
        v3 = slice(curr_state, 5 * iv3, 1, 3);
        curr_area = sphere_tri_area(v1, v2, v3, run_information.radius);
        area[iv1] += curr_area / 3.0;
        area[iv2] += curr_area / 3.0;
        area[iv3] += curr_area / 3.0;
    }
}

void vorticity_initialize(run_config& run_information, vector<double>& dynamics_state, vector<double>& dynamics_areas) {
    // initializes the initial vorticity
    if (run_information.initial_vor_condition == "rh4") {
        rossby_haurwitz_4(run_information, dynamics_state);
    } else if (run_information.initial_vor_condition == "gv") {
        gauss_vortex(run_information, dynamics_state, dynamics_areas);
    }
}

void tracer_initialize(run_config& run_information, vector<double>& dynamics_state) {
    // initializes the tracer
    vector<double> curr_pos, latlon;
    double lat, lon;
    for (int i = 0; i < run_information.dynamics_initial_points; i++) {
        curr_pos = slice(dynamics_state, run_information.info_per_point * i, 1, 3); // gets the position of a particle
        latlon = lat_lon(curr_pos);
        lat = latlon[0];
        lon = latlon[1];
        dynamics_state[run_information.info_per_point * i + 4] = 2 + lat;
    }
}

void fixer_init(run_config& run_information, vector<double>& dynamics_state, vector<double>& dynamics_areas, vector<double>& qmins, vector<double>& qmaxs, vector<double>& target_mass, double omega) {
    double abs_vor;
    qmins.resize(run_information.tracer_count + 1, INT_MAX); // vorticity + tracers
    qmaxs.resize(run_information.tracer_count + 1, INT_MIN); // voriticty + tracers
    target_mass.resize(run_information.tracer_count + 1, 0);
    for (int i = 0; i < run_information.tracer_count; i++) {
        for (int j = 0; j < run_information.dynamics_initial_points; j++) {
            qmins[i + 1] = min(dynamics_state[run_information.info_per_point * j + i + 4], qmins[i + 1]);
            qmaxs[i + 1] = max(dynamics_state[run_information.info_per_point * j + i + 4], qmaxs[i + 1]);
            target_mass[i + 1] += dynamics_areas[j] * dynamics_state[run_information.info_per_point * j + i + 4];
        }
    }
    for (int i = 0; i < run_information.dynamics_initial_points; i++) {
        abs_vor = dynamics_state[run_information.info_per_point * i + 3] + 2 * omega * dynamics_state[run_information.info_per_point * i + 2];
        qmins[0] = min(abs_vor, qmins[0]); // min and max of absolute vorticity
        qmaxs[0] = max(abs_vor, qmaxs[0]); // total vorticity is 0
    }
}

#endif
