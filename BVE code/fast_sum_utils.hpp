#ifndef fast_H
#define fast_H

#include <cmath>
#include <vector>
#include <queue>
#include <chrono>
#include <Accelerate/Accelerate.h>
#include "general_utils.hpp"
#include "structs.hpp"

void point_assign(run_config& run_information, vector<double>& point, vector<vector<double>>& fast_sum_icos_verts,
        vector<vector<vector<int>>>& fast_sum_icos_tri_verts, vector<vector<vector<int>>>& fast_sum_tree_tri_points,
        vector<vector<int>>& fast_sum_tree_point_locs, int point_id) {
    int iv1, iv2, iv3, lb, ub;
    vector<double> v1, v2, v3;

    for (int i = 0; i < run_information.fast_sum_tree_levels; i++) {
        if (i > 0) {
            lb = 4 * fast_sum_tree_point_locs[i-1][point_id]; // utilize tree structure to minimize searching
            ub = lb + 4;
        } else {
            lb = 0;
            ub = 20;
        }
        for (int j = lb; j < ub; j++) {
            iv1 = fast_sum_icos_tri_verts[i][j][0];
            iv2 = fast_sum_icos_tri_verts[i][j][1];
            iv3 = fast_sum_icos_tri_verts[i][j][2];
            v1 = fast_sum_icos_verts[iv1];
            v2 = fast_sum_icos_verts[iv2];
            v3 = fast_sum_icos_verts[iv3];
            if (check_in_tri(v1, v2, v3, point)) {
                fast_sum_tree_point_locs[i][point_id] = j;
                fast_sum_tree_tri_points[i][j].push_back(point_id);
                break;
            }
        }
    }
}

void points_assign(run_config& run_information, vector<double>& dynamics_state, vector<vector<double>>& fast_sum_icos_verts,
        vector<vector<vector<int>>>& fast_sum_icos_tri_verts, vector<vector<vector<int>>>& fast_sum_tree_tri_points,
        vector<vector<int>>& fast_sum_tree_point_locs) {
    vector<double> point;
    for (int i = 0; i < run_information.fast_sum_tree_levels; i++) {
        fast_sum_tree_tri_points[i] = vector<vector<int>> (20 * pow(4, i));
        fast_sum_tree_point_locs[i] = vector<int> (run_information.dynamics_curr_point_count, 0);
    }
    for (int i = 0; i < run_information.dynamics_curr_point_count; i++) {
        point = slice(dynamics_state, run_information.info_per_point * i, 1, 3);
        point_assign(run_information, point, fast_sum_icos_verts, fast_sum_icos_tri_verts, fast_sum_tree_tri_points, fast_sum_tree_point_locs, i);
    }
}

void tree_traverse(run_config& run_information, vector<vector<vector<int>>>& fast_sum_tree_tri_points,
        vector<vector<vector<double>>>& fast_sum_icos_tri_info, vector<interaction_pair>& tree_interactions) {
    int curr_source, curr_target, lev_target, lev_source;
    int particle_count_target, particle_count_source;
    vector<double> center_target, center_source;
    double separation, distance;
    vector<vector<int>> tri_interactions;
    vector<int> curr_interact (4, 0);

    for (int i = 0; i < 20; i++) { // queue of triangle pairs to interact
        for (int j = 0; j < 20; j++) {
            tri_interactions.push_back({i, j, 0, 0});
        }
    }

    tree_interactions.clear();

    while (tri_interactions.size() > 0) {
        curr_interact = tri_interactions.front(); // get triangle pair to interact
        curr_target = curr_interact[0];
        curr_source = curr_interact[1];
        lev_target = curr_interact[2];
        lev_source = curr_interact[3];
        tri_interactions.erase(tri_interactions.begin());
        particle_count_target = fast_sum_tree_tri_points[lev_target][curr_target].size();
        particle_count_source = fast_sum_tree_tri_points[lev_source][curr_source].size();
        if ((particle_count_target == 0) or (particle_count_source == 0)) continue; // if no work, continue to next
        center_target = slice(fast_sum_icos_tri_info[lev_target][curr_target], 0, 1, 3);
        center_source = slice(fast_sum_icos_tri_info[lev_source][curr_source], 0, 1, 3);
        distance = great_circ_dist(center_target, center_source, run_information.radius);
        separation = (fast_sum_icos_tri_info[lev_target][curr_target][3] + fast_sum_icos_tri_info[lev_source][curr_source][3]) / distance;

        if ((distance > 0) and (separation < run_information.fast_sum_theta)) { // triangles are well separated
            interaction_pair new_interact = {lev_target, lev_source, curr_target, curr_source, particle_count_target, particle_count_source, ""};
            if (particle_count_target > run_information.fast_sum_cluster_thresh) {
                new_interact.type += "c";
            } else {
                new_interact.type += "p";
            }
            if (particle_count_source > run_information.fast_sum_cluster_thresh) {
                new_interact.type += "c";
            } else {
                new_interact.type += "p";
            }
            tree_interactions.push_back(new_interact);
        } else {
            if ((particle_count_target < run_information.fast_sum_cluster_thresh) and (particle_count_source < run_information.fast_sum_cluster_thresh)) { // both have few particles, pp
                interaction_pair new_interact = {lev_target, lev_source, curr_target, curr_source, particle_count_target, particle_count_source, "pp"};
                tree_interactions.push_back(new_interact);
            } else if ((lev_target == run_information.fast_sum_tree_levels - 1) and (lev_source == run_information.fast_sum_tree_levels - 1)) { // both are leaves, pp
                interaction_pair new_interact = {lev_target, lev_source, curr_target, curr_source, particle_count_target, particle_count_source, "pp"};
                tree_interactions.push_back(new_interact);
            } else if (lev_target == run_information.fast_sum_tree_levels - 1) { // target is leaf, tree traverse source
                tri_interactions.push_back(vector<int> {curr_target, 4 * curr_source, lev_target, lev_source + 1});
                tri_interactions.push_back(vector<int> {curr_target, 4 * curr_source + 1, lev_target, lev_source + 1});
                tri_interactions.push_back(vector<int> {curr_target, 4 * curr_source + 2, lev_target, lev_source + 1});
                tri_interactions.push_back(vector<int> {curr_target, 4 * curr_source + 3, lev_target, lev_source + 1});
            } else if (lev_source == run_information.fast_sum_tree_levels - 1) { // source is leaf, tree traverse target
                tri_interactions.push_back(vector<int> {4 * curr_target, curr_source, lev_target + 1, lev_source});
                tri_interactions.push_back(vector<int> {4 * curr_target + 1, curr_source, lev_target + 1, lev_source});
                tri_interactions.push_back(vector<int> {4 * curr_target + 2, curr_source, lev_target + 1, lev_source});
                tri_interactions.push_back(vector<int> {4 * curr_target + 3, curr_source, lev_target + 1, lev_source});
            } else { // neither is leaf
                if (particle_count_target >= particle_count_source) { // target has more points, refine target
                    tri_interactions.push_back(vector<int> {4 * curr_target, curr_source, lev_target + 1, lev_source});
                    tri_interactions.push_back(vector<int> {4 * curr_target + 1, curr_source, lev_target + 1, lev_source});
                    tri_interactions.push_back(vector<int> {4 * curr_target + 2, curr_source, lev_target + 1, lev_source});
                    tri_interactions.push_back(vector<int> {4 * curr_target + 3, curr_source, lev_target + 1, lev_source});
                } else { // source has more points, refine source
                    tri_interactions.push_back(vector<int> {curr_target, 4 * curr_source, lev_target, lev_source + 1});
                    tri_interactions.push_back(vector<int> {curr_target, 4 * curr_source + 1, lev_target, lev_source + 1});
                    tri_interactions.push_back(vector<int> {curr_target, 4 * curr_source + 2, lev_target, lev_source + 1});
                    tri_interactions.push_back(vector<int> {curr_target, 4 * curr_source + 3, lev_target, lev_source + 1});
                }
            }
        }
    }
}

void pp(run_config& run_information, vector<double>& modify, vector<double>& curr_state, vector<double>& area,
        interaction_pair& interact, vector<vector<vector<int>>>& fast_sum_tree_tri_points) {
    int target_i, source_j;
    vector<double> particle_i, particle_j, contribution;
    for (int i = 0; i < interact.count_target; i++) {
        target_i = fast_sum_tree_tri_points[interact.lev_target][interact.curr_target][i];
        vector<double> pos_change (3, 0);
        particle_i = slice(curr_state, run_information.info_per_point * target_i, 1, 3);
        for (int j = 0; j < interact.count_source; j++) {
            source_j = fast_sum_tree_tri_points[interact.lev_source][interact.curr_source][j];
            if (target_i != source_j) {
                particle_j = slice(curr_state, run_information.info_per_point * source_j, 1, 3);
                contribution = BVE_gfunc(particle_i, particle_j);
                scalar_mult(contribution, curr_state[run_information.info_per_point * source_j + 3] * area[source_j]);
                vec_add(pos_change, contribution);
            }
        }
        for (int j = 0; j < 3; j++) {
            modify[run_information.info_per_point * target_i + j] += pos_change[j];
        }
    }
}

void pc(run_config& run_information, vector<double>& modify, vector<double>& curr_state, vector<double>& area,
        interaction_pair& interact, vector<vector<vector<int>>>& fast_sum_tree_tri_points, vector<vector<vector<int>>>& fast_sum_icos_tri_verts,
        vector<vector<double>>& fast_sum_icos_verts) {
    vector<double> v1s, v2s, v3s, target_particle, placeholder1, placeholder2, placeholder3, bary_cord, source_particle;
    vector<double> func_vals (3 * run_information.interp_point_count, 0), func_val (3, 0);
    vector<double> alphas_x (run_information.interp_point_count, 0), alphas_y (run_information.interp_point_count, 0), alphas_z (run_information.interp_point_count, 0);
    int iv1s, iv2s, iv3s, point_index;
    double us, vs;
    char trans = 'N';
    int nrhs = 3, dim = run_information.interp_point_count, info;
    vector<double> interp_matrix (run_information.interp_point_count * run_information.interp_point_count, 0);
    vector<int> ipiv (run_information.interp_point_count, 0);
    vector<vector<double>> interp_points (run_information.interp_point_count, vector<double> (3, 0));
    fekete_init(interp_points, run_information.interp_degree);
    interp_mat_init(interp_matrix, interp_points, run_information.interp_degree, run_information.interp_point_count);
    dgetrf_(&dim, &dim, &*interp_matrix.begin(), &dim, &*ipiv.begin(), &info);
    iv1s = fast_sum_icos_tri_verts[interact.lev_source][interact.curr_source][0];
    iv2s = fast_sum_icos_tri_verts[interact.lev_source][interact.curr_source][1];
    iv3s = fast_sum_icos_tri_verts[interact.lev_source][interact.curr_source][2];
    v1s = fast_sum_icos_verts[iv1s];
    v2s = fast_sum_icos_verts[iv2s];
    v3s = fast_sum_icos_verts[iv3s];
    for (int i = 0; i < interact.count_target; i++) {
        point_index = fast_sum_tree_tri_points[interact.lev_target][interact.curr_target][i];
        target_particle = slice(curr_state, run_information.info_per_point * point_index, 1, 3);
        for (int j = 0; j < run_information.interp_point_count; j++) {
            us = interp_points[j][0];
            vs = interp_points[j][1];
            placeholder1 = v1s;
            placeholder2 = v2s;
            placeholder3 = v3s;
            scalar_mult(placeholder1, us);
            scalar_mult(placeholder2, vs);
            scalar_mult(placeholder3, 1.0 - us - vs);
            vec_add(placeholder1, placeholder2);
            vec_add(placeholder1, placeholder3);
            func_val = BVE_gfunc(target_particle, placeholder1);
            for (int k = 0; k < 3; k++) func_vals[j + run_information.interp_point_count * k] = func_val[k];
        }

        dgetrs_(&trans, &dim, &nrhs, &*interp_matrix.begin(), &dim, &*ipiv.begin(), &*func_vals.begin(), &dim, &info);
        if (info > 0) {
            cout << info << endl;
        }

        for (int j = 0; j < run_information.interp_point_count; j++) {
            alphas_x[j] = func_vals[j];
            alphas_y[j] = func_vals[j + run_information.interp_point_count];
            alphas_z[j] = func_vals[j + 2 * run_information.interp_point_count];
        }
        for (int j = 0; j < interact.count_source; j++) {
            point_index = fast_sum_tree_tri_points[interact.lev_source][interact.curr_source][j];
            source_particle = slice(curr_state, run_information.info_per_point * point_index, 1, 3);
            bary_cord = barycoords(v1s, v2s, v3s, source_particle);
            modify[run_information.info_per_point * i] += interp_eval(alphas_x, bary_cord[0], bary_cord[1], run_information.interp_degree) * curr_state[run_information.info_per_point * point_index + 3] * area[point_index];
            modify[run_information.info_per_point * i + 1] += interp_eval(alphas_y, bary_cord[0], bary_cord[1], run_information.interp_degree) * curr_state[run_information.info_per_point * point_index + 3] * area[point_index];
            modify[run_information.info_per_point * i + 2] += interp_eval(alphas_z, bary_cord[0], bary_cord[1], run_information.interp_degree) * curr_state[run_information.info_per_point * point_index + 3] * area[point_index];
        }
    }
}

void cp(run_config& run_information, vector<double>& modify, vector<double>& curr_state, vector<double>& area,
        interaction_pair& interact, vector<vector<vector<int>>>& fast_sum_tree_tri_points, vector<vector<vector<int>>>& fast_sum_icos_tri_verts,
        vector<vector<double>>& fast_sum_icos_verts) {
    int iv1, iv2, iv3, point_index;
    vector<double> v1, v2, v3, placeholder1, placeholder2, placeholder3, source_particle, target_particle, bary_cord;
    double u, v;
    vector<vector<double>> curr_points (run_information.interp_point_count, vector<double> (3, 0));
    vector<double> interptargets (3 * run_information.interp_point_count, 0), func_val (3, 0);
    char trans = 'N';
    int nrhs = 3, dim = run_information.interp_point_count, info;
    vector<double> alphas_x (run_information.interp_point_count, 0), alphas_y (run_information.interp_point_count, 0), alphas_z (run_information.interp_point_count, 0);
    vector<double> interp_matrix (run_information.interp_point_count * run_information.interp_point_count, 0);
    vector<int> ipiv (run_information.interp_point_count, 0);
    vector<vector<double>> interp_points (run_information.interp_point_count, vector<double> (3, 0));
    fekete_init(interp_points, run_information.interp_degree);
    interp_mat_init(interp_matrix, interp_points, run_information.interp_degree, run_information.interp_point_count);
    dgetrf_(&dim, &dim, &*interp_matrix.begin(), &dim, &*ipiv.begin(), &info);
    if (info > 0) {
        cout << info << endl;
    }
    iv1 = fast_sum_icos_tri_verts[interact.lev_target][interact.curr_target][0];
    iv2 = fast_sum_icos_tri_verts[interact.lev_target][interact.curr_target][1];
    iv3 = fast_sum_icos_tri_verts[interact.lev_target][interact.curr_target][2];
    v1 = fast_sum_icos_verts[iv1];
    v2 = fast_sum_icos_verts[iv2];
    v3 = fast_sum_icos_verts[iv3];
    for (int i = 0; i < run_information.interp_point_count; i++) {
        u = interp_points[i][0];
        v = interp_points[i][1];
        placeholder1 = v1;
        placeholder2 = v2;
        placeholder3 = v3;
        scalar_mult(placeholder1, u);
        scalar_mult(placeholder2, v);
        scalar_mult(placeholder3, 1.0 - u - v);
        vec_add(placeholder1, placeholder2);
        vec_add(placeholder1, placeholder3);
        curr_points[i] = placeholder1;
    }

    for (int i = 0; i < run_information.interp_point_count; i++) {
        for (int j = 0; j < interact.count_source; j++) {
            point_index = fast_sum_tree_tri_points[interact.lev_source][interact.curr_source][j];
            source_particle = slice(curr_state, run_information.info_per_point * point_index, 1, 3);
            func_val = BVE_gfunc(curr_points[i], source_particle);

            interptargets[i] += func_val[0] * curr_state[run_information.info_per_point * point_index + 3] * area[point_index];
            interptargets[i + run_information.interp_point_count] += func_val[1] * curr_state[run_information.info_per_point * point_index + 3] * area[point_index];
            interptargets[i + 2 * run_information.interp_point_count] += func_val[2] * curr_state[run_information.info_per_point * point_index + 3] * area[point_index];
        }
    }

    dgetrs_(&trans, &dim, &nrhs, &*interp_matrix.begin(), &dim, &*ipiv.begin(), &*interptargets.begin(), &dim, &info);
    if (info > 0) {
        cout << info << endl;
    }

    for (int i = 0; i < run_information.interp_point_count; i++) {
        alphas_x[i] = interptargets[i];
        alphas_y[i] = interptargets[i + run_information.interp_point_count];
        alphas_z[i] = interptargets[i + 2 * run_information.interp_point_count];
    }

    for (int i = 0; i < interact.count_target; i++) {
        point_index = fast_sum_tree_tri_points[interact.lev_target][interact.curr_target][i];
        target_particle = slice(curr_state, run_information.info_per_point * point_index, 1, 3);
        bary_cord = barycoords(v1, v2, v3, target_particle);
        modify[run_information.info_per_point * point_index] += interp_eval(alphas_x, bary_cord[0], bary_cord[1], run_information.interp_degree);
        modify[run_information.info_per_point * point_index + 1] += interp_eval(alphas_y, bary_cord[0], bary_cord[1], run_information.interp_degree);
        modify[run_information.info_per_point * point_index + 2] += interp_eval(alphas_z, bary_cord[0], bary_cord[1], run_information.interp_degree);
    }
}

void cc(run_config& run_information, vector<double>& modify, vector<double>& curr_state, vector<double>& area,
        interaction_pair& interact, vector<vector<vector<int>>>& fast_sum_tree_tri_points, vector<vector<vector<int>>>& fast_sum_icos_tri_verts,
        vector<vector<double>>& fast_sum_icos_verts) {
    int iv1, iv2, iv3, iv1s, iv2s, iv3s, point_index;
    vector<double> v1, v2, v3, placeholder1, placeholder2, placeholder3, v1s, v2s, v3s, func_vals (3 * run_information.interp_point_count, 0), func_val (3, 0), alphas_x (run_information.interp_point_count, 0), alphas_y (run_information.interp_point_count, 0), alphas_z (run_information.interp_point_count, 0);
    double u, v, us, vs;
    vector<vector<double>> curr_points (run_information.interp_point_count, vector<double> (3, 0));
    int nrhs = 3, dim = run_information.interp_point_count, info;
    char trans = 'N';
    vector<double> bary_cord, target_particle, source_particle;
    vector<double> interptargets (3 * run_information.interp_point_count, 0);
    vector<double> interp_matrix (run_information.interp_point_count * run_information.interp_point_count, 0);
    vector<int> ipiv (run_information.interp_point_count, 0);
    vector<vector<double>> interp_points (run_information.interp_point_count, vector<double> (3, 0));
    fekete_init(interp_points, run_information.interp_degree);
    interp_mat_init(interp_matrix, interp_points, run_information.interp_degree, run_information.interp_point_count);
    dgetrf_(&dim, &dim, &*interp_matrix.begin(), &dim, &*ipiv.begin(), &info);
    iv1 = fast_sum_icos_tri_verts[interact.lev_target][interact.curr_target][0];
    iv2 = fast_sum_icos_tri_verts[interact.lev_target][interact.curr_target][1];
    iv3 = fast_sum_icos_tri_verts[interact.lev_target][interact.curr_target][2];
    v1 = fast_sum_icos_verts[iv1];
    v2 = fast_sum_icos_verts[iv2];
    v3 = fast_sum_icos_verts[iv3];
    for (int i = 0; i < run_information.interp_point_count; i++) { // interpolation points in target triangle
        u = interp_points[i][0];
        v = interp_points[i][1];
        placeholder1 = v1;
        placeholder2 = v2;
        placeholder3 = v3;
        scalar_mult(placeholder1, u);
        scalar_mult(placeholder2, v);
        scalar_mult(placeholder3, 1.0 - u - v);
        vec_add(placeholder1, placeholder2);
        vec_add(placeholder1, placeholder3);
        curr_points[i] = placeholder1;
    }

    iv1s = fast_sum_icos_tri_verts[interact.lev_source][interact.curr_source][0];
    iv2s = fast_sum_icos_tri_verts[interact.lev_source][interact.curr_source][1];
    iv3s = fast_sum_icos_tri_verts[interact.lev_source][interact.curr_source][2];
    v1s = fast_sum_icos_verts[iv1s];
    v2s = fast_sum_icos_verts[iv2s];
    v3s = fast_sum_icos_verts[iv3s];
    for (int i = 0; i < run_information.interp_point_count; i++) { // loop across target interpolation points
        for (int j = 0; j < run_information.interp_point_count; j++) { // loop across source interpolation points
            // for each target interpolation point, interact with the source interpolation points
            us = interp_points[j][0];
            vs = interp_points[j][1];
            placeholder1 = v1s;
            placeholder2 = v2s;
            placeholder3 = v3s;
            scalar_mult(placeholder1, us);
            scalar_mult(placeholder2, vs);
            scalar_mult(placeholder3, 1.0 - us - vs);
            vec_add(placeholder1, placeholder2);
            vec_add(placeholder1, placeholder3);
            func_val = BVE_gfunc(curr_points[i], placeholder1);
            for (int k = 0; k < 3; k++) func_vals[j + run_information.interp_point_count * k] = func_val[k];
        }

        dgetrs_(&trans, &dim, &nrhs, &*interp_matrix.begin(), &dim, &*ipiv.begin(), &*func_vals.begin(), &dim, &info);
        if (info > 0) {
            cout << info << endl;
        }

        for (int j = 0; j < run_information.interp_point_count; j++) {
            alphas_x[j] = func_vals[j];
            alphas_y[j] = func_vals[j + run_information.interp_point_count];
            alphas_z[j] = func_vals[j + 2 * run_information.interp_point_count];
        }

        for (int j = 0; j < interact.count_source; j++) { // interpolate green's function into interior of source triangle
            point_index = fast_sum_tree_tri_points[interact.lev_source][interact.curr_source][j];
            source_particle = slice(curr_state, run_information.info_per_point * point_index, 1, 3);
            bary_cord = barycoords(v1s, v2s, v3s, source_particle);
            interptargets[i] += interp_eval(alphas_x, bary_cord[0], bary_cord[1], run_information.interp_degree) * curr_state[run_information.info_per_point * point_index + 3] * area[point_index];
            interptargets[i + run_information.interp_point_count] += interp_eval(alphas_y, bary_cord[0], bary_cord[1], run_information.interp_degree) * curr_state[run_information.info_per_point * point_index + 3] * area[point_index];
            interptargets[i + 2 * run_information.interp_point_count] += interp_eval(alphas_z, bary_cord[0], bary_cord[1], run_information.interp_degree) * curr_state[run_information.info_per_point * point_index + 3] * area[point_index];
        }
    }

    dgetrs_(&trans, &dim, &nrhs, &*interp_matrix.begin(), &dim, &*ipiv.begin(), &*interptargets.begin(), &dim, &info);

    if (info > 0) {
        cout << info << endl;
    }

    for (int i = 0; i < run_information.interp_point_count; i++) {
        alphas_x[i] = interptargets[i];
        alphas_y[i] = interptargets[i + run_information.interp_point_count];
        alphas_z[i] = interptargets[i + 2 * run_information.interp_point_count];
    }

    for (int i = 0; i < interact.count_target; i++) { // interpolate interaction into target triangle
        point_index = fast_sum_tree_tri_points[interact.lev_target][interact.curr_target][i];
        target_particle = slice(curr_state, run_information.info_per_point * point_index, 1, 3);
        bary_cord = barycoords(v1, v2, v3, target_particle);
        modify[run_information.info_per_point * point_index] += interp_eval(alphas_x, bary_cord[0], bary_cord[1], run_information.interp_degree);
        modify[run_information.info_per_point * point_index + 1] += interp_eval(alphas_y, bary_cord[0], bary_cord[1], run_information.interp_degree);
        modify[run_information.info_per_point * point_index + 2] += interp_eval(alphas_z, bary_cord[0], bary_cord[1], run_information.interp_degree);
    }
}

#endif
