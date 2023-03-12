#ifndef fast_H
#define fast_H

#include <iostream>
#include <cmath>
#include <vector>
#include <queue>
#include <chrono>
#include <Accelerate/Accelerate.h>

#include "helpers.h"
#include "struct_list.h"
#include "icos_utils.h"
#include "interp_utils.h"

void tree_traverse(vector<interaction_pair>& interactions, icos_struct& icos1, double theta, int many_count) {
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

    while (tri_interactions.size() > 0) {
        curr_interact = tri_interactions.front(); // get triangle pair to interact
        curr_target = curr_interact[0];
        curr_source = curr_interact[1];
        lev_target = curr_interact[2];
        lev_source = curr_interact[3];
        // cout << "here" << endl;
        tri_interactions.erase(tri_interactions.begin());
        particle_count_target = icos1.tri_points[lev_target][curr_target].size();
        particle_count_source = icos1.tri_points[lev_source][curr_source].size();
        if ((particle_count_target == 0) or (particle_count_source == 0)) continue; // if no work, continue to next
        center_target = slice(icos1.tri_info[lev_target][curr_target], 0, 1, 3);
        center_source = slice(icos1.tri_info[lev_source][curr_source], 0, 1, 3);
        distance = great_circ_dist(center_target, center_source, icos1.radius);
        separation = (icos1.tri_info[lev_target][curr_target][3] + icos1.tri_info[lev_source][curr_source][3]) / distance;

        if ((distance > 0) and (separation < theta)) { // triangles are well separated
            interaction_pair new_interact = {lev_target, lev_source, curr_target, curr_source, particle_count_target, particle_count_source, ""};
            if (particle_count_target > many_count) {
                new_interact.type += "c";
            } else {
                new_interact.type += "p";
            }
            if (particle_count_source > many_count) {
                new_interact.type += "c";
            } else {
                new_interact.type += "p";
            }
            interactions.push_back(new_interact);
        } else {
            if ((particle_count_target < many_count) and (particle_count_source < many_count)) { // both have few particles, pp
                interaction_pair new_interact = {lev_target, lev_source, curr_target, curr_source, particle_count_target, particle_count_source, "pp"};
                interactions.push_back(new_interact);
            } else if ((lev_target == icos1.levels - 1) and (lev_source == icos1.levels - 1)) { // both are leaves, pp
                interaction_pair new_interact = {lev_target, lev_source, curr_target, curr_source, particle_count_target, particle_count_source, "pp"};
                interactions.push_back(new_interact);
            } else if (lev_target == icos1.levels - 1) { // target is leaf, tree traverse source
                tri_interactions.push_back(vector<int> {curr_target, 4 * curr_source, lev_target, lev_source + 1});
                tri_interactions.push_back(vector<int> {curr_target, 4 * curr_source + 1, lev_target, lev_source + 1});
                tri_interactions.push_back(vector<int> {curr_target, 4 * curr_source + 2, lev_target, lev_source + 1});
                tri_interactions.push_back(vector<int> {curr_target, 4 * curr_source + 3, lev_target, lev_source + 1});
            } else if (lev_source == icos1.levels - 1) { // source is leaf, tree traverse target
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

void pp(vector<double>& modify, vector<double>& curr_state, vector<double>& area, icos_struct icos1, interaction_pair interact) {
// void direct_sum(vector<double>& modify, vector<double>& curr_state, vector<vector<vector<int>>>& tri_points, int lev_target, int lev_source, int curr_target, int curr_source, double omega, vector<double>& area) { // direct sum for P-P interaction
    // int particle_count_target, particle_count_source; // do direct summation for a pair of triangles
    int target_i, source_j;
    // particle_count_target = tri_points[lev_target][curr_target].size();
    // particle_count_source = tri_points[lev_source][curr_source].size();
    vector<double> particle_i, particle_j, contribution;
    for (int i = 0; i < interact.count_target; i++) {
        target_i = icos1.tri_points[interact.lev_target][interact.curr_target][i];
        // vector<double> pos_change {0, 0, 0};
        vector<double> pos_change (3, 0);
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

void pc(vector<double>& modify, vector<double>& curr_state, vector<double>& area, icos_struct icos1, interp_struct interp1, interaction_pair interact) {
    vector<double> v1s, v2s, v3s, target_particle, placeholder1, placeholder2, placeholder3, bary_cord, source_particle;
    vector<double> func_vals (3 * interp1.point_count, 0), func_val (3, 0);
    vector<double> alphas_x (interp1.point_count, 0), alphas_y (interp1.point_count, 0), alphas_z (interp1.point_count, 0);
    int iv1s, iv2s, iv3s, point_index;
    double us, vs;
    char trans = 'N';
    int nrhs = 3, dim = interp1.point_count, info;
    iv1s = icos1.tri_verts[interact.lev_source][interact.curr_source][0];
    iv2s = icos1.tri_verts[interact.lev_source][interact.curr_source][1];
    iv3s = icos1.tri_verts[interact.lev_source][interact.curr_source][2];
    v1s = icos1.verts[iv1s];
    v2s = icos1.verts[iv2s];
    v3s = icos1.verts[iv3s];
    for (int i = 0; i < interact.count_target; i++) {
        point_index = icos1.tri_points[interact.lev_target][interact.curr_target][i];
        target_particle = slice(curr_state, 5 * point_index, 1, 3);
        for (int j = 0; j < interp1.point_count; j++) {
            us = interp1.interp_points[j][0];
            vs = interp1.interp_points[j][1];
            placeholder1 = v1s;
            placeholder2 = v2s;
            placeholder3 = v3s;
            scalar_mult(placeholder1, us);
            scalar_mult(placeholder2, vs);
            scalar_mult(placeholder3, 1.0 - us - vs);
            vec_add(placeholder1, placeholder2);
            vec_add(placeholder1, placeholder3);
            func_val = BVE_gfunc(target_particle, placeholder1);
            for (int k = 0; k < 3; k++) func_vals[j + interp1.point_count * k] = func_val[k];
        }

        dgetrs_(&trans, &dim, &nrhs, &*interp1.matrix.begin(), &dim, &*interp1.ipiv.begin(), &*func_vals.begin(), &dim, &info);
        if (info > 0) {
            cout << info << endl;
        }

        for (int j = 0; j < interp1.point_count; j++) {
            alphas_x[j] = func_vals[j];
            alphas_y[j] = func_vals[j + interp1.point_count];
            alphas_z[j] = func_vals[j + 2 * interp1.point_count];
        }
        for (int j = 0; j < interact.count_source; j++) {
            point_index = icos1.tri_points[interact.lev_source][interact.curr_source][j];
            source_particle = slice(curr_state, 5 * point_index, 1, 3);
            bary_cord = barycoords(v1s, v2s, v3s, source_particle);
            modify[5 * i] += interp_eval(alphas_x, bary_cord[0], bary_cord[1], interp1.degree) * curr_state[5 * point_index + 3] * area[point_index];
            modify[5 * i + 1] += interp_eval(alphas_y, bary_cord[0], bary_cord[1], interp1.degree) * curr_state[5 * point_index + 3] * area[point_index];
            modify[5 * i + 2] += interp_eval(alphas_z, bary_cord[0], bary_cord[1], interp1.degree) * curr_state[5 * point_index + 3] * area[point_index];
        }
    }
}

void cp(vector<double>& modify, vector<double>& curr_state, vector<double>& area, icos_struct icos1, interp_struct interp1, interaction_pair interact) {
    int iv1, iv2, iv3, point_index;
    vector<double> v1, v2, v3, placeholder1, placeholder2, placeholder3, source_particle, target_particle, bary_cord;
    double u, v;
    vector<vector<double>> curr_points (interp1.point_count, vector<double> (3, 0));
    vector<double> interptargets (3 * interp1.point_count, 0), func_val (3, 0);
    char trans = 'N';
    int nrhs = 3, dim = interp1.point_count, info;
    vector<double> alphas_x (interp1.point_count, 0), alphas_y (interp1.point_count, 0), alphas_z (interp1.point_count, 0);
    iv1 = icos1.tri_verts[interact.lev_target][interact.curr_target][0];
    iv2 = icos1.tri_verts[interact.lev_target][interact.curr_target][1];
    iv3 = icos1.tri_verts[interact.lev_target][interact.curr_target][2];
    v1 = icos1.verts[iv1];
    v2 = icos1.verts[iv2];
    v3 = icos1.verts[iv3];
    for (int i = 0; i < interp1.point_count; i++) {
        u = interp1.interp_points[i][0];
        v = interp1.interp_points[i][1];
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

    for (int i = 0; i < interp1.point_count; i++) {
        // for (int j = 0; j < interptargets.size(); j++) interptargets[j] = 0;

        for (int j = 0; j < interact.count_source; j++) {
            point_index = icos1.tri_points[interact.lev_source][interact.curr_source][j];
            source_particle = slice(curr_state, 5 * point_index, 1, 3);
            func_val = BVE_gfunc(curr_points[i], source_particle);
            // cout << func_val[0] << endl;
            // cout << curr_points[0][0] << endl;
            interptargets[i] += func_val[0] * curr_state[5 * point_index + 3] * area[point_index];
            interptargets[i + interp1.point_count] += func_val[1] * curr_state[5 * point_index + 3] * area[point_index];
            interptargets[i + 2 * interp1.point_count] += func_val[2] * curr_state[5 * point_index + 3] * area[point_index];
        }
    }

    dgetrs_(&trans, &dim, &nrhs, &*interp1.matrix.begin(), &dim, &*interp1.ipiv.begin(), &*interptargets.begin(), &dim, &info);
    if (info > 0) {
        cout << info << endl;
    }

    for (int i = 0; i < interp1.point_count; i++) {
        alphas_x[i] = interptargets[i];
        alphas_y[i] = interptargets[i + interp1.point_count];
        alphas_z[i] = interptargets[i + 2 * interp1.point_count];
    }

    for (int i = 0; i < interact.count_target; i++) {
        point_index = icos1.tri_points[interact.lev_target][interact.curr_target][i];
        target_particle = slice(curr_state, 5 * point_index, 1, 3);
        bary_cord = barycoords(v1, v2, v3, target_particle);
        modify[5 * i] += interp_eval(alphas_x, bary_cord[0], bary_cord[1], interp1.degree);
        modify[5 * i + 1] += interp_eval(alphas_y, bary_cord[0], bary_cord[1], interp1.degree);
        modify[5 * i + 2] += interp_eval(alphas_z, bary_cord[0], bary_cord[1], interp1.degree);
    }
}

void cc(vector<double>& modify, vector<double>& curr_state, vector<double>& area, icos_struct icos1, interp_struct interp1, interaction_pair interact) {
    int iv1, iv2, iv3, iv1s, iv2s, iv3s, point_index;
    vector<double> v1, v2, v3, placeholder1, placeholder2, placeholder3, v1s, v2s, v3s, func_vals (3 * interp1.point_count, 0), func_val (3, 0), alphas_x (interp1.point_count, 0), alphas_y (interp1.point_count, 0), alphas_z (interp1.point_count, 0);
    double u, v, us, vs;
    vector<vector<double>> curr_points (interp1.point_count, vector<double> (3, 0));
    int nrhs = 3, dim = interp1.point_count, info;
    char trans = 'N';
    vector<double> bary_cord, target_particle, source_particle;
    vector<double> interptargets (3 * interp1.point_count, 0);
    iv1 = icos1.tri_verts[interact.lev_target][interact.curr_target][0];
    iv2 = icos1.tri_verts[interact.lev_target][interact.curr_target][1];
    iv3 = icos1.tri_verts[interact.lev_target][interact.curr_target][2];
    v1 = icos1.verts[iv1];
    v2 = icos1.verts[iv2];
    v3 = icos1.verts[iv3];
    for (int i = 0; i < interp1.point_count; i++) { // interpolation points in target triangle
        u = interp1.interp_points[i][0];
        v = interp1.interp_points[i][1];
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

    iv1s = icos1.tri_verts[interact.lev_source][interact.curr_source][0];
    iv2s = icos1.tri_verts[interact.lev_source][interact.curr_source][1];
    iv3s = icos1.tri_verts[interact.lev_source][interact.curr_source][2];
    v1s = icos1.verts[iv1s];
    v2s = icos1.verts[iv2s];
    v3s = icos1.verts[iv3s];
    for (int i = 0; i < interp1.point_count; i++) { // loop across target interpolation points
        for (int j = 0; j < interp1.point_count; j++) { // loop across source interpolation points
            // for each target interpolation point, interact with the source interpolation points
            us = interp1.interp_points[j][0];
            vs = interp1.interp_points[j][1];
            placeholder1 = v1s;
            placeholder2 = v2s;
            placeholder3 = v3s;
            scalar_mult(placeholder1, us);
            scalar_mult(placeholder2, vs);
            scalar_mult(placeholder3, 1.0 - us - vs);
            vec_add(placeholder1, placeholder2);
            vec_add(placeholder1, placeholder3);
            func_val = BVE_gfunc(curr_points[i], placeholder1);
            for (int k = 0; k < 3; k++) func_vals[j + interp1.point_count * k] = func_val[k];
        }
        // cout << "here 1" << endl;
        dgetrs_(&trans, &dim, &nrhs, &*interp1.matrix.begin(), &dim, &*interp1.ipiv.begin(), &*func_vals.begin(), &dim, &info);
        // cout << "here 2" << endl;
        if (info > 0) {
            cout << info << endl;
        }

        for (int j = 0; j < interp1.point_count; j++) {
            alphas_x[j] = func_vals[j];
            alphas_y[j] = func_vals[j + interp1.point_count];
            alphas_z[j] = func_vals[j + 2 * interp1.point_count];
        }

        // for (int j = 0; j < interptargets.size(); j++) interptargets[j] = 0;
        // fill(interptargets.begin(), interptargets.end(), 0); // zero out vector

        for (int j = 0; j < interact.count_source; j++) { // interpolate green's function into interior of source triangle
            point_index = icos1.tri_points[interact.lev_source][interact.curr_source][j];
            source_particle = slice(curr_state, 5 * point_index, 1, 3);
            bary_cord = barycoords(v1s, v2s, v3s, source_particle);
            interptargets[i] += interp_eval(alphas_x, bary_cord[0], bary_cord[1], interp1.degree) * curr_state[5 * point_index + 3] * area[point_index];
            interptargets[i + interp1.point_count] += interp_eval(alphas_y, bary_cord[0], bary_cord[1], interp1.degree) * curr_state[5 * point_index + 3] * area[point_index];
            interptargets[i + 2 * interp1.point_count] += interp_eval(alphas_z, bary_cord[0], bary_cord[1], interp1.degree) * curr_state[5 * point_index + 3] * area[point_index];
        }
        // cout << "here 3" << endl;
    }



    dgetrs_(&trans, &dim, &nrhs, &*interp1.matrix.begin(), &dim, &*interp1.ipiv.begin(), &*interptargets.begin(), &dim, &info);
    // cout << "here 4" << endl;
    if (info > 0) {
        cout << info << endl;
    }

    // cout << "here 5" << endl;
    for (int i = 0; i < interp1.point_count; i++) {
        alphas_x[i] = interptargets[i];
        alphas_y[i] = interptargets[i + interp1.point_count];
        alphas_z[i] = interptargets[i + 2 * interp1.point_count];
    }
    // cout << "here 5" << endl;

    for (int i = 0; i < interact.count_target; i++) { // interpolate interaction into target triangle
        point_index = icos1.tri_points[interact.lev_target][interact.curr_target][i];
        target_particle = slice(curr_state, 5 * point_index, 1, 3);
        bary_cord = barycoords(v1, v2, v3, target_particle);
        modify[5 * point_index] += interp_eval(alphas_x, bary_cord[0], bary_cord[1], interp1.degree);
        modify[5 * point_index + 1] += interp_eval(alphas_y, bary_cord[0], bary_cord[1], interp1.degree);
        modify[5 * point_index + 2] += interp_eval(alphas_z, bary_cord[0], bary_cord[1], interp1.degree);
    }
    // cout << "here 6" << endl;
}

#endif
