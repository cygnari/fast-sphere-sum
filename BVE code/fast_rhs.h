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

void tree_traverse(vector<interaction_pair>& interactions, icos_struct& icos1, double theta, int many_count) {
    int curr_source, curr_target, lev_target, lev_source;
    int particle_count_target, particle_count_source;
    vector<double> center_target, center_source;
    double separation, distance;
    vector<vector<int>> tri_interactions;
    vector<int> curr_interact {0, 0, 0, 0};

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

void cp(vector<double>& modify, vector<double>& curr_state, icos_struct icos1, interaction_pair interact, double omega) {
    // stuff
}

void cc(vector<double>& modify, vector<double>& curr_state, icos_struct icos1, interaction_pair interact, double omega) {
    // stuff
}

#endif
