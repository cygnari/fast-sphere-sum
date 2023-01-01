#include <iostream>
#include <cmath>
#include <array>
#include <vector>
#include <Accelerate/Accelerate.h>
#include <fstream>
#include <queue>

#include "helpers.h"

 // speed tests: O3, 1 evaluation of the sum, many_count = 10, theta = 0.7

 // time taken:     345184 microseconds for 2562 particles, 5 levels
 // time taken:     350029 microseconds for 2562 particles, 4 levels
 // time taken:     299788 microseconds for 2562 particles, 3 levels
 // time taken:     113980 microseconds for 2562 particles, 2 levels -- faster than direct sum

 // time taken:    2911411 microseconds for 10242 particles, 4 levels
 // time taken:     660812 microseconds for 10242 particles, 3 levels
 // time taken:    1065472 microseconds for 10242 particles, 2 levels

 // time taken:   40838338 microseconds for 40962 particles, 5 levels
 // time taken:    4730729 microseconds for 40962 particles, 4 levels
 // time taken:    4924663 microseconds for 40962 particles, 3 levels

 // time taken:  635798360 microseconds for 163842 particles, 6 levels
 // time taken:   49719799 microseconds for 163842 particles, 5 levels
 // time taken:   23078307 microseconds for 163842 particles, 4 levels
 // time taken:   66444843 microseconds for 163842 particles, 3 levels

 // time taken:  673310642 microseconds for 655362 particles, 6 levels
 // time taken:  129060434 microseconds for 655362 particles, 5 levels
 // time taken:  279946678 microseconds for 655362 particles, 4 levels
 // time taken: 1043940458 microseconds for 655362 particles, 3 levels


using namespace std;

void direct_sum(vector<double>& modify, vector<double>& curr_state, vector<vector<vector<int>>>& tri_points, int lev_target, int lev_source, int curr_target, int curr_source, double omega, double area) {
    int particle_count_target, particle_count_source; // do direct summation for a pair of triangles
    int target_i, source_j;
    particle_count_target = tri_points[lev_target][curr_target].size();
    particle_count_source = tri_points[lev_source][curr_source].size();
    vector<double> particle_i, particle_j, contribution;
    for (int i = 0; i < particle_count_target; i++) {
        target_i = tri_points[lev_target][curr_target][i];
        vector<double> pos_change {0, 0, 0};
        particle_i = slice(curr_state, 4 * target_i, 1, 3);
        for (int j = 0; j < particle_count_source; j++) {
            source_j = tri_points[lev_source][curr_source][j];
            if (target_i != source_j) {
                particle_j = slice(curr_state, 4 * source_j, 1, 3);
                contribution = BVE_gfunc(particle_i, particle_j);
                scalar_mult(contribution, curr_state[4 * source_j + 3] * area);
                vec_add(pos_change, contribution);
            }
        }
        for (int j = 0; j < 3; j++) {
            modify[4 * target_i + j] += pos_change[j];
        }
    }
}

void BVE_ffunc(vector<double>& modify, vector<double>& curr_state, vector<vector<double>>& vertices, vector<vector<vector<double>>>& tri_info, vector<vector<vector<int>>>& tri_verts, vector<vector<vector<int>>>& tri_points, double t, double delta_t, double omega, double area, int points, int levels, double radius, double theta, int many_count) {
    vector<vector<double>> target_points {{1, 0}, {0.5, 0.5}, {0, 1}, {0, 0.5}, {0, 0}, {0.5, 0}};
    int target_cluster_count = target_points.size();
    vector<vector<double>> source_points {{1, 0}, {0.5, 0.5}, {0, 1}, {0, 0.5}, {0, 0}, {0.5, 0}};
    int source_cluster_count = source_points.size();
    int curr_source, curr_target, lev_target, lev_source;
    int particle_count_target, particle_count_source;
    vector<double> center_target, center_source;
    double u, v, us, vs;
    double separation, distance;
    int iv1, iv2, iv3, iv1s, iv2s, iv3s;
    vector<double> v1, v2, v3, v1s, v2s, v3s;
    vector<double> bary_cord;
    int point_index;
    vector<vector<double>> curr_points (6, vector<double> (3, 0));
    vector<double> interptargets (18, 0);
    vector<double> placeholder1, placeholder2, placeholder3, source_particle, target_particle;
    vector<double> Aimatrix {1, 1, 1, 1, 1, 1, 1, 0.5, 0, 0, 0, 0.5, 0, 0.5, 1, 0.5, 0, 0, 0, 0.25, 0, 0, 0, 0, 1, 0.25, 0, 0, 0, 0.25, 0, 0.25, 1, 0.25, 0, 0};
    vector<double> func_vals (18, 0), func_val (3, 0);
    vector<double> alphas_x (6, 0), alphas_y (6, 0), alphas_z (6, 0);
    int ipiv[6];
    int info, dim = 6, nrhs = 3;
    char trans = 'N';
    dgetrf_(&dim, &dim, &*Aimatrix.begin(), &dim, ipiv, &info); // prefactor Ai matrix
    vector<vector<int>> tri_interactions;
    vector<int> curr_interact {0, 0, 0, 0};
    for (int i = 0; i < 20; i++) { // queue of triangle pairs to interact
        for (int j = 0; j < 20; j++) {
            tri_interactions.push_back({i, j, 0, 0});
        }
    }

    for (int i = 0; i < modify.size(); i++) modify[i] = 0; // zero out update

    while (tri_interactions.size() > 0) {

        curr_interact = tri_interactions.front(); // get triangle pair to interact
        curr_target = curr_interact[0];
        curr_source = curr_interact[1];
        lev_target = curr_interact[2];
        lev_source = curr_interact[3];
        tri_interactions.erase(tri_interactions.begin());
        particle_count_target = tri_points[lev_target][curr_target].size();
        particle_count_source = tri_points[lev_source][curr_source].size();
        if ((particle_count_target == 0) or (particle_count_source == 0)) continue; // if no work, continue to next
        center_target = slice(tri_info[lev_target][curr_target], 0, 1, 3);
        center_source = slice(tri_info[lev_source][curr_source], 0, 1, 3);
        distance = great_circ_dist(center_target, center_source, radius);
        separation = (tri_info[lev_target][curr_target][3] + tri_info[lev_source][curr_source][3]) / distance;
        if ((distance > 0) and (separation < theta)) { // triangles are well separated
            if (particle_count_target > many_count) { // target has many points, do cluster-X interaction
                iv1 = tri_verts[lev_target][curr_target][0];
                iv2 = tri_verts[lev_target][curr_target][1];
                iv3 = tri_verts[lev_target][curr_target][2];
                v1 = vertices[iv1];
                v2 = vertices[iv2];
                v3 = vertices[iv3];
                for (int i = 0; i < 6; i++) {
                    u = target_points[i][0];
                    v = target_points[i][1];
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
                if (particle_count_source > many_count) { // source has many particles, do C-C interaction
                    iv1s = tri_verts[lev_source][curr_source][0];
                    iv2s = tri_verts[lev_source][curr_source][1];
                    iv3s = tri_verts[lev_source][curr_source][2];
                    v1s = vertices[iv1s];
                    v2s = vertices[iv2s];
                    v3s = vertices[iv3s];
                    for (int i = 0; i < 6; i++) {
                        for (int j = 0; j < 6; j++) {
                            us = source_points[j][0];
                            vs = source_points[j][1];
                            placeholder1 = v1s;
                            placeholder2 = v2s;
                            placeholder3 = v3s;
                            scalar_mult(placeholder1, us);
                            scalar_mult(placeholder2, vs);
                            scalar_mult(placeholder3, 1.0 - us - vs);
                            vec_add(placeholder1, placeholder2);
                            vec_add(placeholder1, placeholder3);
                            func_val = BVE_gfunc(curr_points[i], placeholder1);
                            for (int k = 0; k < 3; k++) func_vals[j + 6 * k] = func_val[k];
                        }
                        dgetrs_(&trans, &dim, &nrhs, &*Aimatrix.begin(), &dim, ipiv, &*func_vals.begin(), &dim, &info);
                        for (int j = 0; j < 6; j++) {
                            alphas_x[j] = func_vals[j];
                            alphas_y[j] = func_vals[j+6];
                            alphas_z[j] = func_vals[j+12];
                        }

                        for (int j = 0; j < interptargets.size(); j++) interptargets[j] = 0;

                        for (int j = 0; j < 6; j++) {
                            point_index = tri_points[lev_source][curr_source][j];
                            source_particle = slice(curr_state, 4 * point_index, 1, 3);
                            bary_cord = barycoords(v1s, v2s, v3s, source_particle);
                            interptargets[j] += interp_eval(alphas_x, bary_cord[0], bary_cord[1]) * curr_state[4 * point_index + 3] * area;
                            interptargets[j+6] += interp_eval(alphas_y, bary_cord[0], bary_cord[1]) * curr_state[4 * point_index + 3] * area;
                            interptargets[j+12] += interp_eval(alphas_z, bary_cord[0], bary_cord[1]) * curr_state[4 * point_index + 3] * area;
                        }
                    }
                } else { // source has few particles, do C-P interaction
                    for (int i = 0; i < 6; i++) {
                        for (int j = 0; j < interptargets.size(); j++) interptargets[j] = 0;

                        for (int j = 0; j < particle_count_source; j++) {
                            point_index = tri_points[lev_source][curr_source][j];
                            source_particle = slice(curr_state, 4 * point_index, 1, 3);
                            func_val = BVE_gfunc(curr_points[i], source_particle);
                            interptargets[i] += func_val[0] * curr_state[4 * point_index + 3] * area;
                            interptargets[i+6] += func_val[1] * curr_state[4 * point_index + 3] * area;
                            interptargets[i+12] += func_val[2] * curr_state[4 * point_index + 3] * area;
                        }
                    }
                }
                dgetrs_(&trans, &dim, &nrhs, &*Aimatrix.begin(), &dim, ipiv, &*interptargets.begin(), &dim, &info);

                for (int i = 0; i < 6; i++) {
                    alphas_x[i] = interptargets[i];
                    alphas_y[i] = interptargets[i+6];
                    alphas_z[i] = interptargets[i+12];
                }

                for (int i = 0; i < particle_count_target; i++) {
                    point_index = tri_points[lev_target][curr_target][i];
                    target_particle = slice(curr_state, 4 * point_index, 1, 3);
                    bary_cord = barycoords(v1, v2, v3, target_particle);
                    modify[4 * i] += interp_eval(alphas_x, bary_cord[0], bary_cord[1]);
                    modify[4 * i + 1] += interp_eval(alphas_y, bary_cord[0], bary_cord[1]);
                    modify[4 * i + 2] += interp_eval(alphas_z, bary_cord[0], bary_cord[1]);
                }
            } else { // target has few points, do particle-X interaction
                if (particle_count_source > many_count) { // source has lots of points, P-C interaction
                    iv1s = tri_verts[lev_source][curr_source][0];
                    iv2s = tri_verts[lev_source][curr_source][1];
                    iv3s = tri_verts[lev_source][curr_source][2];
                    v1s = vertices[iv1s];
                    v2s = vertices[iv2s];
                    v3s = vertices[iv3s];
                    for (int i = 0; i < particle_count_target; i++) {
                        point_index = tri_points[lev_target][curr_target][i];
                        target_particle = slice(curr_state, 4 * point_index, 1, 3);
                        for (int j = 0; j < 6; j++) {
                            us = source_points[j][0];
                            vs = source_points[j][1];
                            placeholder1 = v1s;
                            placeholder2 = v2s;
                            placeholder3 = v3s;
                            scalar_mult(placeholder1, us);
                            scalar_mult(placeholder2, vs);
                            scalar_mult(placeholder3, 1.0 - us - vs);
                            vec_add(placeholder1, placeholder2);
                            vec_add(placeholder1, placeholder3);
                            func_val = BVE_gfunc(target_particle, placeholder1);
                            for (int k = 0; k < 3; k++) func_vals[j + 6 * k] = func_val[k];
                        }
                        dgetrs_(&trans, &dim, &nrhs, &*Aimatrix.begin(), &dim, ipiv, &*func_vals.begin(), &dim, &info);

                        for (int j = 0; j < 6; j++) {
                            alphas_x[j] = func_vals[j];
                            alphas_y[j] = func_vals[j+6];
                            alphas_z[j] = func_vals[j+12];
                        }
                        for (int j = 0; j < particle_count_source; j++) {
                            point_index = tri_points[lev_source][curr_source][j];
                            source_particle = slice(curr_state, 4 * point_index, 1, 3);
                            bary_cord = barycoords(v1s, v2s, v3s, source_particle);
                            modify[4 * i] += interp_eval(alphas_x, bary_cord[0], bary_cord[1]) * curr_state[4 * point_index + 3] * area;
                            modify[4 * i + 1] += interp_eval(alphas_y, bary_cord[0], bary_cord[1]) * curr_state[4 * point_index + 3] * area;
                            modify[4 * i + 2] += interp_eval(alphas_z, bary_cord[0], bary_cord[1]) * curr_state[4 * point_index + 3] * area;
                        }
                    }
                } else { // source has few points, P-P interaction
                    direct_sum(modify, curr_state, tri_points, lev_target, lev_source, curr_target, curr_source, omega, area);
                }
            }
        } else { // not well separated
            if ((particle_count_source < many_count) and (particle_count_target < many_count)) { // both have few points, P-P interaction
                direct_sum(modify, curr_state, tri_points, lev_target, lev_source, curr_target, curr_source, omega, area);
            } else if ((lev_target == levels - 1) and (lev_source == levels - 1)) { // both are leaves, do P-P interaction
                direct_sum(modify, curr_state, tri_points, lev_target, lev_source, curr_target, curr_source, omega, area);
            } else if (lev_target == levels - 1) { // target is leaf, source is not, tree traverse
                tri_interactions.push_back(vector<int> {curr_target, 4 * curr_source, lev_target, lev_source + 1});
                tri_interactions.push_back(vector<int> {curr_target, 4 * curr_source + 1, lev_target, lev_source + 1});
                tri_interactions.push_back(vector<int> {curr_target, 4 * curr_source + 2, lev_target, lev_source + 1});
                tri_interactions.push_back(vector<int> {curr_target, 4 * curr_source + 3, lev_target, lev_source + 1});
            } else if (lev_source == levels - 1) { // source is leaf, target is not, tree traverse
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
    scalar_mult(modify, -1.0 / (4.0 * M_PI));
    for (int i = 0; i < points; i++) modify[4 * i + 3] = -2 * omega * modify[4 * i + 2];
}

#define point_count 2562

int main() {

    double delta_t = 0.01, end_t = 1;
    double omega = 2 * M_PI;
    int time_steps = end_t / delta_t;
    double area = (4 * M_PI) / point_count;
    int icos_levels = 2;
    double radius = 1.0;
    double phi = (1 + sqrt(5)) / 2;
    double theta = 0.7;
    int many_count = 10;

    vector<double> curr_state(4 * point_count); // 0 is x_pos, 1 is y_pos, 2 is z_pos, 3 is vorticity
    vector<double> c_1(4 * point_count, 0);

    vector<vector<double>> vertices; // list of mesh vertices
    vector<vector<vector<double>>> triangle_info; // contains triangle centers, radii
    vector<vector<vector<int>>> triangle_verts; // triangle vertex indices
    vector<vector<vector<int>>> tri_points (icos_levels); // points in each triangle
    vector<vector<int>> point_locs (icos_levels, vector<int> (point_count, 0)); // triangle each point is in

    fstream file("../points.csv");
    string line, word;

    ofstream write_out;
    write_out.open("fast_output_test.csv", ofstream::out | ofstream::trunc);

    for (int i = 0; i < point_count; i++) {
        getline(file, line);
        stringstream str(line);
        for (int j = 0; j < 4; j++) {
            getline(str, word, ',');
            curr_state[4 * i + j] = stod(word);
        }
    }

    icos_init(vertices, triangle_info, triangle_verts, radius, icos_levels-1);

    chrono::steady_clock::time_point begin = chrono::steady_clock::now();
    points_assign(triangle_verts, vertices, curr_state, tri_points, point_locs, icos_levels, point_count);
    BVE_ffunc(c_1, curr_state, vertices, triangle_info, triangle_verts, tri_points, 0, delta_t, omega, area, point_count, icos_levels, radius, 0.7, many_count);
    chrono::steady_clock::time_point end = chrono::steady_clock::now();

    for (int i = 0; i < point_count; i++) {
        vector<double> projected = slice(c_1, 4 * i, 1, 3);
        project_to_sphere(projected, 1);
        for (int j = 0; j < 3; j++) c_1[4 * i + j] = projected[j]; // reproject points to surface of sphere
        write_out << c_1[4 * i] << "," << c_1[4 * i + 1] << "," << c_1[4 * i + 2] << "," << c_1[4 * i + 3] << "\n"; // write position
    }

    cout << "time taken: " << chrono::duration_cast<chrono::microseconds>(end - begin).count() << " microseconds" << endl;

    write_out.close();
    return 0;
}
