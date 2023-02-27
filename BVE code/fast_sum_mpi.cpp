#include <iostream>
#include <cmath>
#include <array>
#include <vector>
// #include <Accelerate/Accelerate.h>
#include <fstream>
#include <queue>
#include <sstream>
#include <chrono>
#include <cassert>
#include <mpi.h>

#include "helpers.h"

 // time taken: 1365873 microseconds for 1 time step, 2562 points, O3
 // time taken: 102833420 microseconds for 101 time steps, 2562 points, O3

using namespace std;

void direct_sum(vector<double>& modify, vector<double>& curr_state, vector<vector<vector<int>>>& tri_points, int lev_target, int lev_source, int curr_target, int curr_source, double omega, vector<double>& area) {
    int particle_count_target, particle_count_source; // do direct summation for a pair of triangles
    int target_i, source_j;
    particle_count_target = tri_points[lev_target][curr_target].size();
    particle_count_source = tri_points[lev_source][curr_source].size();
    vector<double> particle_i, particle_j, contribution;
    for (int i = 0; i < particle_count_target; i++) {
        target_i = tri_points[lev_target][curr_target][i];
        vector<double> pos_change {0, 0, 0};
        particle_i = slice(curr_state, 5 * target_i, 1, 3);
        for (int j = 0; j < particle_count_source; j++) {
            source_j = tri_points[lev_source][curr_source][j];
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

void BVE_ffunc(vector<double>& modify, vector<double>& curr_state, vector<vector<double>>& vertices, vector<vector<vector<double>>>& tri_info, vector<vector<vector<int>>>& tri_verts, vector<vector<vector<int>>>& tri_points, double t, double delta_t, double omega, vector<double>& area, int points, int levels, double radius, double theta, int many_count, int tri_lb, int tri_ub) {
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
    for (int i = tri_lb; i < tri_ub; i++) { // queue of triangle pairs to interact
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
                            source_particle = slice(curr_state, 5 * point_index, 1, 3);
                            bary_cord = barycoords(v1s, v2s, v3s, source_particle);
                            interptargets[j] += interp_eval(alphas_x, bary_cord[0], bary_cord[1], 2) * curr_state[5 * point_index + 3] * area[point_index];
                            interptargets[j+6] += interp_eval(alphas_y, bary_cord[0], bary_cord[1], 2) * curr_state[5 * point_index + 3] * area[point_index];
                            interptargets[j+12] += interp_eval(alphas_z, bary_cord[0], bary_cord[1], 2) * curr_state[5 * point_index + 3] * area[point_index];
                        }
                    }
                } else { // source has few particles, do C-P interaction
                    for (int i = 0; i < 6; i++) {
                        for (int j = 0; j < interptargets.size(); j++) interptargets[j] = 0;

                        for (int j = 0; j < particle_count_source; j++) {
                            point_index = tri_points[lev_source][curr_source][j];
                            source_particle = slice(curr_state, 5 * point_index, 1, 3);
                            func_val = BVE_gfunc(curr_points[i], source_particle);
                            interptargets[i] += func_val[0] * curr_state[5 * point_index + 3] * area[point_index];
                            interptargets[i+6] += func_val[1] * curr_state[5 * point_index + 3] * area[point_index];
                            interptargets[i+12] += func_val[2] * curr_state[5 * point_index + 3] * area[point_index];
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
                    target_particle = slice(curr_state, 5 * point_index, 1, 3);
                    bary_cord = barycoords(v1, v2, v3, target_particle);
                    modify[5 * i] += interp_eval(alphas_x, bary_cord[0], bary_cord[1], 2);
                    modify[5 * i + 1] += interp_eval(alphas_y, bary_cord[0], bary_cord[1], 2);
                    modify[5 * i + 2] += interp_eval(alphas_z, bary_cord[0], bary_cord[1], 2);
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
                        target_particle = slice(curr_state, 5 * point_index, 1, 3);
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
                            source_particle = slice(curr_state, 5 * point_index, 1, 3);
                            bary_cord = barycoords(v1s, v2s, v3s, source_particle);
                            modify[5 * i] += interp_eval(alphas_x, bary_cord[0], bary_cord[1], 2) * curr_state[5 * point_index + 3] * area[point_index];
                            modify[5 * i + 1] += interp_eval(alphas_y, bary_cord[0], bary_cord[1], 2) * curr_state[5 * point_index + 3] * area[point_index];
                            modify[5 * i + 2] += interp_eval(alphas_z, bary_cord[0], bary_cord[1], 2) * curr_state[5 * point_index + 3] * area[point_index];
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
    for (int i = 0; i < points; i++) modify[5 * i + 3] = -2 * omega * modify[5 * i + 2];
}

int particle_processor(int points_per_rank, int point_id, int total_ranks) {
    int processor = point_id / points_per_rank;
    // cout << processor << " " << total_ranks << endl;
    if (processor < total_ranks) {
        return processor;
    } else {
        return total_ranks - 1;
    }
}

int processor_particle_count(int id, int total_P, int points) {
    int points_per_rank = points / total_P;
    if (id < total_P - 1) {
        return points_per_rank;
    } else {
        return points - (total_P - 1) * points_per_rank;
    }
}

void sync_buffer(vector<double>& buffer, int point_count, int ID, vector<int>& particle_thread, MPI_Win *win) {
    MPI_Win_fence(0, *win);
    int thread;
    for (int i = 0; i < point_count; i++) {
        thread = particle_thread[i];
        if (thread != ID) {
            MPI_Get(&buffer[5 * i], 5, MPI_DOUBLE, thread, 5 * i, 5, MPI_DOUBLE, *win);
        }
    }
    MPI_Win_fence(0, *win);
}

// #define point_count 655362

int main(int argc, char** argv) {

    MPI_Init(&argc, &argv);
    int P;
    int ID;
    MPI_Status status;
    MPI_Comm_size(MPI_COMM_WORLD, &P);
    MPI_Comm_rank(MPI_COMM_WORLD, &ID);
    MPI_Win win_curr_state;
    MPI_Win win_inter1;
    MPI_Win win_inter2;
    MPI_Win win_inter3;
    MPI_Win win_c1234;

    double delta_t = 0.01, end_t = 1;
    double omega = 2 * M_PI;
    int time_steps = end_t / delta_t;
    // double area = (4 * M_PI) / point_count;
    int icos_levels = 3;
    double radius = 1.0;
    double phi = (1 + sqrt(5)) / 2;
    double theta = 0.7;
    int many_count = 10;
    int point_count = 10242, tri_count = 20480, max_points = 1000000;

    int points_per_rank = point_count / P;
    int tris_per_thread = 20 / P;
    int own_points, own_tris, tri_lb, tri_ub;;

    chrono::steady_clock::time_point begin_time, end_time;

    if (ID < P - 1) {
        own_points = points_per_rank;
        own_tris = tris_per_thread;
        tri_lb = ID * tris_per_thread;
        tri_ub = (ID + 1) * tris_per_thread;
    } else {
        own_points = point_count - (P - 1) * points_per_rank;
        own_tris = 20 - (P - 1) * tris_per_thread;
        tri_lb = 20 - own_tris;
        tri_ub = 20;
    }

    vector<double> curr_state(5 * point_count); // 0 is x_pos, 1 is y_pos, 2 is z_pos, 3 is vorticity
    vector<double> swap_state(5 * point_count); // move particles around so that triangles have contiguous sections
    vector<double> c_1(5 * point_count, 0);
    vector<double> c_2(5 * point_count, 0);
    vector<double> c_3(5 * point_count, 0);
    vector<double> c_4(5 * point_count, 0);
    vector<double> c12(5 * point_count, 0);
    vector<double> c123(5 * point_count, 0);
    vector<double> c1234(5 * point_count, 0);
    vector<double> intermediate_1(5 * point_count);
    vector<double> intermediate_2(5 * point_count);
    vector<double> intermediate_3(5 * point_count);
    vector<double> area(point_count, 0);
    vector<int> state;

    vector<vector<double>> vertices; // list of mesh vertices
    vector<vector<vector<double>>> triangle_info; // contains triangle centers, radii
    vector<vector<vector<int>>> triangle_verts; // triangle vertex indices
    vector<vector<vector<int>>> tri_points (icos_levels); // points in each triangle
    vector<vector<int>> point_locs (icos_levels, vector<int> (point_count, 0)); // triangle each point is in

    vector<vector<int>> triangles(tri_count, vector<int> (3)); // triangles[i] is a length 3 int vector containing the indices of the points of the vertices
    vector<vector<int>> vert_tris(point_count); // vert_tris[i] is the int vector containing the indices of the triangles adjacent to point i
    vector<vector<int>> parent_verts(max_points, vector<int> (2, 0)); //

    vector<int> particle_thread(point_count, 0); // particle_thread[i] is the processor where particle i is located
    // vector<bool> particle_in_curr_thread(point_count, false); // true if particle i is in the current thread
    vector<int> triangle_thread(20); // triangle_thread[i] is the processor where top level triangle i is located
    // vector<bool> tri_in_curr_thread(20, false); // true if triangle i is in the current thread
    vector<int> tri_lbs(P, 0); // lower bound on triangle indices in thread i
    vector<int> tri_ubs(P, 0); // upper bound
    // vector<int> new_order(point_count, 0); // new_order[i] is the old index of the point
    vector<int> new_order(0);
    vector<int> tri_point_lb(P, 0); // low new order index in triangle i
    vector<int> tri_point_ub(P, 0); // high new order index in triangle i
    vector<int> tri_point_count(P, 0); // points in triangle i

    MPI_Win_create(&curr_state[0], 5 * point_count * sizeof(double), sizeof(double), MPI_INFO_NULL, MPI_COMM_WORLD, &win_curr_state);
    MPI_Win_create(&intermediate_1[0], 5 * point_count * sizeof(double), sizeof(double), MPI_INFO_NULL, MPI_COMM_WORLD, &win_inter1);
    MPI_Win_create(&intermediate_2[0], 5 * point_count * sizeof(double), sizeof(double), MPI_INFO_NULL, MPI_COMM_WORLD, &win_inter2);
    MPI_Win_create(&intermediate_3[0], 5 * point_count * sizeof(double), sizeof(double), MPI_INFO_NULL, MPI_COMM_WORLD, &win_inter3);
    MPI_Win_create(&c1234[0], 5 * point_count * sizeof(double), sizeof(double), MPI_INFO_NULL, MPI_COMM_WORLD, &win_c1234);

    // fstream file("../points.csv");
    ifstream file1("../10242points_rh4.csv"); // ifstream = input file stream
    ifstream file2("../10242tris.csv");
    ifstream file3("../10242vert_tris.csv");
    ifstream file4("../10242vert_tri_count.csv");
    string line, word;
    int tri_counts;

    ofstream write_out1("fast_output_mpi.csv", ofstream::out | ofstream::trunc); // ofstream = output file stream
    ofstream write_out2("fast_point_counts_mpi.csv", ofstream::out | ofstream::trunc); // at each time step, write out the number of points


    MPI_File fh;
    MPI_Barrier(MPI_COMM_WORLD);
    MPI_File_open(MPI_COMM_WORLD, "./direct_output_mpi.out", MPI_MODE_WRONLY, MPI_INFO_NULL, &fh);

    int target_loc, local_id, target_points;

    // for (int i = 0; i < point_count; i++) { // sets up particle_thread, in_curr_thread, and own_particles
    //     target_loc = particle_processor(points_per_rank, i, P);
    //     particle_thread[i] = target_loc;
    //     if (target_loc == ID) {
    //         tri_in_curr_thread[i] = true;
    //     }
    // }

    // if (ID == 0) {
    //     for (int i = 0; i < point_count; i++) {
    //         getline(file, line);
    //         stringstream str(line);
    //         for (int j = 0; j < 4; j++) {
    //             getline(str, word, ',');
    //             curr_state[5 * i + j] = stod(word);
    //         }
    //     }
    // }
    // for (int i = 0; i < point_count; i++) {
    //     getline(file1, line);
    //     stringstream str1(line);
    //     for (int j = 0; j < 4; j++) { // read in initial condition of each point
    //         getline(str1, word, ',');
    //         curr_state[5 * i + j] = stod(word);
    //     }
    //
    //     getline(file4, line);
    //     stringstream str4(line);
    //     getline(str4, word, ',');
    //     tri_counts = stod(word); // number of triangles each point borders
    //
    //     vert_tris[i] = vector<int> (tri_counts);
    //     getline(file3, line);
    //     stringstream str3(line);
    //     for (int j = 0; j < tri_counts; j++) { // reads in each points adjacent triangles
    //         getline(str3, word, ',');
    //         vert_tris[i][j] = stod(word);
    //     }
    // }
    if (ID == 0) { // reads in the particles
        vector<double> particle (5, 0);
        vector<double> position;
        double lat;
        for (int i = 0; i < point_count; i++) {
            getline(file1, line);
            stringstream str(line);
            target_loc = particle_thread[i];
            // local_id = i - target_loc * points_per_rank;
            target_points = processor_particle_count(target_loc, P, point_count);
            for (int j = 0; j < 4; j++) {
                getline(str, word, ',');
                particle[j] = stod(word);
            }
            position = slice(particle, 0, 1, 3);
            lat = lat_lon(position)[0];
            particle[4] = lat;
            if (target_loc == 0) {
                // for (int j = 0; j < 5; j++) curr_state[5 * local_id + j] = particle[j];
                for (int j = 0; j < 5; j++) curr_state[5 * i + j] = particle[j];
            } else {
                // cout << "here " << i << endl;
                // MPI_Put(&particle[0], 5, MPI_DOUBLE, target_loc, 5 * local_id, 5, MPI_DOUBLE, win_curr_state);
                MPI_Put(&particle[0], 5, MPI_DOUBLE, target_loc, 5 * i, 5, MPI_DOUBLE, win_curr_state);
            }
        }
        // begin_time = chrono::steady_clock::now();
    }

    for (int i = 0; i < point_count; i++) {
        getline(file4, line);
        stringstream str4(line);
        getline(str4, word, ',');
        tri_counts = stod(word); // number of triangles each point borders

        vert_tris[i] = vector<int> (tri_counts);
        getline(file3, line);
        stringstream str3(line);
        for (int j = 0; j < tri_counts; j++) { // reads in each points adjacent triangles
            getline(str3, word, ',');
            vert_tris[i][j] = stod(word);
        }
    }

    for (int i = 0; i < tri_count; i++) { // reads in triangle vertex information
        getline(file2, line);
        stringstream str2(line);
        for (int j = 0; j < 3; j++) {
            getline(str2, word, ',');
            triangles[i][j] = stod(word);
        }
    }

    MPI_Barrier(MPI_COMM_WORLD);
    // sync_buffer(curr_state, ID, P, lower_bounds, point_counts, &win_curr_state);
    sync_buffer(curr_state, point_count, ID, particle_thread, &win_curr_state);
    MPI_Barrier(MPI_COMM_WORLD);

    // for (int i = 0; i < tri_count; i++) { // reads in triangle vertex information
    //     getline(file2, line);
    //     stringstream str2(line);
    //     for (int j = 0; j < 3; j++) {
    //         getline(str2, word, ',');
    //         triangles[i][j] = stod(word);
    //     }
    // }

    int iv1, iv2, iv3;
    double curr_area;
    vector<double> v1, v2, v3;
    for (int i = 0; i < tri_count; i++) {
        // cout << i << endl;
        iv1 = triangles[i][0];
        iv2 = triangles[i][1];
        iv3 = triangles[i][2];
        // cout << iv1 << " " << iv2 << " " << iv3 << endl;
        v1 = slice(curr_state, 5 * iv1, 1, 3);
        v2 = slice(curr_state, 5 * iv2, 1, 3);
        v3 = slice(curr_state, 5 * iv3, 1, 3);
        curr_area = sphere_tri_area(v1, v2, v3, 1);
        area[iv1] += curr_area / 3.0;
        area[iv2] += curr_area / 3.0;
        area[iv3] += curr_area / 3.0;
    }

    file1.close(); // close all the files we read from
    file2.close();
    file3.close();
    file4.close();

    MPI_Barrier(MPI_COMM_WORLD);


    // vector<double> position;
    // double lat;
    // for (int i = 0; i < point_count; i++) {
    //     position = slice(curr_state, 5 * i, 1, 3);
    //     lat = lat_lon(position)[0];
    //     curr_state[5 * i + 4] = lat; // initial latitude as passive tracer
    // }
    //
    // int iv1, iv2, iv3;
    // double curr_area;
    // vector<double> v1, v2, v3;
    // for (int i = 0; i < tri_count; i++) {
    //     // cout << i << endl;
    //     iv1 = triangles[i][0];
    //     iv2 = triangles[i][1];
    //     iv3 = triangles[i][2];
    //     // cout << iv1 << " " << iv2 << " " << iv3 << endl;
    //     v1 = slice(curr_state, 5 * iv1, 1, 3);
    //     v2 = slice(curr_state, 5 * iv2, 1, 3);
    //     v3 = slice(curr_state, 5 * iv3, 1, 3);
    //     curr_area = sphere_tri_area(v1, v2, v3, 1);
    //     area[iv1] += curr_area / 3.0;
    //     area[iv2] += curr_area / 3.0;
    //     area[iv3] += curr_area / 3.0;
    // }

    MPI_Barrier(MPI_COMM_WORLD);
    // MPI_Win_fence(0, win_curr_state);
    // if (ID != 0) {
    //     MPI_Get(&curr_state[0], 5 * point_count, MPI_DOUBLE, 0, 0, 5 * point_count, MPI_DOUBLE, win_curr_state);
    // }
    // MPI_Win_fence(0, win_curr_state);

    for (int i = 0; i < P; i++) {
        if (i < P - 1) {
            tri_lbs[i] = tris_per_thread * i;
            tri_ubs[i] = tris_per_thread * (i + 1);
        } else {
            tri_lbs[i] = tris_per_thread * i;
            tri_ubs[i] = 20;
        }
    }

    for (int i = 0; i < P; i++) {
        for (int j = tri_lbs[i]; j < tri_ubs[j]; j++) {
            triangle_thread[j] = i;
        }
    }

    MPI_Barrier(MPI_COMM_WORLD);
    icos_init(vertices, triangle_info, triangle_verts, radius, icos_levels-1);
    points_assign(triangle_verts, vertices, curr_state, tri_points, point_locs, icos_levels, point_count);
    // for (int i = 0; i < 20; i++) {
    //
    // }
    // MPI_Barrier(MPI_COMM_WORLD);
    // for (int i = 0; i < point_count; i++) {
    //     int tri = point_locs[0][i];
    //     int thread = triangle_thread[tri];
    //     particle_thread[i] = thread;
    // }
    // MPI_Barrier(MPI_COMM_WORLD);

    if (ID == 0) begin_time = chrono::steady_clock::now();

    if (ID == 0) {
        for (int i = 0; i < point_count; i++) {
            write_out1 << curr_state[5 * i] << "," << curr_state[5 * i + 1] << "," << curr_state[5 * i + 2] << "," << curr_state[5 * i + 3] << "," << curr_state[5 * i + 4] << "," << area[i] << "\n"; // write current state
        }
        write_out2 << point_count << "\n";
        // cout << t << endl;
    }



    // for (int i = 0; i < point_count; i++) {
    //     write_out << curr_state[4 * i] << "," << curr_state[4 * i + 1] << "," << curr_state[4 * i + 2] << "," << curr_state[4 * i + 3] << "\n";
    // }
    //
    // cout << "here" << endl;
    for (int t = 0; t < 1; t++) { // time iterate with RK4
        double curr_time = t * delta_t;
        // points_assign(triangle_verts, vertices, curr_state, tri_points, point_locs, icos_levels, point_count);
        sync_buffer(curr_state, point_count, ID, particle_thread, &win_curr_state);
        MPI_Barrier(MPI_COMM_WORLD);
        BVE_ffunc(c_1, curr_state, vertices, triangle_info, triangle_verts, tri_points, curr_time, delta_t, omega, area, point_count, icos_levels, radius, 0.7, many_count, tri_lb, tri_ub);
        intermediate_1 = c_1;
        scalar_mult(intermediate_1, delta_t / 2);
        vec_add(intermediate_1, curr_state);
        MPI_Barrier(MPI_COMM_WORLD);
        sync_buffer(intermediate_1, point_count, ID, particle_thread, &win_inter1);
        MPI_Barrier(MPI_COMM_WORLD);
    //     // points_assign(triangle_verts, vertices, intermediate_1, tri_points, point_locs, icos_levels, point_count);
        // BVE_ffunc(c_2, intermediate_1, vertices, triangle_info, triangle_verts, tri_points, curr_time + delta_t / 2, delta_t, omega, area, point_count, icos_levels, radius, 0.7, many_count, tri_lb, tri_ub);
        intermediate_2 = c_2;
        scalar_mult(intermediate_2, delta_t / 2);
        vec_add(intermediate_2, curr_state);
        MPI_Barrier(MPI_COMM_WORLD);
        sync_buffer(intermediate_2, point_count, ID, particle_thread, &win_inter2);
        MPI_Barrier(MPI_COMM_WORLD);

    //     // points_assign(triangle_verts, vertices, intermediate_2, tri_points, point_locs, icos_levels, point_count);
        BVE_ffunc(c_3, intermediate_2, vertices, triangle_info, triangle_verts, tri_points, curr_time + delta_t / 2, delta_t, omega, area, point_count, icos_levels, radius, 0.7, many_count, tri_lb, tri_ub);
        intermediate_3 = c_3;
        scalar_mult(intermediate_3, delta_t);
        vec_add(intermediate_3, curr_state);
        MPI_Barrier(MPI_COMM_WORLD);
        sync_buffer(intermediate_3, point_count, ID, particle_thread, &win_inter3);
        MPI_Barrier(MPI_COMM_WORLD);

        // points_assign(triangle_verts, vertices, intermediate_3, tri_points, point_locs, icos_levels, point_count);
        BVE_ffunc(c_4, intermediate_3, vertices, triangle_info, triangle_verts, tri_points, curr_time + delta_t, delta_t, omega, area, point_count, icos_levels, radius, 0.7, many_count, tri_lb, tri_ub);
        c1234 = c_1;
        scalar_mult(c_2, 2);
        vec_add(c1234, c_2);
        scalar_mult(c_3, 2);
        vec_add(c1234, c_3);
        vec_add(c1234, c_4);
        scalar_mult(c1234, delta_t / 6);
        vec_add(c1234, curr_state);
        MPI_Barrier(MPI_COMM_WORLD);
        sync_buffer(c1234, point_count, ID, particle_thread, &win_c1234);
        MPI_Barrier(MPI_COMM_WORLD);

        regrid_point2(c1234, curr_state, triangles, vert_tris, point_count, tri_count, omega, ID, particle_thread);
        cout << "Here 2" << endl;

        // vec_add(curr_state, c1234);
        for (int i = 0; i < point_count; i++) {
            if (particle_thread[i] == ID) {
                vector<double> projected = slice(curr_state, 5 * i, 1, 3);
                project_to_sphere(projected, 1);
                for (int j = 0; j < 3; j++) curr_state[5 * i + j] = projected[j]; // reproject points to surface of sphere
            }
        }
        MPI_Barrier(MPI_COMM_WORLD);
        sync_buffer(curr_state, point_count, ID, particle_thread, &win_curr_state);
        MPI_Barrier(MPI_COMM_WORLD);
        if (ID == 0) {
            for (int i = 0; i < point_count; i++) {
                write_out1 << curr_state[5 * i] << "," << curr_state[5 * i + 1] << "," << curr_state[5 * i + 2] << "," << curr_state[5 * i + 3] << "," << curr_state[5 * i + 4] << "," << area[i] << "\n"; // write current state
            }
            write_out2 << point_count << "\n";
            cout << t << endl;
        }
        MPI_Barrier(MPI_COMM_WORLD);
            // vector<double> projected = slice(curr_state, 4 * i, 1, 3);
            // project_to_sphere(projected, 1);
            // for (int j = 0; j < 3; j++) curr_state[4 * i + j] = projected[j]; // reproject points to surface of sphere
            // write_out << curr_state[4 * i] << "," << curr_state[4 * i + 1] << "," << curr_state[4 * i + 2] << "," << curr_state[4 * i + 3] << "\n"; // write position
        // }
        // cout << t << endl;
    }

    if (ID == 0) {
        chrono::steady_clock::time_point end_time = chrono::steady_clock::now();
        cout << "time taken: " << chrono::duration_cast<chrono::microseconds>(end_time - begin_time).count() << " microseconds" << endl;
    }

    // write_out.close();
    write_out1.close();
    write_out2.close();
    MPI_Finalize();
    return 0;
}
