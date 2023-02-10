#include <iostream>
#include <cmath>
#include <array>
#include <vector>
// #include <Accelerate/Accelerate.h>
#include <fstream>
#include <sstream>
#include <queue>
#include <chrono>

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

void BVE_ffunc(vector<double>& modify, vector<double>& curr_state, vector<vector<double>>& vertices, vector<vector<vector<double>>>& tri_info, vector<vector<vector<int>>>& tri_verts, vector<vector<vector<int>>>& tri_points, double t, double delta_t, double omega, vector<double>& area, int points, int levels, double radius, double theta, int many_count, int degree, vector<double>& interp_matrix, vector<vector<double>>& interp_points, int* ipiv) {
    // vector<vector<double>> cluster_points {{1, 0}, {0.5, 0.5}, {0, 1}, {0, 0.5}, {0, 0}, {0.5, 0}};
    // int target_cluster_count = cluster_points.size();
    // vector<vector<double>> cluster_points {{1, 0}, {0.5, 0.5}, {0, 1}, {0, 0.5}, {0, 0}, {0.5, 0}};
    // int source_cluster_count = cluster_points.size();
    int cluster_count = (degree + 1) * (degree + 2) / 2;
    // cout << cluster_count << endl;
    // cout << degree << endl;
    // vector<vector<double>> cluster_points (cluster_count, vector<double> (3, 0));
    // fekete_init(cluster_points, degree);
    int curr_source, curr_target, lev_target, lev_source;
    int particle_count_target, particle_count_source;
    vector<double> center_target, center_source;
    double u, v, us, vs;
    double separation, distance;
    int iv1, iv2, iv3, iv1s, iv2s, iv3s;
    vector<double> v1, v2, v3, v1s, v2s, v3s;
    vector<double> bary_cord;
    int point_index;
    vector<vector<double>> curr_points (cluster_count, vector<double> (3, 0));
    vector<double> interptargets (3 * cluster_count, 0);
    vector<double> placeholder1, placeholder2, placeholder3, source_particle, target_particle;
    // vector<double> Aimatrix {1, 1, 1, 1, 1, 1, 1, 0.5, 0, 0, 0, 0.5, 0, 0.5, 1, 0.5, 0, 0, 0, 0.25, 0, 0, 0, 0, 1, 0.25, 0, 0, 0, 0.25, 0, 0.25, 1, 0.25, 0, 0};
    // vector<double> Aimatrix (cluster_count * cluster_count, 0);
    // interp_mat_init(Aimatrix, cluster_points, degree, cluster_count);
    vector<double> func_vals (3 * cluster_count, 0), func_val (3, 0);
    vector<double> alphas_x (cluster_count, 0), alphas_y (cluster_count, 0), alphas_z (cluster_count, 0);
    // int ipiv[cluster_count];
    // int info, dim = cluster_count, nrhs = 3;
    int nrhs = 3, dim = cluster_count, info;
    char trans = 'N';
    // dgetrf_(&dim, &dim, &*Aimatrix.begin(), &dim, ipiv, &info); // prefactor Ai matrix
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
        // cout << "Here" << endl;
        if ((distance > 0) and (separation < theta)) { // triangles are well separated
            if (particle_count_target > many_count) { // target has many points, do cluster-X interaction
                // cout << "C-X" << endl;
                iv1 = tri_verts[lev_target][curr_target][0];
                iv2 = tri_verts[lev_target][curr_target][1];
                iv3 = tri_verts[lev_target][curr_target][2];
                v1 = vertices[iv1];
                v2 = vertices[iv2];
                v3 = vertices[iv3];
                for (int i = 0; i < cluster_count; i++) {
                    // cout << i << endl;
                    u = interp_points[i][0];
                    // cout
                    v = interp_points[i][1];
                    // cout << "u: " << u << " v: " << v << endl;
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
                // cout << "Here 2" << endl;
                if (particle_count_source > many_count) { // source has many particles, do C-C interaction
                    // cout << "C-C" << endl;
                    iv1s = tri_verts[lev_source][curr_source][0];
                    iv2s = tri_verts[lev_source][curr_source][1];
                    iv3s = tri_verts[lev_source][curr_source][2];
                    v1s = vertices[iv1s];
                    v2s = vertices[iv2s];
                    v3s = vertices[iv3s];
                    for (int i = 0; i < cluster_count; i++) {
                        for (int j = 0; j < cluster_count; j++) {
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
                            for (int k = 0; k < 3; k++) func_vals[j + cluster_count * k] = func_val[k];
                        }
                        dgetrs_(&trans, &dim, &nrhs, &*interp_matrix.begin(), &dim, ipiv, &*func_vals.begin(), &dim, &info);
                        if (info > 0) {
                            cout << info << endl;
                        }

                        for (int j = 0; j < cluster_count; j++) {
                            alphas_x[j] = func_vals[j];
                            alphas_y[j] = func_vals[j + cluster_count];
                            alphas_z[j] = func_vals[j + 2 * cluster_count];
                        }

                        for (int j = 0; j < interptargets.size(); j++) interptargets[j] = 0;

                        for (int j = 0; j < cluster_count; j++) {
                            point_index = tri_points[lev_source][curr_source][j];
                            source_particle = slice(curr_state, 5 * point_index, 1, 3);
                            bary_cord = barycoords(v1s, v2s, v3s, source_particle);
                            interptargets[j] += interp_eval(alphas_x, bary_cord[0], bary_cord[1], degree) * curr_state[5 * point_index + 3] * area[point_index];
                            interptargets[j + cluster_count] += interp_eval(alphas_y, bary_cord[0], bary_cord[1], degree) * curr_state[5 * point_index + 3] * area[point_index];
                            interptargets[j + 2 * cluster_count] += interp_eval(alphas_z, bary_cord[0], bary_cord[1], degree) * curr_state[5 * point_index + 3] * area[point_index];
                        }
                    }
                } else { // source has few particles, do C-P interaction
                    // cout << "C-P" << endl;
                    for (int i = 0; i < cluster_count; i++) {
                        for (int j = 0; j < interptargets.size(); j++) interptargets[j] = 0;

                        for (int j = 0; j < particle_count_source; j++) {
                            point_index = tri_points[lev_source][curr_source][j];
                            source_particle = slice(curr_state, 5 * point_index, 1, 3);
                            func_val = BVE_gfunc(curr_points[i], source_particle);
                            // cout << func_val[0] << endl;
                            // cout << curr_points[0][0] << endl;
                            interptargets[i] += func_val[0] * curr_state[5 * point_index + 3] * area[point_index];
                            interptargets[i + cluster_count] += func_val[1] * curr_state[5 * point_index + 3] * area[point_index];
                            interptargets[i + 2 * cluster_count] += func_val[2] * curr_state[5 * point_index + 3] * area[point_index];
                        }
                    }
                }
                dgetrs_(&trans, &dim, &nrhs, &*interp_matrix.begin(), &dim, ipiv, &*interptargets.begin(), &dim, &info);
                if (info > 0) {
                    cout << info << endl;
                }

                for (int i = 0; i < cluster_count; i++) {
                    alphas_x[i] = interptargets[i];
                    alphas_y[i] = interptargets[i + cluster_count];
                    alphas_z[i] = interptargets[i + 2 * cluster_count];
                }

                for (int i = 0; i < particle_count_target; i++) {
                    point_index = tri_points[lev_target][curr_target][i];
                    target_particle = slice(curr_state, 5 * point_index, 1, 3);
                    bary_cord = barycoords(v1, v2, v3, target_particle);
                    modify[5 * i] += interp_eval(alphas_x, bary_cord[0], bary_cord[1], degree);
                    modify[5 * i + 1] += interp_eval(alphas_y, bary_cord[0], bary_cord[1], degree);
                    modify[5 * i + 2] += interp_eval(alphas_z, bary_cord[0], bary_cord[1], degree);
                }
            } else { // target has few points, do particle-X interaction
                if (particle_count_source > many_count) { // source has lots of points, P-C interaction
                    // cout << "P-C" << endl;
                    iv1s = tri_verts[lev_source][curr_source][0];
                    iv2s = tri_verts[lev_source][curr_source][1];
                    iv3s = tri_verts[lev_source][curr_source][2];
                    v1s = vertices[iv1s];
                    v2s = vertices[iv2s];
                    v3s = vertices[iv3s];
                    for (int i = 0; i < particle_count_target; i++) {
                        point_index = tri_points[lev_target][curr_target][i];
                        target_particle = slice(curr_state, 5 * point_index, 1, 3);
                        for (int j = 0; j < cluster_count; j++) {
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
                            for (int k = 0; k < 3; k++) func_vals[j + cluster_count * k] = func_val[k];
                        }

                        dgetrs_(&trans, &dim, &nrhs, &*interp_matrix.begin(), &dim, ipiv, &*func_vals.begin(), &dim, &info);
                        if (info > 0) {
                            cout << info << endl;
                        }

                        for (int j = 0; j < cluster_count; j++) {
                            alphas_x[j] = func_vals[j];
                            alphas_y[j] = func_vals[j + cluster_count];
                            alphas_z[j] = func_vals[j + 2 * cluster_count];
                        }
                        for (int j = 0; j < particle_count_source; j++) {
                            point_index = tri_points[lev_source][curr_source][j];
                            source_particle = slice(curr_state, 5 * point_index, 1, 3);
                            bary_cord = barycoords(v1s, v2s, v3s, source_particle);
                            modify[5 * i] += interp_eval(alphas_x, bary_cord[0], bary_cord[1], degree) * curr_state[5 * point_index + 3] * area[point_index];
                            modify[5 * i + 1] += interp_eval(alphas_y, bary_cord[0], bary_cord[1], degree) * curr_state[5 * point_index + 3] * area[point_index];
                            modify[5 * i + 2] += interp_eval(alphas_z, bary_cord[0], bary_cord[1], degree) * curr_state[5 * point_index + 3] * area[point_index];
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

// #define point_count 2562

int main() {
    double delta_t = 0.01, end_t = 1; // end_t = number of days
    int point_count = 2562, tri_count = 5120, time_steps = end_t / delta_t, max_points = 1000000;
    // double delta_t = 0.01, end_t = 1;
    double omega = 2 * M_PI; // coriolis
    // int time_steps = end_t / delta_t;
    // double area = (4 * M_PI) / point_count;
    int icos_levels = 2;
    double radius = 1.0;
    double phi = (1 + sqrt(5)) / 2;
    double theta = 0.7;
    int many_count = 10;
    int degree = 2;
    int cluster_count = (degree + 1) * (degree + 2) / 2;

    vector<double> curr_state(5 * point_count); // 0 is x_pos, 1 is y_pos, 2 is z_pos, 3 is vorticity
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

    vector<vector<double>> cluster_points (cluster_count, vector<double> (3, 0)); // interpolation points
    vector<double> interp_matrix (cluster_count * cluster_count, 0); // interpolation matrix

    // cout << "Here" << endl;
    ifstream file1("../2562points_rh4.csv"); // ifstream = input file stream
    ifstream file2("../2562tris.csv");
    ifstream file3("../2562vert_tris.csv");
    ifstream file4("../2562vert_tri_count.csv");
    string line, word;
    int tri_counts;

    ofstream write_out1("direct_output.csv", ofstream::out | ofstream::trunc); // ofstream = output file stream
    ofstream write_out2("direct_point_counts.csv", ofstream::out | ofstream::trunc); // at each time step, write out the number of points
    // write_out1.open("direct_output.csv", ofstream::out | ofstream::trunc);

    for (int i = 0; i < point_count; i++) {
        getline(file1, line);
        stringstream str1(line);
        for (int j = 0; j < 4; j++) { // read in initial condition of each point
            getline(str1, word, ',');
            curr_state[5 * i + j] = stod(word);
        }

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

    file1.close(); // close all the files we read from
    file2.close();
    file3.close();
    file4.close();

    vector<double> position;
    double lat;
    for (int i = 0; i < point_count; i++) {
        position = slice(curr_state, 5 * i, 1, 3);
        lat = lat_lon(position)[0];
        curr_state[5 * i + 4] = lat; // initial latitude as passive tracer
    }

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


    icos_init(vertices, triangle_info, triangle_verts, radius, icos_levels - 1);
    fekete_init(cluster_points, degree);
    interp_mat_init(interp_matrix, cluster_points, degree, cluster_count);
    int ipiv[cluster_count];
    int info = 0, dim = cluster_count;
    cout << cluster_points[0][0] << endl;
    // char trans = 'N';
    dgetrf_(&dim, &dim, &*interp_matrix.begin(), &dim, ipiv, &info); // prefactor Ai matrix
    // cout << info << endl;
    if (info > 0) {
        cout << "dgetrf: " << info << endl;
    }

    // for (int i = 0; i < point_count; i++) {
    //     write_out << curr_state[4 * i] << "," << curr_state[4 * i + 1] << "," << curr_state[4 * i + 2] << "," << curr_state[4 * i + 3] << "\n";
    // }

    // cout << "Here 3" << endl;

    chrono::steady_clock::time_point begin = chrono::steady_clock::now();


    for (int t = 0; t < 1; t++) { // time iterate with RK4
        double curr_time = t * delta_t;
        points_assign(triangle_verts, vertices, curr_state, tri_points, point_locs, icos_levels, point_count);
        // cout << "Here 4" << endl;
        // cout << c_1[0] << endl;
        BVE_ffunc(c_1, curr_state, vertices, triangle_info, triangle_verts, tri_points, curr_time, delta_t, omega, area, point_count, icos_levels, radius, 0.7, many_count, degree, interp_matrix, cluster_points, ipiv);
        // intermediate_1 = c_1;
        // scalar_mult(intermediate_1, delta_t / 2);
        // vec_add(intermediate_1, curr_state);
        // points_assign(triangle_verts, vertices, intermediate_1, tri_points, point_locs, icos_levels, point_count);
        // BVE_ffunc(c_2, intermediate_1, vertices, triangle_info, triangle_verts, tri_points, curr_time + delta_t / 2, delta_t, omega, area, point_count, icos_levels, radius, 0.7, many_count, 3);
        // intermediate_2 = c_2;
        // scalar_mult(intermediate_2, delta_t / 2);
        // vec_add(intermediate_2, curr_state);
        // points_assign(triangle_verts, vertices, intermediate_2, tri_points, point_locs, icos_levels, point_count);
        // BVE_ffunc(c_3, intermediate_2, vertices, triangle_info, triangle_verts, tri_points, curr_time + delta_t / 2, delta_t, omega, area, point_count, icos_levels, radius, 0.7, many_count);
        // intermediate_3 = c_3;
        // scalar_mult(intermediate_3, delta_t);
        // vec_add(intermediate_3, curr_state);
        // points_assign(triangle_verts, vertices, intermediate_3, tri_points, point_locs, icos_levels, point_count);
        // BVE_ffunc(c_4, intermediate_3, vertices, triangle_info, triangle_verts, tri_points, curr_time + delta_t, delta_t, omega, area, point_count, icos_levels, radius, 0.7, many_count);
        // c1234 = c_1;
        // scalar_mult(c_2, 2);
        // vec_add(c1234, c_2);
        // scalar_mult(c_3, 2);
        // vec_add(c1234, c_3);
        // vec_add(c1234, c_4);
        // scalar_mult(c1234, delta_t / 6);
        // vec_add(c1234, curr_state); // c1234 is new state
        // regrid_points(c1234, curr_state, triangles, vert_tris, point_count, tri_count, omega); // regrids points so that they are regular, modifies curr_state
        // state = amr(curr_state, triangles, vert_tris, area, parent_verts, tri_count, point_count, max_points);
        // vec_add(curr_state, c1234);
        // cout << c_1[0] << endl;
        // for (int i = 0; i < point_count; i++) {
        //     // vector<double> projected = slice(curr_state, 4 * i, 1, 3);
        //     // project_to_sphere(projected, 1);
        //     // for (int j = 0; j < 3; j++) curr_state[4 * i + j] = projected[j]; // reproject points to surface of sphere
        //     // write_out << c_1[5 * i] << "," << c_1[5 * i + 1] << "," << c_1[5 * i + 2] << "," << c_1[5 * i + 3] << "\n"; // write position
        //     write_out1 << curr_state[5 * i] << "," << curr_state[5 * i + 1] << "," << curr_state[5 * i + 2] << "," << curr_state[5 * i + 3] << "," << curr_state[5 * i + 4] << "," << area[i] << "\n"; // write current state
        // }
        // write_out2 << point_count << "\n";
        // cout << t << endl;
    }

    chrono::steady_clock::time_point end = chrono::steady_clock::now();
    cout << "time taken: " << chrono::duration_cast<chrono::microseconds>(end - begin).count() << " microseconds" << endl;

    write_out1.close();
    write_out2.close();
    return 0;
}
