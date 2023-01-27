#include <iostream>
#include <cmath>
#include <array>
#include <vector>
#include <Accelerate/Accelerate.h>
#include <fstream>
#include <sstream>
#include <chrono>
#include <cassert>

#include "helpers.h"

// time taken: 103913287 microseconds / 103.9 seconds for 2562 particles, 100 time steps, flattened + O3

using namespace std;

void BVE_ffunc(vector<double>& modify, vector<double>& curr_state, double t, double delta_t, double omega, vector<double>& area, int points) {
    vector<double> pos_change, particle_i, particle_j, contribution;
    for (int i = 0; i < points; i++) {
        pos_change = {0, 0, 0};
        particle_i = slice(curr_state, 5 * i, 1, 3);
        for (int j = 0; j < points; j++) {
            if (i != j) {
                particle_j = slice(curr_state, 5 * j, 1, 3);
                contribution = BVE_gfunc(particle_i, particle_j);
                scalar_mult(contribution, curr_state[5 * j + 3] * area[j]);
                vec_add(pos_change, contribution);
            }
        }
        scalar_mult(pos_change, -1.0 / (4.0 * M_PI));
        for (int j = 0; j < 3; j++) modify[5 * i + j] = pos_change[j];
        modify[5 * i + 3] = -2 * omega * pos_change[2];
    }
}

int main() {
    double delta_t = 0.01, end_t = 1;
    int point_count = 2562, tri_count = 5120, time_steps = end_t / delta_t;
    double omega = 2 * M_PI; // coriolis factor

    vector<double> curr_state(5 * point_count); // 0 is x_pos, 1 is y_pos, 2 is z_pos, 3 is vorticity, 4 is passive tracer
    vector<double> c_1(5 * point_count, 0);
    vector<double> c_2(5 * point_count, 0);
    vector<double> c_3(5 * point_count, 0);
    vector<double> c_4(5 * point_count, 0);
    vector<double> c1234(5 * point_count, 0);
    vector<double> intermediate_1(5 * point_count);
    vector<double> intermediate_2(5 * point_count);
    vector<double> intermediate_3(5 * point_count);
    vector<double> area(point_count, 4 * M_PI / point_count); // area associated with each particle

    vector<vector<int>> triangles(tri_count, vector<int> (3)); // triangles[i] is a length 3 int vector containing the indices of the points of the vertices
    vector<vector<int>> vert_tris(point_count); // vert_tris[i] is the int vector containing the indices of the triangles adjacent to point i

    ifstream file1("../points.csv"); // ifstream = input file stream
    ifstream file2("../tris.csv");
    ifstream file3("../vert_tris.csv");
    ifstream file4("../vert_tri_count.csv");
    string line, word;
    int tri_counts;

    ofstream write_out; // ofstream = output file stream
    write_out.open("direct_output.csv", ofstream::out | ofstream::trunc);

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

    chrono::steady_clock::time_point begin = chrono::steady_clock::now();

    for (int i = 0; i < point_count; i++) { // write out initial state
        write_out << curr_state[5 * i] << "," << curr_state[5 * i + 1] << "," << curr_state[5 * i + 2] << "," << curr_state[5 * i + 3] << "," << curr_state[5 * i + 4] << "\n";
    }

    for (int t = 0; t < time_steps; t++) { // time iterate with RK4
    // for (int t = 0; t < 1; t++) {
        double curr_time = t * delta_t;
        BVE_ffunc(c_1, curr_state, curr_time, delta_t, omega, area, point_count);
        intermediate_1 = c_1;
        scalar_mult(intermediate_1, delta_t / 2);
        vec_add(intermediate_1, curr_state);
        BVE_ffunc(c_2, intermediate_1, curr_time + delta_t / 2, delta_t, omega, area, point_count);
        intermediate_2 = c_2;
        scalar_mult(intermediate_2, delta_t / 2);
        vec_add(intermediate_2, curr_state);
        BVE_ffunc(c_3, intermediate_2, curr_time + delta_t / 2, delta_t, omega, area, point_count);
        intermediate_3 = c_3;
        scalar_mult(intermediate_3, delta_t);
        vec_add(intermediate_3, curr_state);
        BVE_ffunc(c_4, intermediate_3, curr_time + delta_t, delta_t, omega, area, point_count);
        c1234 = c_1;
        scalar_mult(c_2, 2);
        vec_add(c1234, c_2);
        scalar_mult(c_3, 2);
        vec_add(c1234, c_3);
        vec_add(c1234, c_4);
        scalar_mult(c1234, delta_t / 6);
        vec_add(c1234, curr_state); // c1234 is new state
        regrid_points(c1234, curr_state, triangles, vert_tris, point_count, tri_count); // regrids points so that they are regular, modifies curr_state
        for (int i = 0; i < point_count; i++) {
            vector<double> projected = slice(curr_state, 5 * i, 1, 3);
            project_to_sphere(projected, 1);
            for (int j = 0; j < 3; j++) curr_state[5 * i + j] = projected[j]; // reproject points to surface of sphere
            write_out << curr_state[5 * i] << "," << curr_state[5 * i + 1] << "," << curr_state[5 * i + 2] << "," << curr_state[5 * i + 3] << "," << curr_state[5 * i + 4] << "\n"; // write current state
        }
        cout << t << endl;
    }

    chrono::steady_clock::time_point end = chrono::steady_clock::now();
    cout << "time taken: " << chrono::duration_cast<chrono::microseconds>(end - begin).count() << " microseconds" << endl;

    write_out.close();
    return 0;
}
