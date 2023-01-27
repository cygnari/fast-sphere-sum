#include <iostream>
#include <cmath>
#include <array>
#include <vector>
// #include <Accelerate/Accelerate.h>
#include <fstream>
#include <sstream>
#include <chrono>
#include <cassert>

#include "helpers.h"

// time taken: 103913287 microseconds / 103.9 seconds for 2562 particles, 100 time steps, flattened + O3

using namespace std;

void BVE_ffunc(vector<double>& modify, vector<double>& curr_state, double t, double delta_t, double omega, vector<double>& area, int points) {
    for (int i = 0; i < points; i++) {
        vector<double> pos_change {0, 0, 0};
        vector<double> particle_i = slice(curr_state, 5 * i, 1, 3);
        for (int j = 0; j < points; j++) {
            if (i != j) {
                // cout << i << " " << j << endl;
                vector<double> particle_j = slice(curr_state, 5 * j, 1, 3);
                vector<double> contribution = BVE_gfunc(particle_i, particle_j);
                scalar_mult(contribution, curr_state[5 * j + 3] * area[j]);
                // cout << "contribution" << endl;
                vec_add(pos_change, contribution);
                // for (int k = 0; k < 3; k++) cout << contribution[k] << endl;
            }
        }
        scalar_mult(pos_change, -1.0 / (4.0 * M_PI));
        for (int j = 0; j < 3; j++) modify[5 * i + j] = pos_change[j];
        // for (int j = 0; j < 3; j++) {
        //     modify[4 * i + j] += pos_change[j];
        //     cout << pos_change[j] << endl;
        // }
        modify[5 * i + 3] = -2 * omega * pos_change[2];
    }
    // cout << "modify" << endl;
    // for (int i = 0; i < modify.size(); i++) cout << modify[i] << " ";
    // cout << endl;
}

// #define point_count 2562
// #define tri_count 5120

int main() {
    int point_count = 2562;
    int tri_count = 5120;
    double delta_t = 0.01, end_t = 1;
    double omega = 2 * M_PI;
    int time_steps = end_t / delta_t;
    // double area = (4 * M_PI) / point_count;

    vector<double> curr_state(5 * point_count); // 0 is x_pos, 1 is y_pos, 2 is z_pos, 3 is vorticity, 4 is passive tracer
    vector<double> c_1(5 * point_count, 0);
    vector<double> c_2(5 * point_count, 0);
    vector<double> c_3(5 * point_count, 0);
    vector<double> c_4(5 * point_count, 0);
    // vector<double> c12(4 * point_count);
    // vector<double> c123(4 * point_count);
    vector<double> c1234(5 * point_count, 0);
    vector<double> intermediate_1(5 * point_count);
    vector<double> intermediate_2(5 * point_count);
    vector<double> intermediate_3(5 * point_count);
    vector<double> area(point_count, 4 * M_PI / point_count);

    vector<vector<int>> triangles(tri_count, vector<int> (3));
    vector<vector<int>> vert_tris(point_count);

    fstream file1("../points.csv");
    fstream file2("../tris.csv");
    fstream file3("../vert_tris.csv");
    fstream file4("../vert_tri_count.csv");
    // fstream file1("./points.csv");
    string line, word;
    int tri_counts;

    ofstream write_out;
    write_out.open("direct_output.csv", ofstream::out | ofstream::trunc);

    for (int i = 0; i < point_count; i++) {
        getline(file1, line);
        stringstream str1(line);
        for (int j = 0; j < 4; j++) {
            getline(str1, word, ',');
            curr_state[5 * i + j] = stod(word);
        }

        // getline(file2, line);
        // stringstream str2(line);
        // for (int j = 0; j < 3; j++) {
        //     getline(str2, word, ',');
        //     triangles[i][j] = stod(word);
        // }

        getline(file4, line);
        stringstream str4(line);
        getline(str4, word, ',');
        tri_counts = stod(word);

        vert_tris[i] = vector<int> (tri_counts);
        getline(file3, line);
        stringstream str3(line);
        for (int j = 0; j < tri_counts; j++) {
            getline(str3, word, ',');
            vert_tris[i][j] = stod(word);
        }
    }

    // cout << "here" << endl;

    for (int i = 0; i < tri_count; i++) {
        getline(file2, line);
        stringstream str2(line);
        for (int j = 0; j < 3; j++) {
            getline(str2, word, ',');
            triangles[i][j] = stod(word);
        }
    }

    vector<double> position;
    for (int i = 0; i < point_count; i++) {
        position = slice(curr_state, 5 * i, 1, 3);
    }

    // cout << vert_tris[1][4] << endl;

    // cout << "here" << endl;

    chrono::steady_clock::time_point begin = chrono::steady_clock::now();

    // write out initial state
    for (int i = 0; i < point_count; i++) {
        write_out << curr_state[5 * i] << "," << curr_state[5 * i + 1] << "," << curr_state[5 * i + 2] << "," << curr_state[5 * i + 3] << "\n";
    }

    // for (int t = 0; t < time_steps; t++) { // time iterate with RK4
    for (int t = 0; t < 1; t++) {
        double curr_time = t * delta_t;
        // for (int i = 0; i < point_count; i++) {
        //     // cout << intermediate_1[4 * i] << "," << intermediate_1[4 * i + 1] << "," << intermediate_1[4 * i + 2] << "," << intermediate_1[4 * i + 3] << endl;
        //     write_out << "curr," << curr_state[4 * i] << "," << curr_state[4 * i + 1] << "," << curr_state[4 * i + 2] << "," << curr_state[4 * i + 3] << "\n"; // write position
        // }
        // cout << "1" << endl;
        BVE_ffunc(c_1, curr_state, curr_time, delta_t, omega, area, point_count);
        intermediate_1 = c_1;
        // cout << "intermediate 1" << endl;
        // for (int i = 0; i < point_count; i++) {
        //     // cout << intermediate_1[4 * i] << "," << intermediate_1[4 * i + 1] << "," << intermediate_1[4 * i + 2] << "," << intermediate_1[4 * i + 3] << endl;
        //     write_out << "inter 1," << intermediate_1[4 * i] << "," << intermediate_1[4 * i + 1] << "," << intermediate_1[4 * i + 2] << "," << intermediate_1[4 * i + 3] << "\n"; // write position
        // }
        // cout << "2" << endl;

        scalar_mult(intermediate_1, delta_t / 2);
        // cout << "intermediate 1" << endl;
        // for (int i = 0; i < point_count; i++) {
        //     // cout << intermediate_1[4 * i] << "," << intermediate_1[4 * i + 1] << "," << intermediate_1[4 * i + 2] << "," << intermediate_1[4 * i + 3] << endl;
        //     write_out << "inter 1," << intermediate_1[4 * i] << "," << intermediate_1[4 * i + 1] << "," << intermediate_1[4 * i + 2] << "," << intermediate_1[4 * i + 3] << "\n"; // write position
        // }
        vec_add(intermediate_1, curr_state);

        BVE_ffunc(c_2, intermediate_1, curr_time + delta_t / 2, delta_t, omega, area, point_count);
        intermediate_2 = c_2;
        // cout << "3" << endl;

        // cout << "intermediate 2" << endl;
        // for (int i = 0; i < point_count; i++) {
        //     // cout << intermediate_1[4 * i] << "," << intermediate_1[4 * i + 1] << "," << intermediate_1[4 * i + 2] << "," << intermediate_1[4 * i + 3] << endl;
        //     write_out << "inter 2," << intermediate_2[4 * i] << "," << intermediate_2[4 * i + 1] << "," << intermediate_2[4 * i + 2] << "," << intermediate_2[4 * i + 3] << "\n"; // write position
        // }
        scalar_mult(intermediate_2, delta_t / 2);
        // cout << "intermediate 2" << endl;
        // for (int i = 0; i < point_count; i++) {
        //     // cout << intermediate_1[4 * i] << "," << intermediate_1[4 * i + 1] << "," << intermediate_1[4 * i + 2] << "," << intermediate_1[4 * i + 3] << endl;
        //     write_out << "inter 2," << intermediate_2[4 * i] << "," << intermediate_2[4 * i + 1] << "," << intermediate_2[4 * i + 2] << "," << intermediate_2[4 * i + 3] << "\n"; // write position
        // }
        vec_add(intermediate_2, curr_state);
        // for (int i = 0; i < point_count; i++) {
        //     cout << intermediate_2[4 * i] << "," << intermediate_2[4 * i + 1] << "," << intermediate_2[4 * i + 2] << "," << intermediate_2[4 * i + 3] << endl;
        // }
        BVE_ffunc(c_3, intermediate_2, curr_time + delta_t / 2, delta_t, omega, area, point_count);
        intermediate_3 = c_3;
        // cout << "4" << endl;

        scalar_mult(intermediate_3, delta_t);
        vec_add(intermediate_3, curr_state);
        BVE_ffunc(c_4, intermediate_3, curr_time + delta_t, delta_t, omega, area, point_count);
        // cout << "5" << endl;

        c1234 = c_1;
        scalar_mult(c_2, 2);
        vec_add(c1234, c_2);
        scalar_mult(c_3, 2);
        vec_add(c1234, c_3);
        vec_add(c1234, c_4);
        scalar_mult(c1234, delta_t / 6);
        // vec_add(curr_state, c1234);
        vec_add(c1234, curr_state); // c1234 is new state
        regrid_points(c1234, curr_state, triangles, vert_tris, point_count, tri_count);
        for (int i = 0; i < point_count; i++) {
            vector<double> projected = slice(curr_state, 5 * i, 1, 3);
            // projected = project_to_sphere(projected, 1);
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
