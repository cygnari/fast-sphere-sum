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

// time taken:      274878 microseconds for 2562 particles
// time taken:     4107512 microseconds for 10242 particles
// time taken:    66181297 microseconds for 40962 particles
// time taken:  1064415945 microseconds for 163842 particles
// time taken: 16980804350 microseconds for 655362 particles 

using namespace std;

void BVE_ffunc(vector<double>& modify, vector<double>& curr_state, double t, double delta_t, double omega, double area, int points) {
    for (int i = 0; i < points; i++) {
        vector<double> pos_change {0, 0, 0};
        vector<double> particle_i = slice(curr_state, 4 * i, 1, 3);
        for (int j = 0; j < points; j++) {
            if (i != j) {
                // cout << i << " " << j << endl;
                vector<double> particle_j = slice(curr_state, 4 * j, 1, 3);
                vector<double> contribution = BVE_gfunc(particle_i, particle_j);
                scalar_mult(contribution, curr_state[4 * j + 3] * area);
                // cout << "contribution" << endl;
                vec_add(pos_change, contribution);
                // for (int k = 0; k < 3; k++) cout << contribution[k] << endl;
            }
        }
        scalar_mult(pos_change, -1.0 / (4.0 * M_PI));
        for (int j = 0; j < 3; j++) modify[4 * i + j] = pos_change[j];
        modify[4 * i + 3] = -2 * omega * pos_change[2];
    }
    // cout << "modify" << endl;
    // for (int i = 0; i < modify.size(); i++) cout << modify[i] << " ";
    // cout << endl;
}

#define point_count 655362

int main() {
    double delta_t = 0.01, end_t = 1;
    double omega = 2 * M_PI;
    int time_steps = end_t / delta_t;
    double area = (4 * M_PI) / point_count;

    vector<double> curr_state(4 * point_count); // 0 is x_pos, 1 is y_pos, 2 is z_pos, 3 is vorticity
    vector<double> c_1(4 * point_count, 0);

    fstream file("../points.csv");
    string line, word;

    ofstream write_out;
    write_out.open("direct_output_test.csv", ofstream::out | ofstream::trunc);

    for (int i = 0; i < point_count; i++) {
        getline(file, line);
        stringstream str(line);
        for (int j = 0; j < 4; j++) {
            getline(str, word, ',');
            curr_state[4 * i + j] = stod(word);
        }
    }

    chrono::steady_clock::time_point begin = chrono::steady_clock::now();
    BVE_ffunc(c_1, curr_state, 0, delta_t, omega, area, point_count);
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
