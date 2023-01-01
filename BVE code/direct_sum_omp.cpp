#include <iostream>
#include <cmath>
#include <array>
#include <vector>
// #include <Accelerate/Accelerate.h>
#include <fstream>
#include <sstream>
#include <chrono>
#include <cassert>
#include <omp.h>

#include "helpers_omp.h"

// time taken: 90897188 microseconds for 2562 particles, 1 thread
// time taken: 46454253 microseconds for 2562 particles, 2 threads
// time taken: 23875112 microseconds for 2562 particles, 4 threads
// time taken: 14141746 microseconds for 2562 particles, 10 threads

using namespace std;

void BVE_ffunc(vector<double>& modify, vector<double>& curr_state, double t, double delta_t, double omega, double area, int points, int lb, int ub) {
    // cout << "here " << omp_get_thread_num() << endl;
// #pragma omp parallel for num_threads(10)
    for (int i = lb; i < ub; i++) {
        vector<double> pos_change {0, 0, 0};
        vector<double> particle_i = slice(curr_state, 4 * i, 1, 3);
        // cout << "here " << omp_get_thread_num() << endl;
        for (int j = 0; j < points; j++) {
            if (i != j) {
                // cout << "thread: " << omp_get_thread_num() << " i: " << i << " j: " << j << endl;
                vector<double> particle_j = slice(curr_state, 4 * j, 1, 3);
                vector<double> contribution = BVE_gfunc(particle_i, particle_j);
                scalar_mult(contribution, curr_state[4 * j + 3] * area);
                vec_add(pos_change, contribution);
            }
        }
        scalar_mult(pos_change, -1.0 / (4.0 * M_PI));
        for (int j = 0; j < 3; j++) modify[4 * i + j] = pos_change[j];
        // for (int j = 0; j < 3; j++) {
        //     modify[4 * i + j] += pos_change[j];
        //     cout << pos_change[j] << endl;
        // }
        modify[4 * i + 3] = -2 * omega * pos_change[2];
    }
}

#define point_count 163842

int main() {
    double delta_t = 0.01, end_t = 1;
    double omega = 2 * M_PI;
    int time_steps = end_t / delta_t;
    double area = (4 * M_PI) / point_count;

    vector<double> curr_state(4 * point_count); // 0 is x_pos, 1 is y_pos, 2 is z_pos, 3 is vorticity
    vector<double> c_1(4 * point_count);
    vector<double> c_2(4 * point_count);
    vector<double> c_3(4 * point_count);
    vector<double> c_4(4 * point_count);
    // vector<double> c12(4 * point_count);
    // vector<double> c123(4 * point_count);
    vector<double> c1234(4 * point_count);
    vector<double> intermediate_1(4 * point_count);
    vector<double> intermediate_2(4 * point_count);
    vector<double> intermediate_3(4 * point_count);

    // fstream file("../points.csv");
    fstream file("./points.csv");
    string line, word;

    ofstream write_out;
    write_out.open("direct_output_omp.csv", ofstream::out | ofstream::trunc);

    for (int i = 0; i < point_count; i++) {
        getline(file, line);
        stringstream str(line);
        for (int j = 0; j < 4; j++) {
            getline(str, word, ',');
            curr_state[4 * i + j] = stod(word);
        }
    }

    chrono::steady_clock::time_point begin = chrono::steady_clock::now();

    // write out initial state
    // for (int i = 0; i < point_count; i++) {
    //     write_out << curr_state[4 * i] << "," << curr_state[4 * i + 1] << "," << curr_state[4 * i + 2] << "," << curr_state[4 * i + 3] << "\n";
    // }
#pragma omp parallel
{ // begin parallel
    int thread_count = omp_get_num_threads();
    int thread_id = omp_get_thread_num();

    int points_per_thread = point_count / thread_count;
    int own_points, lb, ub;

    if (thread_id < thread_count - 1) {
        own_points = points_per_thread;
        lb = thread_id * points_per_thread;
        ub = (thread_id + 1) * points_per_thread;
    } else {
        own_points = point_count - (thread_count - 1) * points_per_thread;
        lb = point_count - own_points;
        ub = point_count;
    }
#pragma omp critical
{
    cout << "thread count " << thread_count << " thread id " << thread_id << " lb " << lb << " ub " << ub << " own points " << own_points << endl;
}
#pragma omp barrier
// { // begin parallel
    for (int t = 0; t < 1; t++) { // time iterate with RK4
#pragma omp barrier
        double curr_time = t * delta_t;
        BVE_ffunc(c_1, curr_state, curr_time, delta_t, omega, area, point_count, lb, ub);
#pragma omp barrier
        if (thread_id == 0) intermediate_1 = c_1;
#pragma omp barrier
        vec_scalar_mult_add_omp(intermediate_1, curr_state, delta_t / 2, 4 * lb, 4 * ub);
#pragma omp barrier
        BVE_ffunc(c_2, intermediate_1, curr_time + delta_t / 2, delta_t, omega, area, point_count, lb, ub);
#pragma omp barrier
        if (thread_id == 0) intermediate_2 = c_2;
#pragma omp barrier
        vec_scalar_mult_add_omp(intermediate_2, curr_state, delta_t / 2, 4 * lb, 4 * ub);
#pragma omp barrier
        BVE_ffunc(c_3, intermediate_2, curr_time + delta_t / 2, delta_t, omega, area, point_count, lb, ub);
#pragma omp barrier
        if (thread_id == 0) intermediate_3 = c_3;
#pragma omp barrier
        vec_scalar_mult_add_omp(intermediate_3, curr_state, delta_t / 2, 4 * lb, 4 * ub);
#pragma omp barrier
        BVE_ffunc(c_4, intermediate_3, curr_time + delta_t, delta_t, omega, area, point_count, lb, ub);
#pragma omp barrier
        if (thread_id == 0) c1234 = c_1;
#pragma omp barrier
        vec_add_scalar_mult_omp(c1234, c_2, 2, 4 * lb, 4 * ub);
        vec_add_scalar_mult_omp(c1234, c_3, 2, 4 * lb, 4 * ub);
        vec_add_omp(c1234, c_4, 4 * lb, 4 * ub);
        vec_add_scalar_mult_omp(curr_state, c1234, delta_t / 6, 4 * lb, 4 * ub);
#pragma omp barrier
        for (int i = lb; i < ub; i++) {
            // cout << "i: " << i << endl;
            vector<double> projected = slice(curr_state, 4 * i, 1, 3);
            // projected = project_to_sphere(projected, 1);
            project_to_sphere(projected, 1);
            for (int j = 0; j < 3; j++) curr_state[4 * i + j] = projected[j]; // reproject points to surface of sphere
            // write_out << curr_state[4 * i] << "," << curr_state[4 * i + 1] << "," << curr_state[4 * i + 2] << "," << curr_state[4 * i + 3] << "\n"; // write position
        }
// #pragma omp barrier
// #pragma omp single
// { // thread 0 writes out
//         for (int i = 0; i < point_count; i++) {
//             write_out << curr_state[4 * i] << "," << curr_state[4 * i + 1] << "," << curr_state[4 * i + 2] << "," << curr_state[4 * i + 3] << "\n"; // write position
//         }
//         cout << t << endl;
// }
#pragma omp barrier
    }
} // end parallel

    chrono::steady_clock::time_point end = chrono::steady_clock::now();
    cout << "time taken: " << chrono::duration_cast<chrono::microseconds>(end - begin).count() << " microseconds" << endl;


    write_out.close();
    return 0;
}
