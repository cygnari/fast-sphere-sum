#include <iostream>
#include <cmath>
#include <array>
#include <vector>
// #include <Accelerate/Accelerate.h>
#include <fstream>
#include <sstream>
#include <chrono>
#include <cassert>

#include "helpers_cuda.h"
#include "helpers.h"

// time taken: 103913287 microseconds / 103.9 seconds for 2562 particles, 100 time steps, flattened + O3

using namespace std;

__global__ void BVE_ffunc(double* modify, double* curr_state, double t, double delta_t, double omega, double area, int points) {
    int id = blockIdx.x*blockDim.x + threadIdx.x;
    if (id < points) {
        // vector<double> pos_change {0, 0, 0};
        // double[3] pos_change = {0, 0, 0};
        double pos_change[3] = {0, 0, 0};
        // vector<double> particle_i {0, 0, 0};
        // double[3] particle_i;
        double particle_i[3];
        slice_cuda(curr_state, particle_i, 4 * id, 1, 3);
        for (int j = 0; j < points; j++) {
            if (id != j) {
                // vector<double> particle_j = slice2(curr_state, 4 * j, 1, 3);
                // double[3] particle_j;
                double particle_j[3];
                // vector<double> particle_j {0, 0, 0};
                slice_cuda(curr_state, particle_j, 4 * j, 1, 3);
                // double[3] contribution;
                double contribution[3];
                BVE_gfunc_cuda(particle_i, particle_j, contribution);
                scalar_mult_cuda_dev(contribution, curr_state[4 * j + 3] * area, 3);
                vec_add_cuda_dev(pos_change, contribution, 3);
            }
        }
        scalar_mult_cuda_dev(pos_change, -1.0 / (4.0 * M_PI), 3);
        for (int j = 0; j < 3; j++) modify[4 * id + j] = pos_change[j];
        modify[4 * id + 3] = -2 * omega * pos_change[2];
    }
}

__global__ void projection(double* curr_state, double radius, int points) {
    int id = blockIdx.x*blockDim.x + threadIdx.x;
    if (id < points) {
        // vector<double> projected {0, 0, 0};
        // double[3] projected;
        double projected[3];
        slice_cuda(curr_state, projected, 4 * id, 1, 3);
        project_to_sphere_cuda(projected, radius);
        for (int i = 0; i < 3; i++) curr_state[4 * id + i] = projected[i];
    }
}

#define point_count 655362

int main() {
    double delta_t = 0.01, end_t = 1;
    double omega = 2 * M_PI;
    int time_steps = end_t / delta_t;
    double area = (4 * M_PI) / point_count;

    vector<double> curr_state(4 * point_count); // 0 is x_pos, 1 is y_pos, 2 is z_pos, 3 is vorticity
    vector<double> c_1(4 * point_count, 0);
    vector<double> c_2(4 * point_count, 0);
    vector<double> c_3(4 * point_count, 0);
    vector<double> c_4(4 * point_count, 0);
    vector<double> c1234(4 * point_count, 0);
    vector<double> intermediate_1(4 * point_count);
    vector<double> intermediate_2(4 * point_count);
    vector<double> intermediate_3(4 * point_count);

    // fstream file("../points.csv");
    fstream file("./points.csv");
    string line, word;

    ofstream write_out;
    write_out.open("direct_output.csv", ofstream::out | ofstream::trunc);

    for (int i = 0; i < point_count; i++) {
        getline(file, line);
        stringstream str(line);
        for (int j = 0; j < 4; j++) {
            getline(str, word, ',');
            curr_state[4 * i + j] = stod(word);
        }
    }

    double *d_curr, *d_c1, *d_c2, *d_c3, *d_c4, *d_c1234, *d_inter1, *d_inter2, *d_inter3;

    cudaMalloc(&d_curr, point_count*sizeof(double));
    cudaMalloc(&d_c1, point_count*sizeof(double));
    cudaMalloc(&d_c2, point_count*sizeof(double));
    cudaMalloc(&d_c3, point_count*sizeof(double));
    cudaMalloc(&d_c4, point_count*sizeof(double));
    cudaMalloc(&d_c1234, point_count*sizeof(double));
    cudaMalloc(&d_inter1, point_count*sizeof(double));
    cudaMalloc(&d_inter2, point_count*sizeof(double));
    cudaMalloc(&d_inter3, point_count*sizeof(double));

    cudaMemcpy(d_curr, &curr_state[0], point_count*sizeof(double), cudaMemcpyHostToDevice);
    cudaMemcpy(d_c1, &c_1[0], point_count*sizeof(double), cudaMemcpyHostToDevice);
    cudaMemcpy(d_c2, &c_2[0], point_count*sizeof(double), cudaMemcpyHostToDevice);
    cudaMemcpy(d_c3, &c_3[0], point_count*sizeof(double), cudaMemcpyHostToDevice);
    cudaMemcpy(d_c4, &c_4[0], point_count*sizeof(double), cudaMemcpyHostToDevice);
    cudaMemcpy(d_c1234, &c1234[0], point_count*sizeof(double), cudaMemcpyHostToDevice);
    cudaMemcpy(d_inter1, &intermediate_1[0], point_count*sizeof(double), cudaMemcpyHostToDevice);
    cudaMemcpy(d_inter2, &intermediate_2[0], point_count*sizeof(double), cudaMemcpyHostToDevice);
    cudaMemcpy(d_inter3, &intermediate_3[0], point_count*sizeof(double), cudaMemcpyHostToDevice);

    chrono::steady_clock::time_point begin = chrono::steady_clock::now();

    // write out initial state
    // for (int i = 0; i < point_count; i++) {
    //     write_out << curr_state[4 * i] << "," << curr_state[4 * i + 1] << "," << curr_state[4 * i + 2] << "," << curr_state[4 * i + 3] << "\n";
    // }

    for (int t = 0; t < 1; t++) { // time iterate with RK4
        double curr_time = t * delta_t;
        cudaDeviceSynchronize();
        BVE_ffunc<<<(point_count + 255) / 256, 256>>>(d_c1, d_curr, curr_time, delta_t, omega, area, point_count);
        cudaDeviceSynchronize();
        copy_cuda<<<(point_count + 255) / 256, 256>>>(d_inter1, d_c1, point_count);
        cudaDeviceSynchronize();
        scalar_mult_cuda<<<(point_count + 255) / 256, 256>>>(d_inter1, delta_t / 2, point_count);
        cudaDeviceSynchronize();
        vec_add_cuda<<<(point_count + 255) / 256, 256>>>(d_inter1, d_curr, point_count);
        cudaDeviceSynchronize();
        BVE_ffunc<<<(point_count + 255) / 256, 256>>>(d_c2, d_inter1, curr_time + delta_t / 2, delta_t, omega, area, point_count);
        cudaDeviceSynchronize();
        copy_cuda<<<(point_count + 255) / 256, 256>>>(d_inter2, d_c2, point_count);
        cudaDeviceSynchronize();
        scalar_mult_cuda<<<(point_count + 255) / 256, 256>>>(d_inter2, delta_t / 2, point_count);
        cudaDeviceSynchronize();
        vec_add_cuda<<<(point_count + 255) / 256, 256>>>(d_inter2, d_curr, point_count);
        cudaDeviceSynchronize();
        BVE_ffunc<<<(point_count + 255) / 256, 256>>>(d_c3, d_inter2, curr_time + delta_t / 2, delta_t, omega, area, point_count);
        cudaDeviceSynchronize();
        copy_cuda<<<(point_count + 255) / 256, 256>>>(d_inter3, d_c3, point_count);
        cudaDeviceSynchronize();
        scalar_mult_cuda<<<(point_count + 255) / 256, 256>>>(d_inter3, delta_t, point_count);
        cudaDeviceSynchronize();
        vec_add_cuda<<<(point_count + 255) / 256, 256>>>(d_inter3, d_curr, point_count);
        cudaDeviceSynchronize();
        BVE_ffunc<<<(point_count + 255) / 256, 256>>>(d_c4, d_inter3, curr_time + delta_t, delta_t, omega, area, point_count);
        cudaDeviceSynchronize();
        copy_cuda<<<(point_count + 255) / 256, 256>>>(d_c1234, d_c1, point_count);
        cudaDeviceSynchronize();
        scalar_mult_cuda<<<(point_count + 255) / 256, 256>>>(d_c2, 2, point_count);
        cudaDeviceSynchronize();
        vec_add_cuda<<<(point_count + 255) / 256, 256>>>(d_c1234, d_c2, point_count);
        cudaDeviceSynchronize();
        scalar_mult_cuda<<<(point_count + 255) / 256, 256>>>(d_c3, 2, point_count);
        cudaDeviceSynchronize();
        vec_add_cuda<<<(point_count + 255) / 256, 256>>>(d_c1234, d_c3, point_count);
        cudaDeviceSynchronize();
        vec_add_cuda<<<(point_count + 255) / 256, 256>>>(d_c1234, d_c4, point_count);
        cudaDeviceSynchronize();
        scalar_mult_cuda<<<(point_count + 255) / 256, 256>>>(d_c1234, delta_t / 6, point_count);
        cudaDeviceSynchronize();
        vec_add_cuda<<<(point_count + 255) / 256, 256>>>(d_curr, d_c1234, point_count);
        cudaDeviceSynchronize();
        projection<<<(point_count + 255) / 256, 256>>>(d_curr, 1.0, point_count);
        // for (int i = 0; i < point_count; i++) {
        //     vector<double> projected = slice(curr_state, 4 * i, 1, 3);
        //     // projected = project_to_sphere(projected, 1);
        //     project_to_sphere(projected, 1);
        //     for (int j = 0; j < 3; j++) curr_state[4 * i + j] = projected[j]; // reproject points to surface of sphere
        //     // write_out << curr_state[4 * i] << "," << curr_state[4 * i + 1] << "," << curr_state[4 * i + 2] << "," << curr_state[4 * i + 3] << "\n"; // write position
        // }
        // cout << t << endl;
    }

    chrono::steady_clock::time_point end = chrono::steady_clock::now();
    cout << "time taken: " << chrono::duration_cast<chrono::microseconds>(end - begin).count() << " microseconds" << endl;

    cudaFree(d_curr);
    cudaFree(d_c1);
    cudaFree(d_c2);
    cudaFree(d_c3);
    cudaFree(d_c4);
    cudaFree(d_c1234);
    cudaFree(d_inter1);
    cudaFree(d_inter2);
    cudaFree(d_inter3);

    write_out.close();
    return 0;
}
