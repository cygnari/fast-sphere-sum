#include <iostream>
#include <cmath>
#include <array>
#include <vector>
// #include <Accelerate/Accelerate.h>
#include <fstream>
#include <sstream>
#include <chrono>
#include <cassert>
#include <mpi.h>

#include "helpers.h"

// time taken: 103913287 microseconds / 103.9 seconds for 2562 particles, 100 time steps, flattened + O3

using namespace std;

void BVE_ffunc(vector<double>& modify, vector<double>& curr_state, double t, double delta_t, double omega, vector<double>& area, int total_points, int lb, int ub) {

    for (int i = lb; i < ub; i++) {
        vector<double> pos_change (3, 0);
        vector<double> particle_i = slice(curr_state, 5 * i, 1, 3);
        // cout << "particle i " << particle_i[0] << " " << particle_i[1] << " " << particle_i[2] << endl;
        for (int j = 0; j < total_points; j++) {
            // int global_id = own_particles[i];
            if (i != j) {
                // cout << own_particles[i] << " " << j << endl;
                // int local_id_j = local_ids[j];
                vector<double> particle_j = slice(curr_state, 5 * j, 1, 3);
                // if (in_curr_thread[j]) {
                //     particle_j = slice(curr_state, 4 * local_id_j, 1, 4);
                // } else {
                //     int source_thread = particle_thread[j];
                //     int status = MPI_Get(&particle_j[0], 4, MPI_DOUBLE, source_thread, 4 * local_id_j, 4, MPI_DOUBLE, window);
                //     MPI_Win_fence(0, window);
                // }
                // vector<double> pos_j = slice(particle_j, 0, 1, 3);

                vector<double> contribution = BVE_gfunc(particle_i, particle_j);
                scalar_mult(contribution, curr_state[5 * j + 3] * area[j]);
                // cout << "contribution" << endl;
                vec_add(pos_change, contribution);
                // for (int k = 0; k < 3; k++) cout << contribution[k] << endl;
            }
        }
        scalar_mult(pos_change, -1.0 / (4.0 * M_PI));
        for (int j = 0; j < 3; j++) modify[5 * i + j] = pos_change[j];
        modify[5 * i + 3] = -2 * omega * pos_change[2];
    }
    // MPI_Barrier(MPI_COMM_WORLD);
    // MPI_Win_fence(0, window);
    // cout << "modify" << endl;
    // for (int i = 0; i < modify.size(); i++) cout << modify[i] << " ";
    // cout << endl;
}

void sync_buffer(vector<double>& buffer, int ID, int P, vector<int>& lb, vector<int>& counts, MPI_Win *win) {
    MPI_Win_fence(0, *win);
    // cout << buffer.size() << endl;
    // cout << "buff 8 " << buffer[7] << endl;
    for (int i = 0; i < P; i++) {
        cout << i << endl;

        if (i != ID) {
            // cout << "i: " << i << " ID: " << ID << " P: " << P << " lb: " << lb[i] << " num: " << counts[i] << endl;
            MPI_Get(&buffer[5 * lb[i]], 5 * counts[i], MPI_DOUBLE, i, 5 * lb[i], 5 * counts[i], MPI_DOUBLE, *win);
            // MPI_Get(&buffer[0], 4 * counts[i], MPI_DOUBLE, i, 0, 4 * counts[i], MPI_DOUBLE, *win);
            // cout << "id: " << ID << " success" << endl;
        }
    }
    MPI_Win_fence(0, *win);
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

// #define point_count 163842

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

    double delta_t = 0.01, end_t = 1; // end_t = number of days
    int point_count = 10242, tri_count = 20480, time_steps = end_t / delta_t, max_points = 1000000;

    chrono::steady_clock::time_point begin_time;

    int points_per_rank = point_count / P;
    int own_points;

    if (ID < P - 1) {
        own_points = points_per_rank;
    } else {
        own_points = point_count - (P - 1) * points_per_rank;
    }
    cout << "P: " << P << " ID: " << ID << " points_per_rank: " << points_per_rank << " own_points: " << own_points << endl;
    // cout << "P: " << P << " ID: " << ID << endl;


    // int point_count = 40962, tri_count = 81920, time_steps = end_t / delta_t, max_points = 1000000;
    // double delta_t = 0.01, end_t = 1;
    double omega = 2 * M_PI;
    // int time_steps = end_t / delta_t;
    // double area = (4 * M_PI) / point_count;

    vector<double> curr_state(5 * point_count, 0); // 0 is x_pos, 1 is y_pos, 2 is z_pos, 3 is vorticity
    vector<double> c_1(5 * point_count, 0);
    vector<double> c_2(5 * point_count, 0);
    vector<double> c_3(5 * point_count, 0);
    vector<double> c_4(5 * point_count, 0);
    vector<double> c1234(5 * point_count, 0);
    vector<double> intermediate_1(5 * point_count, 0);
    vector<double> intermediate_2(5 * point_count, 0);
    vector<double> intermediate_3(5 * point_count, 0);
    vector<double> area(point_count, 0);
    vector<int> state;

    vector<double> particle_thread(point_count, 0); // particle_thread[i] is the processor where particle i is located
    // vector<bool> in_curr_thread(point_count, false); // true if particle i is in the current thread
    // vector<int> own_particles(own_points); // own_particles[i] is the global id of the particle
    // vector<int> local_ids(point_count); // local id of particle i
    vector<int> lower_bounds (P, 0); // lower bound of each processor particles
    vector<int> upper_bounds (P, 0); // upper bound of each processor particles
    vector<int> point_counts (P, 0); // points in each processor

    MPI_Win_create(&curr_state[0], 5 * point_count * sizeof(double), sizeof(double), MPI_INFO_NULL, MPI_COMM_WORLD, &win_curr_state);
    MPI_Win_create(&intermediate_1[0], 5 * point_count * sizeof(double), sizeof(double), MPI_INFO_NULL, MPI_COMM_WORLD, &win_inter1);
    MPI_Win_create(&intermediate_2[0], 5 * point_count * sizeof(double), sizeof(double), MPI_INFO_NULL, MPI_COMM_WORLD, &win_inter2);
    MPI_Win_create(&intermediate_3[0], 5 * point_count * sizeof(double), sizeof(double), MPI_INFO_NULL, MPI_COMM_WORLD, &win_inter3);

    MPI_Win_fence(0, win_curr_state);
    MPI_File fh;

    // fstream file("../points.csv");
    // fstream file("./points.csv");
    // string line, word;
    ifstream file1("../10242points_rh4.csv"); // ifstream = input file stream
    ifstream file2("../10242tris.csv");
    ifstream file3("../10242vert_tris.csv");
    ifstream file4("../10242vert_tri_count.csv");
    string line, word;
    int tri_counts;

    // ofstream write_out;
    ofstream write_out1("direct_output_mpi.csv", ofstream::out | ofstream::trunc); // ofstream = output file stream
    ofstream write_out2("direct_point_counts_mpi.csv", ofstream::out | ofstream::trunc); // at each time step, write out the number of points

    // write_out1.open("./direct_output_mpi.out", ofstream::out | ofstream::trunc);
    write_out1.close();
    write_out2.close();

    MPI_Barrier(MPI_COMM_WORLD);
    // cout << "here" << endl;
    // MPI_File_open(MPI_COMM_WORLD, "./direct_output_mpi.out", MPI_MODE_WRONLY, MPI_INFO_NULL, &fh);
    // cout << "here" << endl;

    int target_loc, local_id, target_points;

    for (int i = 0; i < point_count; i++) { // sets up particle_thread, in_curr_thread, and own_particles
        target_loc = particle_processor(points_per_rank, i, P);
        // local_id = i - target_loc * points_per_rank;
        // target_points = processor_particle_count(target_loc, P, point_count);
        particle_thread[i] = target_loc;
        // if (target_loc == ID) {
        //     in_curr_thread[i] = true;
        //     own_particles[local_id] = i;
        //     local_ids[i] = local_id;
            // cout << "particle i " << i << " thread " << ID << " local id " << local_id << endl;
        // }
    }

    for (int i = 0; i < P; i++) {
        if (i < P - 1) {
            lower_bounds[i] = points_per_rank * i;
            upper_bounds[i] = points_per_rank * (i + 1);
            // point_count[i] = points_per_rank;
        } else {
            lower_bounds[i] = points_per_rank * i;
            upper_bounds[i] = point_count;
            // point_count[i] =
        }
        point_counts[i] = upper_bounds[i] - lower_bounds[i];
    }


    if (ID == 0) { // reads in the particles
        vector<double> particle (5, 0);
        for (int i = 0; i < point_count; i++) {
            getline(file1, line);
            stringstream str(line);
            target_loc = particle_thread[i];
            local_id = i - target_loc * points_per_rank;
            target_points = processor_particle_count(target_loc, P, point_count);
            for (int j = 0; j < 5; j++) {
                getline(str, word, ',');
                particle[j] = stod(word);
            }
            if (target_loc == 0) {
                for (int j = 0; j < 5; j++) curr_state[5 * local_id + j] = particle[j];
            } else {
                // cout << "here " << i << endl;
                MPI_Put(&particle[0], 5, MPI_DOUBLE, target_loc, 5 * local_id, 5, MPI_DOUBLE, win_curr_state);
            }
        }
        begin_time = chrono::steady_clock::now();

        // for (int i = 0; i < point_count; i++) { // write out the initial state
        //     cout << "particle " << i << endl;
        //     if (in_curr_thread[i]) {
        //         local_id = local_ids[i];
        //         write_out << curr_state[4 * local_id] << "," << curr_state[4 * local_id + 1] << "," << curr_state[4 * local_id + 2] << "," << curr_state[4 * local_id + 3] << "\n";
        //     } else {
        //         particle = {0, 0, 0, 0};
        //         target_loc = particle_thread[i];
        //         local_id = i - target_loc * points_per_rank;
        //         cout << "target " << target_loc << " local id " << local_id << endl;
        //         // MPI_Get(&particle[0], 4, MPI_DOUBLE, target_loc, 4 * local_id, 4, MPI_DOUBLE, win_curr_state);
        //         write_out << particle[0] << "," << particle[1] << "," << particle[2] << "," << particle[3] << "\n";
        //     }
        // }
    }
    // MPI_Barrier(MPI_COMM_WORLD);
    //
    // for (int i = 0; i < point_count; i++) { // write out the initial state
    //     if (in_curr_thread[i]) {
    //         local_id = i - ID * points_per_rank;
    //         int offset = sizeof(double) * i * 4;
    //         // cout << "offset: " << offset << endl;
    //         // MPI_File_write_at_all(fh, offset, &curr_state[4*local_id], 4, MPI_DOUBLE, &status);
    //     }
    // }
    MPI_Barrier(MPI_COMM_WORLD);
    // MPI_Win_fence(0, win_curr_state);
    // vector<double> particle (3, 0);
    // MPI_Get(&particle[0], 3, MPI_DOUBLE, 1, 0, 3, MPI_DOUBLE, win_curr_state);
    // MPI_Get(&particle[0], 3, MPI_DOUBLE, 0, 0, 3, MPI_DOUBLE, win_curr_state);
    // MPI_Win_fence(0, win_curr_state);
    // cout << ID << " particle " << particle[0]  << " " << particle[1] << " " << particle[2] << endl;


    //
    cout << "here" << endl;
    for (int t = 0; t < 1; t++) { // time iterate with RK4
        double curr_time = t * delta_t;
        // MPI_Win_fence(0, win_curr_state);
        sync_buffer(curr_state, ID, P, lower_bounds, point_counts, &win_curr_state);
        MPI_Barrier(MPI_COMM_WORLD);
        BVE_ffunc(c_1, curr_state, curr_time, delta_t, omega, area, point_count, lower_bounds[ID], upper_bounds[ID]);
        // intermediate_1 = c_1;
        // scalar_mult(intermediate_1, delta_t / 2);
        // vec_add(intermediate_1, curr_state);
        // MPI_Barrier(MPI_COMM_WORLD);
        // sync_buffer(intermediate_1, ID, P, lower_bounds, point_counts, &win_inter1);
        // MPI_Barrier(MPI_COMM_WORLD);
        // BVE_ffunc(c_2, intermediate_1, curr_time, delta_t, omega, area, point_count, lower_bounds[ID], upper_bounds[ID]);
        // intermediate_2 = c_2;
        // scalar_mult(intermediate_2, delta_t / 2);
        // vec_add(intermediate_2, curr_state);
        // MPI_Barrier(MPI_COMM_WORLD);
        // sync_buffer(intermediate_2, ID, P, lower_bounds, point_counts, &win_inter2);
        // MPI_Barrier(MPI_COMM_WORLD);
        // BVE_ffunc(c_3, intermediate_2, curr_time, delta_t, omega, area, point_count, lower_bounds[ID], upper_bounds[ID]);
        // intermediate_3 = c_3;
        // scalar_mult(intermediate_3, delta_t);
        // vec_add(intermediate_3, curr_state);
        // MPI_Barrier(MPI_COMM_WORLD);
        // sync_buffer(intermediate_3, ID, P, lower_bounds, point_counts, &win_inter3);
        // MPI_Barrier(MPI_COMM_WORLD);
        // BVE_ffunc(c_4, intermediate_3, curr_time, delta_t, omega, area, point_count, lower_bounds[ID], upper_bounds[ID]);
        // c1234 = c_1;
        // scalar_mult(c_2, 2);
        // vec_add(c1234, c_2);
        // scalar_mult(c_3, 2);
        // vec_add(c1234, c_3);
        // vec_add(c1234, c_4);
        // scalar_mult(c1234, delta_t / 6);
        // vec_add(curr_state, c1234);
        // for (int i = lower_bounds[ID]; i < upper_bounds[ID]; i++) {
        //     vector<double> projected = slice(curr_state, 4 * i, 1, 3);
        //     // projected = project_to_sphere(projected, 1);
        //     project_to_sphere(projected, 1);
        //     for (int j = 0; j < 3; j++) curr_state[4 * i + j] = projected[j]; // reproject points to surface of sphere
        //     // write_out << curr_state[4 * i] << "," << curr_state[4 * i + 1] << "," << curr_state[4 * i + 2] << "," << curr_state[4 * i + 3] << "\n"; // write position
        //     // int offset = ((t + 1) * point_count + own_particles[i]) * 4 * sizeof(double);
        //     // cout << "offset: " << offset << endl;
        //     // MPI_File_write_at_all(fh, offset, &curr_state[4*i], 4, MPI_DOUBLE, &status);
        // }
        MPI_Barrier(MPI_COMM_WORLD);
        // for (int i = 0; i < point_count; i++) { // write out the initial state
        //     if (in_curr_thread[i]) {
        //         // local_id = i - ID * points_per_rank;
        //         int local_id = local_ids[i];
        //         int offset = ((t + 1) * point_count + i) * 4 * sizeof(double);
        //         // cout << "offset: " << offset << endl;
        //         MPI_File_write_at_all(fh, offset, &curr_state[4*local_id], 4, MPI_DOUBLE, &status);
        //     }
        // }
        // MPI_Barrier(MPI_COMM_WORLD);

        if (ID == 0) cout << t << endl;
    }

    if (ID == 0) {
        chrono::steady_clock::time_point end_time = chrono::steady_clock::now();
        cout << "time taken: " << chrono::duration_cast<chrono::microseconds>(end_time - begin_time).count() << " microseconds" << endl;
    }

    MPI_Barrier(MPI_COMM_WORLD);

    // MPI_File_close(&fh);

    MPI_Win_fence(0, win_curr_state);
    // cout << "here" << endl;
    MPI_Finalize();
    return 0;
}
