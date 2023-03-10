#include <iostream>
#include <cmath>
#include <array>
#include <vector>
#include <Accelerate/Accelerate.h>
#include <fstream>
#include <queue>
#include <sstream>
#include <chrono>
#include <cassert>
#include <mpi.h>

#include "helpers.h"
#include "vorticity_functions.h"
#include "rhs_eval.h"
#include "icos_utils.h"
#include "amr.h"
#include "mpi_utils.h"
#include "remesh.h"
#include "fast_rhs.h"
#include "interp_utils.h"

// compile with mpic++ -O3 -std=c++11 -framework Accelerate new_dynamics.cpp -o new_dynamics

int main(int argc, char** argv) { // run the dynamics
    // argv[1] is whether to use direct or fast summation, d for direct, f for fast
    // argv[2] is whether to use amr or not, n for no, a for amr
    // argv[3] is whether ot not to use mpi, n for no, m for mpi
    // argv[4] is the initial vorticity distribution
    // argv[5] is whether to remesh or not, n for no, r for remesh
    // argv[6] is the number of points, options 642, 2562, 10242, 40962, 163842

    double delta_t = 0.01; // time step in days
    double end_t = 3.0; //end time in days
    // int point_count = 2562; // number of points
    int point_count = stoi(argv[6]);
    int old_point_count; // previous number of points if AMR is being used
    int tri_count = 2 * (point_count - 2); // number of triangle faces
    int time_steps = end_t / delta_t; // number of time steps, for real runs
    // int time_steps = 1; // testing configuration
    double omega = 2 * M_PI; // coriolis factor
    double radius = 1.0; // radius of the sphere
    bool use_mpi = false, use_amr = false, use_remesh = false;
    int mpi_ranks;

    vector<double> curr_state(5 * point_count); // 0 is x_pos, 1 is y_pos, 2 is z_pos, 3 is vorticity, 4 is passive tracer
    vector<double> c_1(5 * point_count, 0);
    vector<double> c_2(5 * point_count, 0);
    vector<double> c_3(5 * point_count, 0);
    vector<double> c_4(5 * point_count, 0);
    vector<double> c1234(5 * point_count, 0);
    vector<double> intermediate_1(5 * point_count);
    vector<double> intermediate_2(5 * point_count);
    vector<double> intermediate_3(5 * point_count);
    vector<double> area(point_count, 0);

    vector<vector<int>> triangles(tri_count, vector<int> (3)); // triangles[i] is a length 3 int vector containing the indices of the points of the vertices
    vector<vector<int>> vert_tris(point_count); // vert_tris[i] is the int vector containing the indices of the triangles adjacent to point i

    // read the points, replace eventually with inline point generation
    read_points(point_count, tri_count, curr_state, vert_tris, triangles);

    // initialize areas
    area_init(curr_state, area, triangles, tri_count);

    // initialize vorticity and tracer
    vorticity_init(curr_state, area, point_count, argv[4]);

    if (argv[1][0] == 'f') { // fast summation init
        int icos_levels = 5; // how much to refine the matrix
        double theta = 0.7; // cluster acceptance threshold
        int many_count = 10; // cluster acceptance count
        int degree = 2; // degree of interpolation for clusters
        int cluster_count = (degree + 1) * (degree + 2) / 2; // number of interpolation points for clusters
        vector<vector<int>> point_locs (icos_levels, vector<int> (point_count, 0)); // triangle each point is in
        icos_struct icos1 = icosahedron_init(icos_levels, radius);
        interp_struct interp1 = interp_init(degree);
        points_assign(icos1, curr_state, point_locs, point_count);

        vector<interaction_pair> interaction_list;
    }


    if (argv[2][0] == 'a') { // initialize for amr
        int max_points = 10 * point_count; // limit on amr
        vector<vector<int>> parent_verts(max_points, vector<int> (2, 0)); // the two points that a new point comes from
        use_amr = true;
    }

    if (argv[3][0] == 'm') { // use MPI
        use_mpi = true;
        MPI_Init(&argc, &argv);
        int P; // number of MPI_ranks
        int ID; // process ID
        MPI_Status status;
        MPI_Comm_size(MPI_COMM_WORLD, &P);
        MPI_Comm_rank(MPI_COMM_WORLD, &ID);
        MPI_Win win_curr_state;
        MPI_Win win_inter1;
        MPI_Win win_inter2;
        MPI_Win win_inter3;
        MPI_Win win_c1234;
    }

    if (argv[5][0] == 'r') {
        use_remesh = true;
    }

    chrono::steady_clock::time_point begin = chrono::steady_clock::now();

    string run_config = to_string(point_count) + "_" + argv[4] + "_" + argv[1]; // point count, configuration, and direct/fast
    if (use_amr) {
        run_config += "_";
        run_config += argv[2];
    }
    if (use_mpi) {
        run_config += "_";
        run_config += argv[3];
    }
    if (use_remesh) {
        run_config += "_";
        run_config += argv[5];
    }
    run_config += "_" + to_string(time_steps);
    ofstream write_out1("output_" + run_config + ".csv", ofstream::out | ofstream::trunc); // ofstream = output file stream
    ofstream write_out2("point_counts_" + run_config + ".csv", ofstream::out | ofstream::trunc); // at each time step, write out the number of points

    for (int i = 0; i < point_count; i++) { // write out initial state
        write_out1 << curr_state[5 * i] << "," << curr_state[5 * i + 1] << "," << curr_state[5 * i + 2] << "," << curr_state[5 * i + 3] << "," << curr_state[5 * i + 4] <<  "," << area[i] << "\n";
    }
    write_out2 << point_count << "\n";

    for (int t = 0; t < time_steps; t++) {
        double curr_time = t * delta_t;
        // vector<double> c_1(5 * point_count, 0);
        // vector<double> c_2(5 * point_count, 0);
        // vector<double> c_3(5 * point_count, 0);
        // vector<double> c_4(5 * point_count, 0);
        // vector<double> c1234(5 * point_count, 0);
        // vector<double> intermediate_1(5 * point_count);
        // vector<double> intermediate_2(5 * point_count);
        // vector<double> intermediate_3(5 * point_count);
        // cout << " state curr " << curr_state[5 * 1281] << endl;
        // BVE_ffunc(c_1, curr_state, curr_time, delta_t, omega, area, point_count);
        rhs_func(c_1, curr_state, area, omega, point_count, argv[1][0]);
        intermediate_1 = c_1;
        scalar_mult(intermediate_1, delta_t / 2);
        vec_add(intermediate_1, curr_state);
        // cout << " state inter1 " << intermediate_1[5 * 1281] << endl;
        // BVE_ffunc(c_2, intermediate_1, curr_time + delta_t / 2, delta_t, omega, area, point_count);
        rhs_func(c_2, intermediate_1, area, omega, point_count, argv[1][0]);
        intermediate_2 = c_2;
        scalar_mult(intermediate_2, delta_t / 2);
        vec_add(intermediate_2, curr_state);
        // BVE_ffunc(c_3, intermediate_2, curr_time + delta_t / 2, delta_t, omega, area, point_count);
        rhs_func(c_3, intermediate_2, area, omega, point_count, argv[1][0]);
        intermediate_3 = c_3;
        scalar_mult(intermediate_3, delta_t);
        vec_add(intermediate_3, curr_state);
        // BVE_ffunc(c_4, intermediate_3, curr_time + delta_t, delta_t, omega, area, point_count);
        rhs_func(c_4, intermediate_3, area, omega, point_count, argv[1][0]);
        c1234 = c_1;
        scalar_mult(c_2, 2);
        vec_add(c1234, c_2);
        scalar_mult(c_3, 2);
        vec_add(c1234, c_3);
        vec_add(c1234, c_4);
        scalar_mult(c1234, delta_t / 6);
        if (use_remesh) {
            vec_add(c1234, curr_state);
            regrid_points2(c1234, curr_state, triangles, vert_tris, point_count, tri_count, omega, 0, point_count, 0); // regrids points so that they are regular, modifies curr_state
        } else {
            vec_add(curr_state, c1234);
        }
        // state = amr(curr_state, triangles, vert_tris, area, parent_verts, tri_count, point_count, max_points);
        // point_count = state[1];
        // tri_count = state[0];
        for (int i = 0; i < point_count; i++) {
            vector<double> projected = slice(curr_state, 5 * i, 1, 3);
            project_to_sphere(projected, 1);
            for (int j = 0; j < 3; j++) curr_state[5 * i + j] = projected[j]; // reproject points to surface of sphere
            write_out1 << curr_state[5 * i] << "," << curr_state[5 * i + 1] << "," << curr_state[5 * i + 2] << "," << curr_state[5 * i + 3] << "," << curr_state[5 * i + 4] << "," << area[i] << "\n"; // write current state
        }
        write_out2 << point_count << "\n";
        cout << "t: " << t << " point_count " << point_count << endl;
    }

    chrono::steady_clock::time_point end = chrono::steady_clock::now();
    cout << "time taken: " << chrono::duration_cast<chrono::microseconds>(end - begin).count() << " microseconds" << endl;

    write_out1.close();
    write_out2.close();

    return 0;
}
