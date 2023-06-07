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
#include "struct_list.h"
#include "vorticity_functions.h"
#include "rhs_eval.h"
#include "icos_utils.h"
#include "amr.h"
#include "mpi_utils.h"
#include "remesh.h"
#include "fast_rhs.h"
#include "interp_utils.h"
#include "run_config.h"

// compile with mpic++ -O3 -std=c++11 -framework Accelerate new_dynamics.cpp -o new_dynamics

int main(int argc, char** argv) { // run the dynamics

    run_information config1;
    config_init(argc, argv, config1, 1.0);

    double delta_t = 0.01; // time step in days
    // double end_t = 3.0; //end time in days
    // double end_t = 1.0;
    // int old_point_count; // previous number of points if AMR is being used
    // int tri_count = 2 * (config1.point_count - 2); // number of triangle faces
    // config1.tri_count = 2 * (config1.point_count - 2);
    int time_steps;
    if (config1.testing) time_steps = 1; // testing configuration
    else time_steps = config1.end_time / delta_t;
    // int time_steps = config1.end_time / delta_t; // number of time steps, for real runs
    // int time_steps = 1; // testing configuration
    double omega = 2 * M_PI; // coriolis factor
    // double radius = 1.0; // radius of the sphere

    vector<double> curr_state(5 * config1.point_count); // 0 is x_pos, 1 is y_pos, 2 is z_pos, 3 is vorticity, 4 is passive tracer
    config1.dynamics_state = curr_state;
    // vector<double> c_1(5 * config1.point_count, 0);
    // vector<double> c_2(5 * config1.point_count, 0);
    // vector<double> c_3(5 * config1.point_count, 0);
    // vector<double> c_4(5 * config1.point_count, 0);
    vector<double> inter_state(5 * config1.point_count); // use this for the intermediate RK4 steps
    vector<double> c1234(5 * config1.point_count, 0); // accumulate the RK4 update here
    // vector<double> intermediate_1(5 * config1.point_count);
    // vector<double> intermediate_2(5 * config1.point_count);
    // vector<double> intermediate_3(5 * config1.point_count);
    vector<double> area(config1.point_count, 0);
    config1.dynamics_areas = area;

    vector<vector<int>> triangles(config1.tri_count, vector<int> (3)); // triangles[i] is a length 3 int vector containing the indices of the points of the vertices
    vector<vector<int>> vert_tris(config1.point_count); // vert_tris[i] is the int vector containing the indices of the triangles adjacent to point i
    config1.amr_vert_tris = vert_tris;
    config1.dynamics_tris = triangles;

    // read the points, replace eventually with inline point generation
    read_points(config1, curr_state);

    // initialize areas
    area_init(curr_state, area, config1.dynamics_tris, config1.tri_count);

    // initialize vorticity and tracer
    vorticity_init(curr_state, area, config1.point_count, config1.init_cond);



    if (config1.use_fast) { // fast summation init
        // int icos_levels = 3; // how much to refine the matrix
        config1.many_count = 10; // cluster acceptance count
        config1.degree = 2; // degree of interpolation for clusters
        vector<vector<int>> point_locs (config1.tree_levels, vector<int> (config1.point_count, 0)); // triangle each point is in
        // icos_struct icos1;
        // interp_struct interp1;
        config1.tree_point_locs = point_locs;
        icosahedron_init(config1, 1.0);
        // interp1 = interp_init(degree);
        interp_init(config1);
        points_assign(config1, curr_state);
        tree_traverse(config1);

        // config1.icos = icos1;
        // config1.interp = interp1;
        // config1.many_count = many_count;
        // config1.theta = theta;
    }

    if (config1.use_amr) { // initialize for amr
        // amr_struct amr1;
        // config1.max_count = 10 * config1.point_count;
        vector<vector<int>> parent_verts(config1.max_count, vector<int> (2, 0)); // the two points that a new point comes from
        // amr1.parent_verts = parent_verts;
        // config1.amr = amr1;
        // config1.max_count = 10 * config1.point_count;
        config1.amr_parent_verts = parent_verts;
    }

    if (config1.use_mpi) { // use MPI
        // use_mpi = true;
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
        config1.mpi_ranks = P;
    }

    chrono::steady_clock::time_point begin = chrono::steady_clock::now();

    string out_config = to_string(config1.point_count) + "_" + config1.init_cond + "_"; // point count, configuration
    if (config1.use_fast) {
        out_config += "fast_" + to_string(config1.tree_levels);
    }
    else out_config += "direct";
    if (config1.use_amr) out_config += "_amr";
    if (config1.use_mpi) out_config += "_mpi";
    if (config1.use_remesh) out_config += "_remesh";

    out_config += "_" + to_string(time_steps);
    ofstream write_out1("output_" + out_config + ".csv", ofstream::out | ofstream::trunc); // ofstream = output file stream
    ofstream write_out2("point_counts_" + out_config + ".csv", ofstream::out | ofstream::trunc); // at each time step, write out the number of points

    for (int i = 0; i < config1.point_count; i++) { // write out initial state
        write_out1 << curr_state[5 * i] << "," << curr_state[5 * i + 1] << "," << curr_state[5 * i + 2] << "," << curr_state[5 * i + 3] << "," << curr_state[5 * i + 4] <<  "," << area[i] << "\n";
    }
    write_out2 << config1.point_count << "\n";

    for (int t = 0; t < time_steps; t++) {
        double curr_time = t * delta_t;
        // if (config1.use_amr) {
        //     cout << 5 * config1.point_count << endl;
        //     vector<double> c_1(5 * config1.point_count, 0);
        //     vector<double> c_2(5 * config1.point_count, 0);
        //     vector<double> c_3(5 * config1.point_count, 0);
        //     vector<double> c_4(5 * config1.point_count, 0);
        //     vector<double> c1234(5 * config1.point_count, 0);
        //     vector<double> intermediate_1(5 * config1.point_count);
        //     vector<double> intermediate_2(5 * config1.point_count);
        //     vector<double> intermediate_3(5 * config1.point_count);
        //     cout << intermediate_1.size() << endl;
        // }
        // cout << curr_state.size() << " " << intermediate_1.size() << " " << c_1.size() << endl;
        // vector<double> c_1(5 * config1.point_count, 0);
        // vector<double> c_2(5 * config1.point_count, 0);
        // vector<double> c_3(5 * config1.point_count, 0);
        // vector<double> c_4(5 * config1.point_count, 0);
        vector<double> c1234(5 * config1.point_count, 0);
        // vector<double> intermediate_1(5 * config1.point_count);
        // vector<double> intermediate_2(5 * config1.point_count);
        // vector<double> intermediate_3(5 * config1.point_count);
        vector<double> inter_state(5 * config1.point_count); // use this for the intermediate RK4 steps
        // cout << " state curr " << curr_state[5 * 1281] << endl;
        // BVE_ffunc(c_1, curr_state, curr_time, delta_t, omega, area, point_count);
        // rhs_func(c_1, curr_state, area, omega, point_count, argv[1][0]);
        // cout << "here 1" << endl;
        rhs_func(inter_state, curr_state, area, omega, config1);
        // intermediate_1 = c_1;
        c1234 = inter_state;
        // cout << "here 2" << endl;
        scalar_mult(intermediate_1, delta_t / 2);
        // cout << "here 2.5" << endl;
        // cout << curr_state.size() << " " << intermediate_1.size() << " " << c_1.size() << endl;
        vec_add(intermediate_1, curr_state);
        // cout << "here 3" << endl;
        // points_assign(config1.icos, intermediate_1, config1.icos.point_locs, config1.point_count);
        // tree_traverse(config1.icos.interactions, config1.icos, config1.theta, config1.many_count);
        // cout << "here" << endl;
        // cout << " state inter1 " << intermediate_1[5 * 1281] << endl;
        // BVE_ffunc(c_2, intermediate_1, curr_time + delta_t / 2, delta_t, omega, area, point_count);
        project_points(intermediate_1, config1.point_count, omega);
        rhs_func(c_2, intermediate_1, area, omega, config1);
        // cout << "here" << endl;
        intermediate_2 = c_2;
        scalar_mult(intermediate_2, delta_t / 2);
        vec_add(intermediate_2, curr_state);
        // cout << "here" << endl;
        // points_assign(config1.icos, intermediate_2, config1.icos.point_locs, config1.point_count);
        // tree_traverse(config1.icos.interactions, config1.icos, config1.theta, config1.many_count);
        // cout << "here" << endl;
        // BVE_ffunc(c_3, intermediate_2, curr_time + delta_t / 2, delta_t, omega, area, point_count);
        project_points(intermediate_2, config1.point_count, omega);
        rhs_func(c_3, intermediate_2, area, omega, config1);
        intermediate_3 = c_3;
        scalar_mult(intermediate_3, delta_t);
        vec_add(intermediate_3, curr_state);
        // points_assign(config1.icos, intermediate_3, config1.icos.point_locs, config1.point_count);
        // tree_traverse(config1.icos.interactions, config1.icos, config1.theta, config1.many_count);
        // BVE_ffunc(c_4, intermediate_3, curr_time + delta_t, delta_t, omega, area, point_count);
        project_points(intermediate_3, config1.point_count, omega);
        rhs_func(c_4, intermediate_3, area, omega, config1);
        c1234 = c_1;
        scalar_mult(c_2, 2);
        vec_add(c1234, c_2);
        scalar_mult(c_3, 2);
        vec_add(c1234, c_3);
        vec_add(c1234, c_4);
        scalar_mult(c1234, delta_t / 6);
        // project_points(curr_state, config1.point_count);
        if (config1.use_remesh) {
            vec_add(c1234, curr_state);
            project_points(c1234, config1.point_count, omega);
            regrid_points2(c1234, curr_state, config1.dynamics_tris, config1.amr_vert_tris, config1.point_count, config1.tri_count, omega, 0, config1.point_count, 0); // regrids points so that they are regular, modifies curr_state
        } else vec_add(curr_state, c1234);
        if (config1.use_amr) {
            // stuff
            // old_point_count = config1.point_count;
            // config1.point_count = amr(curr_state, triangles, vert_tris, area, config1.amr, config1.tri_count, config1.point_count);
            amr(curr_state, area, config1);
            if (config1.use_fast) {
                points_assign(config1, curr_state);
            }
        }
        // point_count = state[1];
        // tri_count = state[0];
        project_points(curr_state, config1.point_count, omega);
        for (int i = 0; i < config1.point_count; i++) {
            // vector<double> projected = slice(curr_state, 5 * i, 1, 3);
            // project_to_sphere(projected, 1);
            // for (int j = 0; j < 3; j++) curr_state[5 * i + j] = projected[j]; // reproject points to surface of sphere
            write_out1 << curr_state[5 * i] << "," << curr_state[5 * i + 1] << "," << curr_state[5 * i + 2] << "," << curr_state[5 * i + 3] << "," << curr_state[5 * i + 4] << "," << area[i] << "\n"; // write current state
        }
        write_out2 << config1.point_count << "\n";
        cout << "t: " << t << " point_count " << config1.point_count << endl;
    }

    chrono::steady_clock::time_point end = chrono::steady_clock::now();
    cout << "time taken: " << chrono::duration_cast<chrono::microseconds>(end - begin).count() << " microseconds" << endl;

    write_out1.close();
    write_out2.close();

    return 0;
}
