#include <cmath>
#include <vector>
#include <Accelerate/Accelerate.h> // lapack on mac
#include <fstream>
#include <queue>
#include <sstream>
#include <chrono>
#include <cassert>
#include <mpi.h>

#include "general_utils.hpp"
// #include "interp_utils.hpp"
#include "init_utils.hpp"
#include "structs.hpp"
#include "input_utils.hpp"
#include "io_utils.hpp"
#include "rhs_utils.hpp"

using namespace std;

double omega = 2 * M_PI;

int main() {

    run_config run_information;
    read_run_config("namelist.txt", run_information); // reads in run configuration information

    vector<double> dynamics_state; // list of points and other information in a flattened array
    vector<vector<vector<int>>> dynamics_triangles; // at level i, each entry is a vector which contains the 3 vertices and the refinement level of the triangle
    vector<vector<vector<int>>> dynamics_points_adj_triangles; // for point i, level j, the triangles it touches at that level
    vector<vector<int>> dynamics_parent_triangles; // for triangle i, the triangle one level above that it comes from
    vector<vector<int>> dynamics_child_triangles; // for triangle i, the four triangles coming from it
    vector<vector<int>> dynamics_points_parents; // for point i, the two points that it is the midpoint of

    // dynamics_points_initialize(run_information, dynamics_state, dynamics_triangles, dynamics_points_adj_triangles, dynamics_parent_triangles, dynamics_child_triangles, dynamics_points_parents);
    dynamics_points_initialize(run_information, dynamics_state, dynamics_triangles, dynamics_points_parents);

    vector<double> dynamics_areas (run_information.dynamics_initial_points, 0);
    area_initialize(run_information, dynamics_state, dynamics_triangles, dynamics_areas); // finds areas for each point
    vorticity_initialize(run_information, dynamics_state, dynamics_areas); // initializes vorticity values for each point
    tracer_initialize(run_information, dynamics_state); // initializes tracer values for each point


    string output_filename = to_string(run_information.dynamics_initial_points) + "_" + run_information.initial_vor_condition + "_";
    if (run_information.use_fast) {
        output_filename += "fast_" + to_string(run_information.fast_sum_tree_levels) + "_" + to_string(run_information.fast_sum_theta);
    } else output_filename += "direct";
    if (run_information.use_amr) output_filename += "_amr";
    if (run_information.use_mpi) output_filename += "_mpi";
    if (run_information.use_remesh) output_filename += "_remesh";
    stringstream ss;
    int precision;
    if (run_information.end_time >= 1.0) precision = 0;
    else precision = -log10(run_information.end_time);
    ss << fixed << setprecision(precision) << run_information.end_time;
    output_filename += "_" + ss.str();

    ofstream write_out1("./run-output/output_" + output_filename + ".csv", ofstream::out | ofstream::trunc); // ofstream = output file stream
    ofstream write_out2("./run-output/point_counts_" + output_filename + ".csv", ofstream::out | ofstream::trunc); // at each time step, write out the number of points

    write_state(run_information, dynamics_state, dynamics_areas, write_out1, write_out2);

    chrono::steady_clock::time_point begin = chrono::steady_clock::now();

    for (int t = 0; t < run_information.time_steps; t++ ) {
        vector<double> c1234(run_information.info_per_point * run_information.dynamics_curr_point_count, 0);
        vector<double> inter_state(run_information.info_per_point * run_information.dynamics_curr_point_count, 0);
        vector<double> rhs_update(run_information.info_per_point * run_information.dynamics_curr_point_count, 0);
        rhs_func(run_information, rhs_update, dynamics_state, dynamics_areas, omega); // RK4 k_1
        c1234 = rhs_update; // k_1
        inter_state = dynamics_state; // k_1
        scalar_mult(inter_state, run_information.delta_t / 2.0); // delta_t/2*k1
        vec_add(inter_state, rhs_update); // x+delta_t/2*k1
        rhs_func(run_information, rhs_update, inter_state, dynamics_areas, omega); // RK4 k_2
        inter_state = rhs_update; // k_2
        scalar_mult(rhs_update, 2.0); // 2k_2
        vec_add(c1234, rhs_update); // k_1 + 2k_2
        scalar_mult(inter_state, run_information.delta_t / 2.0); // delta_t/2 * k_2
        vec_add(inter_state, dynamics_state); // x+delta_t/2*k2
        rhs_func(run_information, rhs_update, inter_state, dynamics_areas, omega); // RK4 k_3
        inter_state = rhs_update; // k_3
        scalar_mult(rhs_update, 2.0); // 2k_3
        vec_add(c1234, rhs_update); // k_1 + 2k_2 + 2k_3
        scalar_mult(inter_state, run_information.delta_t); // delta_t * k_3
        vec_add(inter_state, dynamics_state); // x + delta_t * k_3
        rhs_func(run_information, rhs_update, inter_state, dynamics_areas, omega); // RK4 k_4
        vec_add(c1234, rhs_update); // k_1 + 2k_2 + 2k_3 + k_4
        scalar_mult(c1234, run_information.delta_t / 6.0); // RK4 update
        vec_add(dynamics_state, c1234);

        write_state(run_information, dynamics_state, dynamics_areas, write_out1, write_out2);

        // progress the dynamics
    }

    chrono::steady_clock::time_point end = chrono::steady_clock::now();
    cout << "time taken: " << chrono::duration_cast<chrono::microseconds>(end - begin).count() << " microseconds" << endl;
    // cout << dynamics_state[0] << " " << dynamics_state[1] << " " <<  dynamics_state[2] << " " << dynamics_state[3] << endl;
    // cout << run_information.dynamics_curr_point_count << endl;
    // cout << run_information.dynamics_initial_points << endl;
    return 0;
}
