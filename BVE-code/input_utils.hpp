#ifndef input_H
#define input_H

#include <fstream>
#include <sstream>
#include "structs.hpp"
#include <iostream>

using namespace std;

void read_run_config(string file_name, run_config& run_information) {
    // reads run information of file_name
    ifstream config_file(file_name);
    string line, word1, word2;

    while (true) {
        getline(config_file, line);
        stringstream str1(line);
        getline(str1, word1, '=');
        getline(str1, word2);
        if (word1 == "use_amr") {
            if (stoi(word2) == 1) {
                run_information.use_amr = true;
                run_information.use_remesh = true; // amr cannot happen without remeshing
            }
        } else if (word1 == "use_remesh") {
            if (stoi(word2) == 1) {
                run_information.use_remesh = true;
            }
        } else if (word1 == "use_fast") {
            if (stoi(word2) == 1) {
                run_information.use_fast = true;
            }
        } else if (word1 == "force_conservative") {
            if (stoi(word2) == 1) {
                run_information.use_fixer = true;
            }
        } else if (word1 == "vor_fix") {
            if (stoi(word2) >= 1) {
                run_information.vor_fix = true;
            }
            if (stoi(word2) >= 2) {
                run_information.vor_limiter = true;
            }
        } else if (word1 == "out_path") {
            run_information.out_path = word2;
        } else if (word1 == "write_output") {
            if (stoi(word2) == 1) {
                run_information.write_output = true;
            }
        } else if (word1 == "write_precision") {
            run_information.write_precision = stoi(word2);
        } else if (word1 == "write_tris") {
            if (stoi(word2) == 1) {
                run_information.write_tris = true;
            }
        } else if (word1 == "write_stream") {
            if (stoi(word2) == 1) {
                run_information.write_stream = true;
            }
        } else if (word1 == "radius") {
            run_information.radius = stod(word2);
        } else if (word1 == "end_time") {
            run_information.end_time = stod(word2);
        } else if (word1 == "delta_t") {
            run_information.delta_t = stod(word2);
        } else if (word1 == "dynamics_levels_min") {
            run_information.dynamics_levels_min = stoi(word2);
        } else if (word1 == "amr_levels_max") {
            if (run_information.use_amr) {
                run_information.amr_levels = stoi(word2);
                run_information.dynamics_levels_max = run_information.amr_levels + run_information.dynamics_levels_min;
            } else {
                run_information.dynamics_levels_max = run_information.dynamics_levels_min;
            }
        } else if (word1 == "initial_vor_condition") {
            run_information.initial_vor_condition = word2;
        } else if (word1 == "vor_force") {
            run_information.vor_forcing = word2;
        } else if (word1 == "icp1") {
            run_information.init_cond_param1 = stoi(word2);
        } else if (word1 == "icp2") {
            run_information.init_cond_param2 = stod(word2);
        } else if (word1 == "frp1") {
            run_information.forcing_param1 = stoi(word2);
        } else if (word1 == "frp2") {
            run_information.forcing_param2 = stod(word2);
        } else if (word1 == "fast_sum_cluster_thresh") {
            run_information.fast_sum_cluster_thresh = stoi(word2);
        } else if (word1 == "fast_sum_tree_levels") {
            run_information.fast_sum_tree_levels = stoi(word2);
            if (run_information.fast_sum_tree_levels == -1) {
                run_information.fast_sum_tree_levels = run_information.dynamics_levels_min - 3;
            }
        } else if (word1 == "fast_sum_theta") {
            run_information.fast_sum_theta = stod(word2);
        } else if (word1 == "fast_sum_rotate") {
            if (stoi(word2) == 0) {
                run_information.fast_sum_rotate = false;
            }
        } else if (word1 == "interp_degree") {
            run_information.interp_degree = stoi(word2);
        } else if (word1 == "tracer_count") {
            run_information.tracer_count = stoi(word2);
            run_information.info_per_point = run_information.tracer_count + 4;
        } else if (word1 == "amr_vor_thresh") {
            run_information.amr_vor_thresh = stod(word2);
        } else if (word1 == "amr_circ_thresh") {
            run_information.amr_circ_thresh = stod(word2);
        } else if (word1 == "fast_rot_alph") {
            run_information.fast_sum_rotate_alph = stod(word2);
        } else if (word1 == "fast_rot_beta") {
            run_information.fast_sum_rotate_beta = stod(word2);
        } else if (word1 == "fast_rot_gamm") {
            run_information.fast_sum_rotate_gamm = stod(word2);
        } else {
            run_information.time_steps = int(run_information.end_time / run_information.delta_t);
            run_information.dynamics_initial_points = 10 * pow(4, run_information.dynamics_levels_min - 1) + 2;
            run_information.dynamics_initial_triangles = 20 * pow(4, run_information.dynamics_levels_min - 1);
            run_information.dynamics_max_points = 10 * pow(4, run_information.dynamics_levels_max - 1) + 2;
            run_information.dynamics_max_triangles = 20 * pow(4, run_information.dynamics_levels_max - 1);
            run_information.dynamics_curr_point_count = run_information.dynamics_initial_points;
            run_information.dynamics_curr_tri_count = run_information.dynamics_initial_triangles;
            run_information.interp_point_count = int((1 + run_information.interp_degree) * (2 + run_information.interp_degree) / 2);
            if (run_information.write_stream) run_information.info_per_point += 1;
            return;
        }
    }
}


#endif
