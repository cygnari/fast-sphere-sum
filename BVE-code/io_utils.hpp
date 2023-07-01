#ifndef io_H
#define io_H

#include <fstream>
#include <sstream>
#include "structs.hpp"
#include <filesystem>

void write_state(run_config& run_information, vector<double>& dynamics_state, vector<double>& dynamics_area, ofstream& file_writer1, ofstream& file_writer2) {
    for (int i = 0; i < run_information.dynamics_curr_point_count; i++) { // write out initial state
        for (int j = 0; j < run_information.info_per_point; j++) {
            file_writer1 << setprecision(run_information.write_precision) << dynamics_state[run_information.info_per_point * i + j] << ",";
        }
        file_writer1 << setprecision(run_information.write_precision) << dynamics_area[i] << "\n";
    }
    file_writer2 << run_information.dynamics_curr_point_count << "\n";
}

void write_triangles(run_config& run_information, vector<vector<vector<int>>>& dynamics_triangles,
            vector<vector<bool>>& dynamics_triangles_is_leaf, ofstream& file_writer3, ofstream& file_writer4) {
    for (int i = 0; i < run_information.dynamics_levels_max; i++) {
        for (int j = 0; j < 20 * pow(4, i); j++) {
            if (dynamics_triangles_is_leaf[i][j]) {
                for (int k = 0; k < 3; k++) {
                    file_writer3 << dynamics_triangles[i][j][k] << ",";
                }
                file_writer3 << "\n";
            }
        }
    }
    file_writer4 << run_information.dynamics_curr_tri_count << "\n";
}

string create_config(run_config& run_information) {
    stringstream ss1, ss2, ss3;
    int precision;
    string output_filename = to_string(run_information.dynamics_initial_points) + "_" + run_information.initial_vor_condition + "_";
    if (run_information.init_cond_param1 > 0) {
        output_filename += to_string(run_information.init_cond_param1) + "_";
    }
    if (run_information.init_cond_param2 > 0) {
        precision = max(int(ceil(-log10(run_information.init_cond_param2))), 0);
        ss1 << fixed << setprecision(precision) << run_information.init_cond_param2;
        output_filename += ss1.str() + "_";
    }
    if (run_information.vor_forcing != "none") {
        output_filename += run_information.vor_forcing + "_";
        if (run_information.forcing_param1 > 0) {
            output_filename += to_string(run_information.forcing_param1) + "_";
        }
        if (run_information.forcing_param2 > 0) {
            precision = max(int(ceil(-log10(run_information.forcing_param2))), 0);
            ss3 << fixed << setprecision(precision) << run_information.forcing_param2;
            output_filename += ss3.str() + "_";
        }
    }
    if (run_information.use_fast) {
        output_filename += "fast_" + to_string(run_information.fast_sum_tree_levels) + "_" + to_string(run_information.fast_sum_theta).substr(0, 3);
    } else output_filename += "direct";
    if (run_information.use_amr) output_filename += "_amr_" + to_string(run_information.amr_levels);
    if (run_information.use_remesh) output_filename += "_remesh";
    if (run_information.use_fixer) output_filename += "_fixer";

    precision = max(int(ceil(-log10(run_information.end_time))), 0);
    ss2 << fixed << setprecision(precision) << run_information.end_time;
    output_filename += "_" + ss2.str();
    return output_filename;
}

#endif
