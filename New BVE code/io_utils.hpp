#ifndef io_H
#define io_H

#include <fstream>
#include <sstream>
#include "structs.hpp"

void write_state(run_config& run_information, vector<double>& dynamics_state, vector<double>& dynamics_area, ofstream& file_writer1, ofstream& file_writer2) {
    for (int i = 0; i < run_information.dynamics_curr_point_count; i++) { // write out initial state
        for (int j = 0; j < run_information.info_per_point; j++) {
            file_writer1 << setprecision(15) << dynamics_state[run_information.info_per_point * i + j] << ",";
        }
        file_writer1 << setprecision(15) << dynamics_area[i] << "\n";
    }
    file_writer2 << run_information.dynamics_curr_point_count << "\n";
}

void write_triangles(run_config& run_information, vector<vector<vector<int>>>& dynamics_triangles,
            vector<vector<bool>>& dynamics_triangles_is_leaf, ofstream& file_writer3, ofstream& file_writer4) {
    // for (int i = 0; i < 20 * pow(4, run_information.dynamics_levels_min - 1); i++) {
    //     for (int j = 0; j < 3; j++) {
    //         file_writer3 << dynamics_triangles[run_information.dynamics_levels_min-1][i][j] << ",";
    //     }
    //     file_writer3 << "\n";
    // }
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

#endif
