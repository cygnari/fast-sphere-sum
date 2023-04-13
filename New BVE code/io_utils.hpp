#ifndef io_H
#define io_H

#include <fstream>
#include <sstream>
#include "structs.hpp"

void write_state(run_config& run_information, vector<double>& dynamics_state, vector<double>& dynamics_area, ofstream& file_writer1, ofstream& file_writer2) {
    for (int i = 0; i < run_information.dynamics_curr_point_count; i++) { // write out initial state
        for (int j = 0; j < run_information.info_per_point; j++) {
            file_writer1 << dynamics_state[run_information.info_per_point * i + j] << ",";
        }
        file_writer1 << dynamics_area[i] << "\n";
    }
    file_writer2 << run_information.dynamics_curr_point_count << "\n";
}

#endif
