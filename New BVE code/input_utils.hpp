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
        if (word1 == "use_mpi") {
            if (stoi(word2) == 1) {
                run_information.use_mpi = true;
            }
        } else if (word1 == "use_amr") {
            if (stoi(word2) == 1) {
                run_information.use_amr = true;
            }
        } else if (word1 == "use_remesh") {
            if (stoi(word2) == 1) {
                run_information.use_remesh = true;
            }
        } else if (word1 == "use_fast") {
            if (stoi(word2) == 1) {
                run_information.use_fast = true;
            }
        } else if (word1 == "radius") {
            run_information.radius = stod(word2);
        } else if (word1 == "end_time") {
            run_information.end_time = stod(word2);
        } else if (word1 == "delta_t") {
            run_information.delta_t = stod(word2);
        } else if (word1 == "dynamics_levels_min") {
            run_information.dynamics_levels_min = stoi(word2);
        } else if (word1 == "dynamics_levels_max") {
            run_information.dynamics_levels_max = stoi(word2);
        } else if (word1 == "initial_vor_condition") {
            run_information.initial_vor_condition = word2;
        } else if (word1 == "fast_sum_cluster_thresh") {
            run_information.fast_sum_cluster_thresh = stoi(word2);
        } else if (word1 == "fast_sum_tree_levels") {
            run_information.fast_sum_tree_levels = stoi(word2);
        } else if (word1 == "fast_sum_theta") {
            run_information.fast_sum_theta = stod(word2);
        } else if (word1 == "interp_degree") {
            run_information.interp_degree = stoi(word2);
        } else if (word1 == "info_per_point") {
            run_information.info_per_point = stoi(word2);
        } else {
            run_information.time_steps = int(run_information.end_time / run_information.delta_t);
            run_information.dynamics_initial_points = 10 * pow(4, run_information.dynamics_levels_min) + 2;
            run_information.dynamics_initial_triangles = 20 * pow(4, run_information.dynamics_levels_min);
            run_information.dynamics_max_points = 10 * pow(4, run_information.dynamics_levels_max) + 2;
            run_information.dynamics_max_triangles = 20 * pow(4, run_information.dynamics_levels_max);
            return;
        }
    }

    // getline(config_file, line);
    // // cout << line << endl;
    // stringstream str1(line);
    //
    // getline(str1, word, '=');
    // if (word == "mpi") {
    //     getline(str1, word);
    //     if (stoi(word) == 1) {
    //         run_information.use_mpi = true;
    //     }
    // }
    // getline(config_file, line);
    // // stringstream str2(line);
    // cout << line << endl;
    // if (line == "") {
    //     cout << "here" << endl;
    //     return;
    // }
    // cout << "here2" << endl;
    // switch (word) {
    //     case "mpi":
    //         cout << "here" << endl;
    //         break;
    // }
    // cout << word << endl;
    // getline(str1, word);
    // cout << word << endl;
}


#endif
