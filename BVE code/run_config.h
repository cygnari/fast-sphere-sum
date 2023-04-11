#ifndef run_H
#define run_H

#include "struct_list.h"

void config_init(int argc, char** argv, run_information& config1, double radius) {
    int opt;
    config1.radius = radius;
    while ((opt = getopt(argc, argv, "marfi:c:t:l:xh:b:")) != -1) {
        // -m for mpi, -a for amr, -r for remeshing, -f for fast summation, -i = rh4/gv for initial condition,
        // -c = 1/2/3/4/5/6/7/8 for initial particle refinement levels, -t = end time in days, -l = icos levels for fast sum,
        // -x for only one time step, -h = 0.5/0.7/0.9 for MAC value, -b = 3/4/5 for amr levels;
        switch (opt) {
            case 'm':
                config1.use_mpi = true;
                break;
            case 'a':
                config1.use_amr = true;
                break;
            case 'r':
                config1.use_remesh = true;
                break;
            case 'f':
                config1.use_fast = true;
                break;
            case 'i': // initial condition
                config1.init_cond = string(optarg);
                break;
            case 'c':
                config1.start_level = stoi(optarg);
                config1.point_count = pow(4, config1.start_level) * 10 + 2;
                config1.tri_count = 2 * (config1.point_count - 2);
                break;
            case 't':
                config1.end_time = stoi(optarg);
                break;
            case 'l':
                config1.tree_levels = stoi(optarg);
                break;
            case 'x':
                config1.testing = true;
                break;
            case 'h':
                config1.theta = atof(optarg);
                break;
            case 'b':
                config1.max_level = config1.start_level + stoi(optarg);
                config1.max_count = pow(4, config1.max_level) * 10 + 2;
                break;
        }
    }
}

#endif
