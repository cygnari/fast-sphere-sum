#ifndef run_H
#define run_H

#include "struct_list.h"

run_config config_init(int argc, char** argv, double radius) {
    int opt;
    run_config config1;
    config1.radius = radius;
    while ((opt = getopt(argc, argv, "marfi:c:")) != -1) {
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
                config1.point_count = stoi(optarg);
                break;
        }
    }
    return config1;
}

#endif