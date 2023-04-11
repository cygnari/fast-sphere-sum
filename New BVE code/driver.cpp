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

using namespace std;

int main() {
  
    run_config run_information;
    read_run_config("namelist.txt", run_information); // reads in run configuration information

    vector<double> dynamics_state; // list of points and other information in a flattened array
    vector<vector<int>> dynamics_triangles; // each entry is a vector which contains the 3 vertices and the refinement level of the triangle
    dynamics_points_initialize(run_information, dynamics_state, dynamics_triangles);

    return 0;
}
