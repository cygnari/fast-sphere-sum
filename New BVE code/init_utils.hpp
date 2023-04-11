#ifndef init_H
#define init_H

#include <cmath>
#include <vector>
#include "general_utils.hpp"
#include "structs.hpp"

void dynamics_points_initialize(run_config& run_information, vector<double>& dynamics_state, vector<vector<int>>& dynamics_triangles) {
    // creates all the dynamics points and corresponding triangles
    double phi = (1 + sqrt(5)) / 2;
    dynamics_state.resize(run_information.dynamics_initial_points * run_information.info_per_point, 0);
    dynamics_triangles.resize(run_information.dynamics_initial_triangles, vector<int> (4, 0));
    run_information.dynamics_curr_point_count = 12;
    vector_copy2(dynamics_state, project_to_sphere_2(vector<double> {0, 1, phi}, run_information.radius), 0, 3);
    vector_copy2(dynamics_state, project_to_sphere_2(vector<double> {0, -1, phi}, run_information.radius), run_information.info_per_point, 3);
    vector_copy2(dynamics_state, project_to_sphere_2(vector<double> {0, 1, -phi}, run_information.radius), 2 * run_information.info_per_point, 3);
    vector_copy2(dynamics_state, project_to_sphere_2(vector<double> {0, -1, -phi}, run_information.radius), 3 * run_information.info_per_point, 3);
    vector_copy2(dynamics_state, project_to_sphere_2(vector<double> {1, phi, 0}, run_information.radius), 4 * run_information.info_per_point, 3);
    vector_copy2(dynamics_state, project_to_sphere_2(vector<double> {1, -phi, 0}, run_information.radius), 5 * run_information.info_per_point, 3);
    vector_copy2(dynamics_state, project_to_sphere_2(vector<double> {-1, phi, 0}, run_information.radius), 6 * run_information.info_per_point, 3);
    vector_copy2(dynamics_state, project_to_sphere_2(vector<double> {-1, -phi, 0}, run_information.radius), 7 * run_information.info_per_point, 3);
    vector_copy2(dynamics_state, project_to_sphere_2(vector<double> {phi, 0, 1}, run_information.radius), 8 * run_information.info_per_point, 3);
    vector_copy2(dynamics_state, project_to_sphere_2(vector<double> {phi, 0, -1}, run_information.radius), 9 * run_information.info_per_point, 3);
    vector_copy2(dynamics_state, project_to_sphere_2(vector<double> {-phi, 0, 1}, run_information.radius), 10 * run_information.info_per_point, 3);
    vector_copy2(dynamics_state, project_to_sphere_2(vector<double> {-phi, 0, -1}, run_information.radius), 11 * run_information.info_per_point, 3);
    dynamics_triangles[0].insert(dynamics_triangles[0].begin(), {1, 2, 9, 0}); // 0, 1, 2 are indices of the three vertices
    dynamics_triangles[1].insert(dynamics_triangles[1].begin(), {1, 2, 11, 0}); // 20 starting faces
    dynamics_triangles[2].insert(dynamics_triangles[2].begin(), {1, 5, 7, 0});
    dynamics_triangles[3].insert(dynamics_triangles[3].begin(), {1, 5, 9, 0});
    dynamics_triangles[4].insert(dynamics_triangles[4].begin(), {1, 7, 11, 0});
    dynamics_triangles[5].insert(dynamics_triangles[5].begin(), {2, 6, 8, 0});
    dynamics_triangles[6].insert(dynamics_triangles[6].begin(), {2, 6, 9, 0});
    dynamics_triangles[7].insert(dynamics_triangles[7].begin(), {2, 8, 11, 0});
    dynamics_triangles[8].insert(dynamics_triangles[8].begin(), {3, 4, 10, 0});
    dynamics_triangles[9].insert(dynamics_triangles[9].begin(), {3, 4, 12, 0});
    dynamics_triangles[10].insert(dynamics_triangles[10].begin(), {3, 5, 7, 0});
    dynamics_triangles[11].insert(dynamics_triangles[11].begin(), {3, 5, 10, 0});
    dynamics_triangles[12].insert(dynamics_triangles[12].begin(), {3, 7, 12, 0});
    dynamics_triangles[13].insert(dynamics_triangles[13].begin(), {4, 6, 8, 0});
    dynamics_triangles[14].insert(dynamics_triangles[14].begin(), {4, 6, 10, 0});
    dynamics_triangles[15].insert(dynamics_triangles[15].begin(), {4, 8, 12, 0});
    dynamics_triangles[16].insert(dynamics_triangles[16].begin(), {5, 9, 10, 0});
    dynamics_triangles[17].insert(dynamics_triangles[17].begin(), {6, 9, 10, 0});
    dynamics_triangles[18].insert(dynamics_triangles[18].begin(), {7, 11, 12, 0});
    dynamics_triangles[19].insert(dynamics_triangles[19].begin(), {8, 11, 12, 0});

}

void vorticity_initialize(run_config& run_information, vector<double>& dynamics_state) {
    // initializes the initial vorticity
}

void tracer_initialize(run_config& run_information, vector<double>& dynamics_state) {
    // initializes the tracer
}

#endif
