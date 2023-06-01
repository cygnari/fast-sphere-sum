#ifndef amr_H
#define amr_H

#include <cmath>
#include <vector>
#include <iostream>
#include "general_utils.hpp"
#include "structs.hpp"
#include "interp_utils.hpp"
#include "init_utils.hpp"

void amr(run_config& run_information, vector<double>& new_dynamics_state, vector<double>& old_dynamics_state,
        vector<vector<vector<int>>>& new_dynamics_triangles, vector<vector<vector<int>>>& old_dynamics_triangles,
        vector<vector<bool>>& new_dynamics_triangles_is_leaf, vector<vector<bool>>& old_dynamics_triangles_is_leaf,
        vector<double>& dynamics_areas, double omega) {

    int iv, iv1, iv2, iv3, iv12, iv23, iv31, iv1n, iv2n, iv3n, iv4n, iv5n, iv6n, tri_level, tri_index, super_tri_index;
    vector<double> v, v1, v2, v3, v12, v23, v31, v11, v22, v33;
    double vormin, vormax, vor1, vor2, vor3, tri_area, vor;
    int old_point_count = run_information.dynamics_curr_point_count;

    new_dynamics_state = slice(old_dynamics_state, 0, 1, run_information.info_per_point * run_information.dynamics_initial_points);
    new_dynamics_triangles = old_dynamics_triangles;
    new_dynamics_triangles.resize(run_information.dynamics_levels_min);
    new_dynamics_triangles.resize(run_information.dynamics_levels_max);
    dynamics_areas.clear();
    area_initialize(run_information, new_dynamics_state, new_dynamics_triangles, dynamics_areas);
    new_dynamics_triangles_is_leaf.resize(run_information.dynamics_levels_min);

    run_information.dynamics_curr_point_count = run_information.dynamics_initial_points;
    run_information.dynamics_curr_tri_count = run_information.dynamics_initial_triangles;
    new_dynamics_triangles_is_leaf[run_information.dynamics_levels_min-1] = vector<bool> (run_information.dynamics_initial_triangles, true);
    new_dynamics_triangles_is_leaf.resize(run_information.dynamics_levels_max);
    for (int i = run_information.dynamics_levels_min; i < run_information.dynamics_levels_max; i++) {
        new_dynamics_triangles_is_leaf[i] = vector<bool> (20 * pow(4, i), false);
        new_dynamics_triangles[i] = vector<vector<int>> (20 * pow(4, i), vector<int> (4, 0));
    }

    for (int i = run_information.dynamics_levels_min - 1; i < run_information.dynamics_levels_max - 1; i++) {
        for (int j = 0; j < 20 * pow(4, i); j++) {
            if (new_dynamics_triangles_is_leaf[i][j]) {
                iv1 = new_dynamics_triangles[i][j][0];
                iv2 = new_dynamics_triangles[i][j][1];
                iv3 = new_dynamics_triangles[i][j][2];
                v1 = slice(new_dynamics_state, run_information.info_per_point * iv1, 1, 3);
                v2 = slice(new_dynamics_state, run_information.info_per_point * iv2, 1, 3);
                v3 = slice(new_dynamics_state, run_information.info_per_point * iv3, 1, 3);
                tri_area = sphere_tri_area(v1, v2, v3, run_information.radius);
                vor1 = new_dynamics_state[run_information.info_per_point * iv1 + 3];
                vor2 = new_dynamics_state[run_information.info_per_point * iv2 + 3];
                vor3 = new_dynamics_state[run_information.info_per_point * iv3 + 3];
                vor = abs((vor1 + vor2 + vor3) / 3.0);
                vormax = max(vor1, max(vor2, vor3));
                vormin = min(vor1, min(vor2, vor3));
                if ((tri_area * vor > run_information.amr_circ_thresh) or (vormax - vormin > run_information.amr_vor_thresh)) {
                    // refine
                    dynamics_areas[iv1] -= tri_area / 6.0;
                    dynamics_areas[iv2] -= tri_area / 6.0;
                    dynamics_areas[iv3] -= tri_area / 6.0;

                    new_dynamics_triangles_is_leaf[i][j] = false;

                    v12 = v1;
                    v23 = v2;
                    v31 = v3;
                    vec_add(v12, v2);
                    vec_add(v23, v3);
                    vec_add(v31, v1);
                    scalar_mult(v12, 0.5);
                    scalar_mult(v23, 0.5);
                    scalar_mult(v31, 0.5);

                    iv12 = check_point_exist2(new_dynamics_state, v12, run_information.dynamics_curr_point_count, pow(10, -14), run_information.info_per_point);
                    iv23 = check_point_exist2(new_dynamics_state, v23, run_information.dynamics_curr_point_count, pow(10, -14), run_information.info_per_point);
                    iv31 = check_point_exist2(new_dynamics_state, v31, run_information.dynamics_curr_point_count, pow(10, -14), run_information.info_per_point);

                    v11 = slice(new_dynamics_state, run_information.info_per_point * iv1, 1, run_information.info_per_point);
                    v22 = slice(new_dynamics_state, run_information.info_per_point * iv2, 1, run_information.info_per_point);
                    v33 = slice(new_dynamics_state, run_information.info_per_point * iv3, 1, run_information.info_per_point);

                    if (iv12 == -1) {
                        // not yet refined to at this time step
                        iv12 = check_point_exist2(old_dynamics_state, v12, old_point_count, pow(10, -14), run_information.info_per_point);
                        if (iv12 == -1) {
                            // did not exist at previous time step
                            tie(tri_level, tri_index) = find_leaf_tri(v12, old_dynamics_state, old_dynamics_triangles, old_dynamics_triangles_is_leaf, run_information.info_per_point, run_information.dynamics_levels_max);
                            super_tri_index = floor(tri_index / 4.0);
                            iv1n = old_dynamics_triangles[tri_level-1][super_tri_index][0];
                            iv2n = old_dynamics_triangles[tri_level-1][super_tri_index][1];
                            iv3n = old_dynamics_triangles[tri_level-1][super_tri_index][2];
                            iv4n = old_dynamics_triangles[tri_level][4*super_tri_index+3][0];
                            iv5n = old_dynamics_triangles[tri_level][4*super_tri_index+3][1];
                            iv6n = old_dynamics_triangles[tri_level][4*super_tri_index+3][2];
                            v12 = biquadratic_interp(run_information, v12, iv1n, iv2n, iv3n, iv4n, iv5n, iv6n, old_dynamics_state);
                        } else {
                            // existed at previous time step
                            v12 = slice(old_dynamics_state, run_information.info_per_point * iv12, 1, run_information.info_per_point);
                        }
                        iv12 = run_information.dynamics_curr_point_count;
                        run_information.dynamics_curr_point_count += 1;
                        new_dynamics_state.insert(new_dynamics_state.end(), v12.begin(), v12.end());
                        dynamics_areas.insert(dynamics_areas.end(), tri_area / 6.0);
                    } else {
                        // existing point
                        dynamics_areas[iv12] += tri_area / 6.0;
                    }
                    if (iv23 == -1) {
                        iv23 = check_point_exist2(old_dynamics_state, v23, old_point_count, pow(10, -14), run_information.info_per_point);
                        if (iv23 == -1) {
                            tie(tri_level, tri_index) = find_leaf_tri(v23, old_dynamics_state, old_dynamics_triangles, old_dynamics_triangles_is_leaf, run_information.info_per_point, run_information.dynamics_levels_max);
                            super_tri_index = floor(tri_index / 4.0);
                            iv1n = old_dynamics_triangles[tri_level-1][super_tri_index][0];
                            iv2n = old_dynamics_triangles[tri_level-1][super_tri_index][1];
                            iv3n = old_dynamics_triangles[tri_level-1][super_tri_index][2];
                            iv4n = old_dynamics_triangles[tri_level][4*super_tri_index+3][0];
                            iv5n = old_dynamics_triangles[tri_level][4*super_tri_index+3][1];
                            iv6n = old_dynamics_triangles[tri_level][4*super_tri_index+3][2];
                            v23 = biquadratic_interp(run_information, v23, iv1n, iv2n, iv3n, iv4n, iv5n, iv6n, old_dynamics_state);
                        } else {
                            v23 = slice(old_dynamics_state, run_information.info_per_point * iv23, 1, run_information.info_per_point);
                        }
                        iv23 = run_information.dynamics_curr_point_count;
                        run_information.dynamics_curr_point_count += 1;
                        new_dynamics_state.insert(new_dynamics_state.end(), v23.begin(), v23.end());
                        dynamics_areas.insert(dynamics_areas.end(), tri_area / 6.0);
                    } else {
                        dynamics_areas[iv23] += tri_area / 6.0;
                    }
                    if (iv31 == -1) {
                        iv31 = check_point_exist2(old_dynamics_state, v31, old_point_count, pow(10, -14), run_information.info_per_point);
                        if (iv31 == -1) {
                            tie(tri_level, tri_index) = find_leaf_tri(v31, old_dynamics_state, old_dynamics_triangles, old_dynamics_triangles_is_leaf, run_information.info_per_point, run_information.dynamics_levels_max);
                            super_tri_index = floor(tri_index / 4.0);
                            iv1n = old_dynamics_triangles[tri_level-1][super_tri_index][0];
                            iv2n = old_dynamics_triangles[tri_level-1][super_tri_index][1];
                            iv3n = old_dynamics_triangles[tri_level-1][super_tri_index][2];
                            iv4n = old_dynamics_triangles[tri_level][4*super_tri_index+3][0];
                            iv5n = old_dynamics_triangles[tri_level][4*super_tri_index+3][1];
                            iv6n = old_dynamics_triangles[tri_level][4*super_tri_index+3][2];
                            v31 = biquadratic_interp(run_information, v31, iv1n, iv2n, iv3n, iv4n, iv5n, iv6n, old_dynamics_state);
                        } else {
                            v31 = slice(old_dynamics_state, run_information.info_per_point * iv31, 1, run_information.info_per_point);
                        }
                        iv31 = run_information.dynamics_curr_point_count;
                        run_information.dynamics_curr_point_count += 1;
                        new_dynamics_state.insert(new_dynamics_state.end(), v31.begin(), v31.end());
                        dynamics_areas.insert(dynamics_areas.end(), tri_area / 6.0);
                    } else {
                        dynamics_areas[iv31] += tri_area / 6.0;
                    }
                    new_dynamics_triangles[i+1][4*j].insert(new_dynamics_triangles[i+1][4*j].begin(), {iv1, iv12, iv31, i + 1});
                    new_dynamics_triangles[i+1][4*j+1].insert(new_dynamics_triangles[i+1][4*j+1].begin(), {iv2, iv23, iv12, i + 1});
                    new_dynamics_triangles[i+1][4*j+2].insert(new_dynamics_triangles[i+1][4*j+2].begin(), {iv3, iv31, iv23, i + 1});
                    new_dynamics_triangles[i+1][4*j+3].insert(new_dynamics_triangles[i+1][4*j+3].begin(), {iv12, iv23, iv31, i + 1});
                    new_dynamics_triangles_is_leaf[i+1][4*j] = true;
                    new_dynamics_triangles_is_leaf[i+1][4*j+1] = true;
                    new_dynamics_triangles_is_leaf[i+1][4*j+2] = true;
                    new_dynamics_triangles_is_leaf[i+1][4*j+3] = true;
                    run_information.dynamics_curr_tri_count += 3;
                }
            }
        }
    }
}

void amr_wrapper(run_config& run_information, vector<double>& dynamics_state, vector<vector<vector<int>>>& dynamics_triangles,
        vector<vector<bool>>& dynamics_triangles_is_leaf, vector<double>& dynamics_areas, double omega) {
    // wraps amr routine
    vector<double> new_dynamics_state;
    vector<vector<vector<int>>> new_dynamics_triangles;
    vector<vector<bool>> new_dynamics_triangles_is_leaf;
    amr(run_information, new_dynamics_state, dynamics_state, new_dynamics_triangles, dynamics_triangles, new_dynamics_triangles_is_leaf, dynamics_triangles_is_leaf, dynamics_areas, omega);
    dynamics_state.assign(new_dynamics_state.begin(), new_dynamics_state.end());
    dynamics_triangles = new_dynamics_triangles;
    dynamics_triangles_is_leaf = new_dynamics_triangles_is_leaf;
}

#endif
