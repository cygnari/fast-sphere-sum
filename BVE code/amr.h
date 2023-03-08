#ifndef amr_H
#define amr_H

#include <iostream>
#include <cmath>
#include <array>
#include <vector>
#include <Accelerate/Accelerate.h>
#include <cassert>
#include <algorithm>
#include "helpers.h"

int amr(vector<double>& curr_state, vector<vector<int>>& triangles, vector<vector<int>>& vert_tris, vector<double>& areas, vector<vector<int>>& parent_verts, int tri_count, int point_count, int max_points) { // adaptive mesh refinement
    double circulation_tol = 0.0025 / 4.0; // for 2562 particles, reduce for more starting particles
    double vorticity_tol = 0.2 / 2.0; // for 2562 particles, reduce for more starting particles
    int iv1, iv2, iv3, iv12, iv23, iv31,  itriv1v12v31, itriv2v12v23, itriv3v31v23, itriv12v23v31, old_tri_count;
    double vor1, vor2, vor3, max_val, min_val, vor1n, vor2n, vor3n, circulation, tri_area;
    vector<double> v1, v2, v3, v12, v23, v31, p1, p2, p3;
    vector<int> return_values {tri_count, point_count}; // return new tri_count and point_count;
    vector<int> iv1tris, iv2tris, iv3tris;
    if (point_count >= max_points) return point_count;
    else {
        old_tri_count = tri_count;
        for (int i = 0; i < old_tri_count; i++) { // work on each triangle
            // cout << "triangle: " << i << endl;
            iv1 = triangles[i][0];
            iv2 = triangles[i][1];
            iv3 = triangles[i][2];
            vor1 = curr_state[5 * iv1 + 3];
            vor2 = curr_state[5 * iv2 + 3];
            vor3 = curr_state[5 * iv3 + 3];
            max_val = max(vor1, max(vor2, vor3));
            min_val = min(vor1, min(vor2, vor3));
            // tri_area = (areas[iv1] + areas[iv2] + areas[iv3]) / 6.0;
            p1 = slice(curr_state, 5 * iv1, 1, 3);
            p2 = slice(curr_state, 5 * iv2, 1, 3);
            p3 = slice(curr_state, 5 * iv3, 1, 3);
            tri_area = sphere_tri_area(p1, p2, p3, 1.0);
            circulation =  tri_area * abs(vor1 + vor2 + vor3) / 3.0;
            if ((((max_val - min_val) > vorticity_tol) or (circulation > circulation_tol)) and (point_count < max_points)) { // if over threshold, refine
                // cout << "amr" << endl;
                areas[iv1] -= tri_area / 6.0;
                areas[iv2] -= tri_area / 6.0;
                areas[iv3] -= tri_area / 6.0;
                v1 = slice(curr_state, 5 * iv1, 1, 5);
                v2 = slice(curr_state, 5 * iv2, 1, 5);
                v3 = slice(curr_state, 5 * iv3, 1, 5);
                v12 = v1;
                v23 = v2;
                v31 = v3;
                vec_add(v12, v2);
                vec_add(v23, v3);
                vec_add(v31, v1);
                scalar_mult(v12, 0.5);
                scalar_mult(v23, 0.5);
                scalar_mult(v31, 0.5);
                // cout << v12[0] << " " << curr_state[0] << endl;
                // cout << "point_count " << point_count << endl;
                // cout << 1 << endl;
                // iv12 = check_in_vec2(curr_state, v12, point_count);
                iv12 = check_in_vec3(parent_verts, min(iv1, iv2), max(iv1, iv2));
                iv23 = check_in_vec3(parent_verts, min(iv3, iv2), max(iv3, iv2));
                iv31 = check_in_vec3(parent_verts, min(iv1, iv3), max(iv1, iv3));
                // cout << 2 << endl;
                // iv23 = check_in_vec2(curr_state, v23, point_count);
                // cout << 3 << endl;
                // iv31 = check_in_vec2(curr_state, v31, point_count);
                // cout << 4 << endl;
                if (iv12 == -1) { // v12 is a new point
                    iv12 = point_count;
                    curr_state.insert(curr_state.end(), v12.begin(), v12.end());
                    point_count += 1;
                    areas.insert(areas.end(), tri_area / 6.0);
                    vert_tris.push_back(vector<int> (0, 0));
                    parent_verts[iv12][0] = min(iv1, iv2);
                    parent_verts[iv12][1] = max(iv1, iv2);
                    // if (iv12 == 3907) {
                    //     cout << "iv12 3907: " << iv1 << " " << iv2 << endl;
                    // }
                    // if (iv12 == 3451) {
                    //     cout << "iv12 3451: " << iv1 << " " << iv2 << endl;
                    // }
                } else { // v12 is not new
                    areas[iv12] += tri_area / 6.0;
                }
                if (iv23 == -1) { // v23 is a new point
                    iv23 = point_count;
                    curr_state.insert(curr_state.end(), v23.begin(), v23.end());
                    point_count += 1;
                    areas.insert(areas.end(), tri_area / 6.0);
                    vert_tris.push_back(vector<int> (0, 0));
                    parent_verts[iv23][0] = min(iv2, iv3);
                    parent_verts[iv23][1] = max(iv2, iv3);
                    // if (iv23 == 3907) {
                    //     cout << "iv23 3907: " << iv2 << " " << iv3 << endl;
                    // }
                    // if (iv23 == 3451) {
                    //     cout << "iv23 3451: " << iv2 << " " << iv3 << endl;
                    //     cout << "v23: " << v23[0] << " " << v23[1] << " " << v23[2] << endl;
                    //     cout << "v2: " << v2[0] << " " << v2[1] << " " << v2[2] << endl;
                    //     cout << "v3: " << v3[0] << " " << v3[1] << " " << v3[2] << endl;
                    // }
                } else { // v23 is not a new point
                    areas[iv23] += tri_area / 6.0;
                }
                if (iv31 == -1) { // v31 is a new point
                    iv31 = point_count;
                    curr_state.insert(curr_state.end(), v31.begin(), v31.end());
                    point_count += 1;
                    areas.insert(areas.end(), tri_area / 6.0);
                    vert_tris.push_back(vector<int> (0, 0));
                    parent_verts[iv31][0] = min(iv1, iv3);
                    parent_verts[iv31][1] = max(iv1, iv3);
                    // if (iv31 == 3907) {
                    //     cout << "iv31 3907: " << iv1 << " " << iv3 << endl;
                    //     cout << "v31: " << v31[0] << " " << v31[1] << " " << v31[2] << endl;
                    //     cout << "v1: " << v1[0] << " " << v1[1] << " " << v1[2] << endl;
                    //     cout << "v3: " << v3[0] << " " << v3[1] << " " << v3[2] << endl;
                    //     cout << "diff: " << v31[0] - curr_state[5 * 3451] << " " << v31[1] - curr_state[5 * 3451 + 1] << " " << v31[2] - curr_state[5 * 3451 + 2] << endl;
                    //     // cout << "diff2: " <<
                    // }
                    // if (iv31 == 3451) {
                    //     cout << "iv31 3451: " << iv1 << " " << iv3 << endl;
                    // }
                }  else {
                    areas[iv31] += tri_area / 6.0;
                }
                // cout << 5 << endl;
                triangles[i] = {iv1, iv12, iv31};
                itriv1v12v31 = i;
                triangles.push_back({iv2, iv12, iv23});
                itriv2v12v23 = tri_count;
                tri_count += 1;
                triangles.push_back({iv3, iv31, iv23});
                itriv3v31v23 = tri_count;
                tri_count += 1;
                triangles.push_back({iv12, iv23, iv31});
                itriv12v23v31 = tri_count;
                tri_count += 1;
                // cout << 6 << endl;
                // cout << "iv31 " << iv31 << " new triangles " << itriv1v12v31 << " " << itriv3v31v23 << " " << itriv12v23v31 << endl;
                replace(vert_tris[iv2], i, itriv2v12v23);
                replace(vert_tris[iv3], i, itriv3v31v23);
                // cout << 7 << endl;
                vert_tris[iv12].insert(vert_tris[iv12].end(), i);
                // cout << 7.1 << endl;
                vert_tris[iv12].insert(vert_tris[iv12].end(), itriv2v12v23);
                // cout << 7.2 << endl;
                vert_tris[iv12].insert(vert_tris[iv12].end(), itriv12v23v31);
                // cout << 7.3 << endl;
                vert_tris[iv23].insert(vert_tris[iv23].end(), itriv2v12v23);
                // cout << 7.4 << endl;
                vert_tris[iv23].insert(vert_tris[iv23].end(), itriv3v31v23);
                vert_tris[iv23].insert(vert_tris[iv23].end(), itriv12v23v31);
                vert_tris[iv31].insert(vert_tris[iv31].end(), itriv1v12v31);
                vert_tris[iv31].insert(vert_tris[iv31].end(), itriv3v31v23);
                vert_tris[iv31].insert(vert_tris[iv31].end(), itriv12v23v31);
                // cout << 8 << endl;
            } else { // do not refine if not over threshold
                continue;
            }
        }
    }
    // return_values = {tri_count, point_count};
    // return return_values;
    return point_count;
}

#endif
