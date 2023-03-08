#ifndef remesh_H
#define remesh_H

#include <iostream>
#include <cmath>
#include <array>
#include <vector>
#include <Accelerate/Accelerate.h>
#include <cassert>
#include <algorithm>
#include "helpers.h"

void regrid_point(vector<double>& curr_state, vector<double>& target_points, vector<vector<int>>& triangles, vector<vector<int>>& vert_tris, int point_count, int tri_count, double omega, int ID, vector<int> particle_thread) { // remesh back to original particle locations
    vector<double> curr_target, curr_pos, v1, v2, v3; // need to fix vorticity remeshing
    vector<int> poss_tris;
    // double curr_vor;
    int test_count, iv1, iv2, iv3;
    int success = 0, failure = 0, bad = 0;
    bool found;
    // for (int i = lb; i < ub; i++) {
    for (int i = 0; i < point_count; i++) {
        if (particle_thread[i] == ID) {
            // cout << "regrid: " << i << endl;
            found = false;
            curr_target = slice(target_points, 5 * i, 1, 3);
            curr_pos = slice(curr_state, 5 * i, 1, 3);
            // curr_vor = curr_state[5 * i + 3];
            poss_tris = vert_tris[i];
            test_count = poss_tris.size();
            // cout << "test count " << test_count << endl;
            // for (int j = 0; j < test_count; j++) cout << poss_tris[j] << endl;
            for (int j = 0; j < test_count; j++) { // first check triangles adjacent to initial point
                // cout << "test triangle " << j << endl;
                // cout << "possible triangle j " << poss_tris[j] << endl;
                iv1 = triangles[poss_tris[j]][0];
                // cout << "iv1 " << iv1 << endl;
                iv2 = triangles[poss_tris[j]][1];
                // cout << "iv2 " << iv2 << endl;
                iv3 = triangles[poss_tris[j]][2];
                // cout << "iv3 " << iv3 << endl;
                v1 = slice(curr_state, 5 * iv1, 1, 3);
                v2 = slice(curr_state, 5 * iv2, 1, 3);
                v3 = slice(curr_state, 5 * iv3, 1, 3);
                if (check_in_tri(v1, v2, v3, curr_target)) {
                    tri_interp(iv1, iv2, iv3, v1, v2, v3, curr_state, target_points, curr_target, i, omega);
                    success += 1;
                    found = true;
                    break;
                }
            }
            if (found) {
                continue;
            }
            failure += 1;
            for (int j = 0; j < tri_count; j++) { // if point not found in adjacent triangles, search all of them
                iv1 = triangles[j][0];
                iv2 = triangles[j][1];
                iv3 = triangles[j][2];
                v1 = slice(curr_state, 5 * iv1, 1, 3);
                v2 = slice(curr_state, 5 * iv2, 1, 3);
                v3 = slice(curr_state, 5 * iv3, 1, 3);
                if (check_in_tri(v1, v2, v3, curr_target)) {
                    tri_interp(iv1, iv2, iv3, v1, v2, v3, curr_state, target_points, curr_target, i, omega);
                    found = true;
                    // cout << "point: " << i << " in triangle: " << j << " bary cords: " << bary[0] << " " << bary[1] << " " << bary[2] << endl;
                    // cout << "curr vor: " << curr_vor << " vert vors " << curr_state[5 * iv1 + 3] << " " << curr_state[5 * iv2 + 3] << " " << curr_state[5 * iv3 + 3] <<  endl;
                    break;
                }
            }
            if (not found) { // hopefully not
                bad += 1;
            }
        }
    }
    // cout << "success: " << success << " failure: " << failure << " bad: " << bad << endl;
    // cout << "success: " << success << " failure: " << failure << " bad: " << bad << " ID " << ID << endl;
}

// void regrid_points(vector<double>& curr_state, vector<double>& target_points, vector<vector<int>>& triangles, vector<vector<int>>& vert_tris, int point_count, int tri_count, double omega, int lb, int ub, int ID) { // remesh back to original particle locations
//     vector<double> curr_target, curr_pos, v1, v2, v3; // need to fix vorticity remeshing
//     vector<int> poss_tris;
//     // double curr_vor;
//     int test_count, iv1, iv2, iv3;
//     int success = 0, failure = 0, bad = 0;
//     bool found;
//     for (int i = lb; i < ub; i++) {
//         // cout << "regrid: " << i << endl;
//         found = false;
//         curr_target = slice(target_points, 5 * i, 1, 3);
//         curr_pos = slice(curr_state, 5 * i, 1, 3);
//         // curr_vor = curr_state[5 * i + 3];
//         poss_tris = vert_tris[i];
//         test_count = poss_tris.size();
//         // cout << "test count " << test_count << endl;
//         // for (int j = 0; j < test_count; j++) cout << poss_tris[j] << endl;
//         for (int j = 0; j < test_count; j++) { // first check triangles adjacent to initial point
//             // cout << "test triangle " << j << endl;
//             // cout << "possible triangle j " << poss_tris[j] << endl;
//             iv1 = triangles[poss_tris[j]][0];
//             // cout << "iv1 " << iv1 << endl;
//             iv2 = triangles[poss_tris[j]][1];
//             // cout << "iv2 " << iv2 << endl;
//             iv3 = triangles[poss_tris[j]][2];
//             // cout << "iv3 " << iv3 << endl;
//             v1 = slice(curr_state, 5 * iv1, 1, 3);
//             v2 = slice(curr_state, 5 * iv2, 1, 3);
//             v3 = slice(curr_state, 5 * iv3, 1, 3);
//             if (check_in_tri(v1, v2, v3, curr_target)) {
//                 tri_interp(iv1, iv2, iv3, v1, v2, v3, curr_state, target_points, curr_target, i, omega);
//                 success += 1;
//                 found = true;
//                 break;
//             }
//         }
//         if (found) {
//             continue;
//         }
//         failure += 1;
//         for (int j = 0; j < tri_count; j++) { // if point not found in adjacent triangles, search all of them
//             iv1 = triangles[j][0];
//             iv2 = triangles[j][1];
//             iv3 = triangles[j][2];
//             v1 = slice(curr_state, 5 * iv1, 1, 3);
//             v2 = slice(curr_state, 5 * iv2, 1, 3);
//             v3 = slice(curr_state, 5 * iv3, 1, 3);
//             if (check_in_tri(v1, v2, v3, curr_target)) {
//                 tri_interp(iv1, iv2, iv3, v1, v2, v3, curr_state, target_points, curr_target, i, omega);
//                 found = true;
//                 // cout << "point: " << i << " in triangle: " << j << " bary cords: " << bary[0] << " " << bary[1] << " " << bary[2] << endl;
//                 // cout << "curr vor: " << curr_vor << " vert vors " << curr_state[5 * iv1 + 3] << " " << curr_state[5 * iv2 + 3] << " " << curr_state[5 * iv3 + 3] <<  endl;
//                 break;
//             }
//         }
//         if (not found) { // hopefully not
//             bad += 1;
//         }
//     }
//     // cout << "success: " << success << " failure: " << failure << " bad: " << bad << endl;
//     // cout << "success: " << success << " failure: " << failure << " bad: " << bad << " ID " << ID << endl;
// }

#endif
