#ifndef icos_H
#define icos_H

#include <iostream>
#include <cmath>
#include <array>
#include <vector>
#include <Accelerate/Accelerate.h>
#include <cassert>
#include <algorithm>
#include "helpers.h"

struct icos_struct {
    int levels;
    double radius;
    vector<vector<double>> verts; // list of icosahedron vertices
    vector<vector<vector<double>>> tri_info; // information about each triangle forming the faces
    vector<vector<vector<int>>> tri_verts; // length 3 int vector listing the vertex indices
    vector<vector<vector<int>>> tri_points; // points inside each triangle
};

void icos_init(icos_struct& icos, int levels, double radius) {
    vector<vector<double>> verts;
    vector<vector<vector<double>>> tri_info;
    vector<vector<vector<int>>> tri_verts;

    double phi = (1 + sqrt(5)) / 2;
    vector<double> center, v1, v2, v3, v12, v23, v31;
    int iv1, iv2, iv3, iv12, iv23, iv13;
    verts.push_back(project_to_sphere_2(vector<double> {0, 1, phi}, radius)); // 12 starting points
    verts.push_back(project_to_sphere_2(vector<double> {0, -1, phi}, radius));
    verts.push_back(project_to_sphere_2(vector<double> {0, 1, -phi}, radius));
    verts.push_back(project_to_sphere_2(vector<double> {0, -1, -phi}, radius));
    verts.push_back(project_to_sphere_2(vector<double> {1, phi, 0}, radius));
    verts.push_back(project_to_sphere_2(vector<double> {1, -phi, 0}, radius));
    verts.push_back(project_to_sphere_2(vector<double> {-1, phi, 0}, radius));
    verts.push_back(project_to_sphere_2(vector<double> {-1, -phi, 0}, radius));
    verts.push_back(project_to_sphere_2(vector<double> {phi, 0, 1}, radius));
    verts.push_back(project_to_sphere_2(vector<double> {phi, 0, -1}, radius));
    verts.push_back(project_to_sphere_2(vector<double> {-phi, 0, 1}, radius));
    verts.push_back(project_to_sphere_2(vector<double> {-phi, 0, -1}, radius));
    tri_verts.push_back(vector<vector<int>> (20, vector<int> (0)));
    tri_info.push_back(vector<vector<double>> (20, vector<double> (0)));
    tri_verts[0][0].insert(tri_verts[0][0].end(), {1, 2, 9}); // 0, 1, 2 are indices of the three vertices
    tri_verts[0][1].insert(tri_verts[0][1].end(), {1, 2, 11}); // 20 starting faces
    tri_verts[0][2].insert(tri_verts[0][2].end(), {1, 5, 7});
    tri_verts[0][3].insert(tri_verts[0][3].end(), {1, 5, 9});
    tri_verts[0][4].insert(tri_verts[0][4].end(), {1, 7, 11});
    tri_verts[0][5].insert(tri_verts[0][5].end(), {2, 6, 8});
    tri_verts[0][6].insert(tri_verts[0][6].end(), {2, 6, 9});
    tri_verts[0][7].insert(tri_verts[0][7].end(), {2, 8, 11});
    tri_verts[0][8].insert(tri_verts[0][8].end(), {3, 4, 10});
    tri_verts[0][9].insert(tri_verts[0][9].end(), {3, 4, 12});
    tri_verts[0][10].insert(tri_verts[0][10].end(), {3, 5, 7});
    tri_verts[0][11].insert(tri_verts[0][11].end(), {3, 5, 10});
    tri_verts[0][12].insert(tri_verts[0][12].end(), {3, 7, 12});
    tri_verts[0][13].insert(tri_verts[0][13].end(), {4, 6, 8});
    tri_verts[0][14].insert(tri_verts[0][14].end(), {4, 6, 10});
    tri_verts[0][15].insert(tri_verts[0][15].end(), {4, 8, 12});
    tri_verts[0][16].insert(tri_verts[0][16].end(), {5, 9, 10});
    tri_verts[0][17].insert(tri_verts[0][17].end(), {6, 9, 10});
    tri_verts[0][18].insert(tri_verts[0][18].end(), {7, 11, 12});
    tri_verts[0][19].insert(tri_verts[0][19].end(), {8, 11, 12});
    for (int i = 0; i < 20; i++) { // info about the first 20 faces
        for (int j = 0; j < 3; j++) tri_verts[0][i][j] -= 1;
        iv1 = tri_verts[0][i][0];
        iv2 = tri_verts[0][i][1];
        iv3 = tri_verts[0][i][2];
        v1 = verts[iv1];
        v2 = verts[iv2];
        v3 = verts[iv3];
        center = tri_center(v1, v2, v3, radius);
        tri_info[0][i].insert(tri_info[0][i].end(), center.begin(), center.end()); // index 0 1 2 is the triangle center
        tri_info[0][i].push_back(tri_radius(v1, v2, v3, center)); // index 3 is the triangle radius
    }
    for (int i = 0; i < levels; i++) { // iterative refinement
        tri_info.push_back(vector<vector<double>> (20 * pow(4, i + 1), vector<double> (0)));
        tri_verts.push_back(vector<vector<int>> (20 * pow(4, i + 1), vector<int> (0)));
        for (int j = 0; j < 20 * pow(4, i); j++) {
            iv1 = tri_verts[i][j][0];
            iv2 = tri_verts[i][j][1];
            iv3 = tri_verts[i][j][2];
            v1 = verts[iv1];
            v2 = verts[iv2];
            v3 = verts[iv3];
            v12 = v1;
            v23 = v2;
            v31 = v3;
            vec_add(v12, v2); // v12 halfway between v1 and v2
            vec_add(v23, v3);
            vec_add(v31, v1);
            scalar_mult(v12, 0.5);
            scalar_mult(v23, 0.5);
            scalar_mult(v31, 0.5);
            project_to_sphere(v12, radius);
            project_to_sphere(v23, radius);
            project_to_sphere(v31, radius);
            iv12 = check_in_vec(verts, v12); // check if v12 already exists
            iv13 = check_in_vec(verts, v31);
            iv23 = check_in_vec(verts, v23);
            if (iv12 == -1) {
                iv12 = verts.size();
                verts.push_back(v12);
            }
            if (iv13 == -1) {
                iv13 = verts.size();
                verts.push_back(v31);
            }
            if (iv23 == -1) {
                iv23 = verts.size();
                verts.push_back(v23);
            }
            tri_verts[i+1][4*j].insert(tri_verts[i+1][4*j].end(), {iv1, iv13, iv12}); // 4 children triangles
            tri_verts[i+1][4*j+1].insert(tri_verts[i+1][4*j+1].end(),{iv3, iv23, iv13});
            tri_verts[i+1][4*j+2].insert(tri_verts[i+1][4*j+2].end(),{iv2, iv12, iv23});
            tri_verts[i+1][4*j+3].insert(tri_verts[i+1][4*j+3].end(),{iv12, iv13, iv23});

            center = tri_center(v1, v12, v31, radius);
            tri_info[i+1][4*j].insert(tri_info[i+1][4*j].end(), center.begin(), center.end());
            tri_info[i+1][4*j].push_back(tri_radius(v1, v12, v31, center));

            center = tri_center(v3, v23, v31, radius);
            tri_info[i+1][4*j+1].insert(tri_info[i+1][4*j+1].end(), center.begin(), center.end());
            tri_info[i+1][4*j+1].push_back(tri_radius(v3, v23, v31, center));

            center = tri_center(v2, v12, v23, radius);
            tri_info[i+1][4*j+2].insert(tri_info[i+1][4*j+2].end(), center.begin(), center.end());
            tri_info[i+1][4*j+2].push_back(tri_radius(v2, v12, v23, center));

            center = tri_center(v12, v31, v23, radius);
            tri_info[i+1][4*j+3].insert(tri_info[i+1][4*j+3].end(), center.begin(), center.end());
            tri_info[i+1][4*j+3].push_back(tri_radius(v12, v31, v23, center));
        }
    }
    icos.verts = verts;
    icos.tri_info = tri_info;
    icos.tri_verts = tri_verts;
}

icos_struct icosahedron_init(int levels, double radius) {
    icos_struct icos1;
    icos1.levels = levels;
    icos1.radius = radius;
    icos_init(icos1, levels, radius);
    return icos1;
}

void point_assign(icos_struct& icos, vector<double>& point, vector<vector<int>>& point_locs, int point_id) {
    int iv1, iv2, iv3, lb, ub;
    vector<double> v1, v2, v3;
    for (int i = 0; i < icos.levels; i++) {
        if (i > 0) {
            lb = 4 * point_locs[i-1][point_id]; // utilize tree structure to minimize searching
            ub = lb + 4;
        } else {
            lb = 0;
            ub = 20;
        }
        for (int j = lb; j < ub; j++) {
            iv1 = icos.tri_verts[i][j][0];
            iv2 = icos.tri_verts[i][j][1];
            iv3 = icos.tri_verts[i][j][2];
            v1 = icos.verts[iv1];
            v2 = icos.verts[iv2];
            v3 = icos.verts[iv3];
            if (check_in_tri(v1, v2, v3, point)) {
                point_locs[i][point_id] = j;
                icos.tri_points[i][j].push_back(point_id);
                break;
            }
        }
    }
}

void points_assign(icos_struct& icos, vector<double>& points, vector<vector<int>>& point_locs, int point_count) {
// void points_assign(vector<vector<vector<int>>>& tri_verts, vector<vector<double>>& verts, vector<double>& points, vector<vector<vector<int>>>& tri_points, vector<vector<int>>& point_locs, int levels, int point_count) { // finds which icosahedron triangle each dynamics point is in
    vector<vector<vector<int>>> tri_points (icos.levels);
    icos.tri_points = tri_points;
    vector<double> point;
    // vector<double> v1, v2, v3, point;
    // int iv1, iv2, iv3, lb, ub;
    // point_locs.clear();
    // tri_points.clear();
    // point_locs.push_back(vector<int> (point_count, 0));
    // icos.tri_points[0] = vector<vector<int>> (20);
    for (int i = 0; i < icos.levels; i++) {
        point_locs.push_back(vector<int> (point_count, 0));
        icos.tri_points[i] = vector<vector<int>> (20 * pow(4, i));
    }
    for (int i = 0; i < point_count; i++) {
        point = slice(points, 5 * i, 1, 3);
        point_assign(icos, point, point_locs, i);
    }
    // for (int i = 0; i < point_count; i++) { // finds the placement of each point in the initial icosahedron
    //     point = slice(points, 5 * i, 1, 3);
    //     for (int j = 0; j < 20; j++) {
    //         iv1 = icos.tri_verts[0][j][0];
    //         iv2 = icos.tri_verts[0][j][1];
    //         iv3 = icos.tri_verts[0][j][2];
    //         v1 = icos.verts[iv1];
    //         v2 = icos.verts[iv2];
    //         v3 = icos.verts[iv3];
    //         if (check_in_tri(v1, v2, v3, point)) {
    //             point_locs[0][i] = j;
    //             icos.tri_points[0][j].push_back(i);
    //             break;
    //         }
    //     }
    // }
    // for (int i = 1; i < icos.levels; i++) { // finds point loc in refined faces
    //     // point_locs.push_back(vector<int> (point_count, 0));
    //     // tri_points[i] = vector<vector<int>> (20 * pow(4, i));
    //     for (int j = 0; j < point_count; j++) {
    //         point = slice(points, 5 * j, 1, 3);
    //         lb = 4 * point_locs[i-1][j]; // utilize tree structure to minimize searching
    //         ub = lb + 4;
    //         for (int k = lb; k < ub; k++) {
    //             iv1 = icos.tri_verts[i][k][0];
    //             iv2 = icos.tri_verts[i][k][1];
    //             iv3 = icos.tri_verts[i][k][2];
    //             v1 = icos.verts[iv1];
    //             v2 = icos.verts[iv2];
    //             v3 = icos.verts[iv3];
    //             if (check_in_tri(v1, v2, v3, point)) {
    //                 point_locs[i][j] = k;
    //                 icos.tri_points[i][k].push_back(j);
    //                 break;
    //             }
    //         }
    //     }
    // }
}

#endif
