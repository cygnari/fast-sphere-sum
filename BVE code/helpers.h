// helper functions

#ifndef F_H
#define F_H

#include <iostream>
#include <cmath>
#include <array>
#include <vector>
#include <Accelerate/Accelerate.h>
#include <cassert>
#include <algorithm>

extern "C" {
    extern int dgesv_(int*,int*,double*,int*,int*,double*,int*,int*);
    extern int dgetrf_(int*,int*,double*,int*,int*,int*);
    extern int dgetrs_(char*,int*,int*,double*,int*,int*,double*,int*,int*);
}

using namespace std;

double dot_prod(vector<double>& x, vector<double>& y) { // dot product of vectors x and y
    double sum = 0;
    assert (x.size() == y.size());
    for (int i = 0; i < x.size(); i++) sum += x[i] * y[i];
    return sum;
}

double vec_norm(vector<double>& x) { // L2 norm of vector x
    return sqrt(dot_prod(x, x));
}

void vec_minus(vector<double>& x, vector<double>& y) { // subtracts y from x in place, modifies x
    assert(x.size() == y.size());
    for (int i = 0; i < x.size(); i++) x[i] -= y[i];
}

void vec_add(vector<double>& x, vector<double>& y) { // adds y to x in place, modifies x
  assert(x.size() == y.size());
  for (int i = 0; i < x.size(); i++) x[i] += y[i];
}

void scalar_mult(vector<double>& x, double scalar) { // multiplies x by scalar in place, modifies x
    for (int i = 0; i < x.size(); i++) x[i] *= scalar;
}

vector<double> slice(vector<double>& x, int start, int stride, int length) { // returns a slice of length from x, starting at start, every stride elements
    vector<double> out(length);
    for (int i = 0; i < length; i++) out[i] = x[start + i * stride];
    return out;
}

vector<double> cross_product(vector<double>& x, vector<double>& y) { // cross product of x and y
    assert (x.size() >= 3);
    assert (y.size() >= 3);
    vector<double> vals(3);
    vals[0] = x[1] * y[2] - x[2] * y[1];
    vals[1] = x[2] * y[0] - x[0] * y[2];
    vals[2] = x[0] * y[1] - x[1] * y[0];
    return vals;
}

vector<double> cart_to_sphere(vector<double>& p1) { // turns cartesian coordinates to spherical coordinates
    assert (p1.size() >= 3);
    vector<double> sphere(3);
    double x = p1[0], y = p1[1], z = p1[2];
    sphere[0] = sqrt(x * x + y * y + z * z);
    sphere[1] = atan2(sqrt(x * x + y * y), z);
    sphere[2] = atan2(y, x);
    return sphere;
}

vector<double> sphere_to_cart(double radius, double colat, double lon) { // turns spherical coordinates to cartesian coordinates
    vector<double> cart(3);
    cart[0] = radius * sin(colat) * cos(lon);
    cart[1] = radius * sin(colat) * sin(lon);
    cart[2] = radius * cos(colat);
    return cart;
}

void project_to_sphere(vector<double>& p1, double radius) { // projects a point to the surface of a sphere of radius, modifies p1
    double norm = vec_norm(p1);
    for (int i = 0; i < p1.size(); i++) p1[i] *= radius / norm;
}

vector<double> project_to_sphere_2(vector<double> p1, double radius) { // projects a point to the surface of a sphere of radius, returns the new coordinates
    double norm = vec_norm(p1);
    vector<double> newcords(p1.size());
    for (int i = 0; i < p1.size(); i++) newcords[i] = radius * p1[i] / norm;
    return newcords;
}

double sphere_tri_area(vector<double>& p1, vector<double>& p2, vector<double>& p3, double radius) { // finds the area of a spherical triangle
    assert (p1.size() == p2.size());
    assert (p1.size() == p3.size());
    vector<double> p1n, p2n, p3n;
    double a, b, c, s, z, area;
    p1n = p1; // don't modify p1, modify p1n
    p2n = p2;
    p3n = p3;
    vec_minus(p2n, p3); // p2n = p2 - p3
    vec_minus(p3n, p1); // p3n = p3 - p1
    vec_minus(p1n, p2); // p1n = p1- p2
    // double a = acos(1 - 0.5 * vec_norm(p2n));
    // double b = acos(1 - 0.5 * vec_norm(p3n));
    // double c = acos(1 - 0.5 * vec_norm(p1n));
    a = acos(1 - 0.5 * dot_prod(p2n, p2n));
    b = acos(1 - 0.5 * dot_prod(p3n, p3n));
    c = acos(1 - 0.5 * dot_prod(p1n, p1n));
    s = (a + b + c) / 2;
    z = tan(s / 2) * tan((s - a) / 2) * tan((s - b) / 2) * tan((s - c) / 2);
    area = 4 * radius * radius * atan(sqrt(z));
    return area;
}

vector<double> lat_lon(vector<double>& p1) { // finds latitude and longitude of cartesian point
    assert (p1.size() >= 3);
    vector<double> latlon(2);
    double x = p1[0], y = p1[1], z = p1[2];
    latlon[0] = M_PI / 2 - atan2(sqrt(x * x + y * y), z); // latitude
    latlon[1] = atan2(y, x); // longitude
    return latlon;
}

vector<double> three_three_solve(vector<vector<double>> Amat, vector<double> b) { // solve 3x3 linear system explicitly, not in use
    vector<double> solution (3, 0);
    double a11 = Amat[0][0], a12 = Amat[0][1], a13 = Amat[0][2], a21 = Amat[1][0], a22 = Amat[1][1], a23 = Amat[1][2], a31 = Amat[2][0], a32 = Amat[2][1], a33 = Amat[2][2];
    double b1 = b[0], b2 = b[1], b3 = b[2];
    double det = -a13 * a22 * a31 + a12 * a23 * a31 + a13 * a21 * a32 - a11 * a32 * a23 - a12 * a21 * a33 + a11 * a22 * a33;
    solution[0] = ((a22 * a33 - a23 * a32) * b1 + (a13 * a32 - a12 * a33) * b2 + (a12 * a23 - a13 * a22) * b3) / det;
    solution[1] = ((a23 * a31 - a21 * a33) * b1 + (a11 * a33 - a13 * a31) * b2 + (a13 * a21 - a11 * a23) * b3) / det;
    solution[2] = ((a21 * a32 - a22 * a31) * b1 + (a12 * a31 - a11 * a32) * b2 + (a11 * a22 - a12 * a21) * b3) / det;
    return solution;
}

vector<double> barycoords(vector<double>& p1, vector<double>& p2, vector<double>& p3, vector<double>& p) { // finds triangle barycentric coordinates of point p
    assert (p1.size() == 3);
    assert (p2.size() == 3);
    assert (p3.size() == 3);
    assert (p.size() == 3);
    vector<double> coords(p.begin(), p.end());
    int dim = p.size();
    int nrhs = 1;
    vector<double> mat {p1[0], p1[1], p1[2], p2[0], p2[1], p2[2], p3[0], p3[1], p3[2]};
    vector<int> ipiv(3);
    int info;
    dgesv_(&dim, &nrhs, &*mat.begin(), &dim, &*ipiv.begin(), &*coords.begin(), &dim, &info);
    return coords;
}

vector<double> norm_barycoords(vector<double>& p1, vector<double>& p2, vector<double>& p3, vector<double>& p) { // returns normalized barycentric coordinates so b1 + b2 + b3 = 1
    vector<double> coords;
    coords = barycoords(p1, p2, p3, p);
    scalar_mult(coords, 1.0 / (coords[0] + coords[1] + coords[2]));
    return coords;
}

bool check_in_tri(vector<double>& p1, vector<double>& p2, vector<double>& p3, vector<double>& p) { // checks if point p is in triangle
    vector<double> bary_coord = barycoords(p1, p2, p3, p);
    if ((bary_coord[0] >= 0) and (bary_coord[1] >= 0) and (bary_coord[2] >= 0)) return true;
    else return false;
}

vector<double> tri_sides(vector<double>& p1, vector<double>& p2, vector<double>& p3) {// finds side lengths of triangle
    vector<double> sides(3);
    vector<double> p23 = p2;
    vector<double> p31 = p3;
    vector<double> p12 = p1;
    vec_minus(p23, p3);
    vec_minus(p31, p1);
    vec_minus(p12, p2);
    sides[0] = vec_norm(p23);
    sides[1] = vec_norm(p31);
    sides[2] = vec_norm(p12);
    return sides;
}

double circum_poly(double a, double b, double c) { // polynomial for circumcenter
    return a * a + b * b - c * c;
}

vector<double> tri_center(vector<double>& p1, vector<double>& p2, vector<double>& p3, double radius) { // finds triangle circumcenter x y z coords
    vector<double> tri_cent(3);
    vector<double> trisides = tri_sides(p1, p2, p3);
    double a = trisides[0], b = trisides[1], c = trisides[2];
    double scaling = (a + b + c) * (a + b - c) * (a - b + c) * (b + c - a);
    double coeff1 = a * a * circum_poly(b, c, a) / scaling;
    double coeff2 = b * b * circum_poly(c, a, b) / scaling;
    double coeff3 = c * c * circum_poly(a, b, c) / scaling;
    tri_cent[0] = coeff1 * p1[0] + coeff2 * p2[0] + coeff3 * p3[0];
    tri_cent[1] = coeff1 * p1[1] + coeff2 * p2[1] + coeff3 * p3[1];
    tri_cent[2] = coeff1 * p1[2] + coeff2 * p2[2] + coeff3 * p3[2];
    project_to_sphere(tri_cent, radius);
    return tri_cent;
}

double tri_radius(vector<double>& p1, vector<double>& p2, vector<double>& p3, vector<double>& center) { // finds triangle circumradius
    vector<double> p1c = p1;
    vector<double> p2c = p2;
    vector<double> p3c = p3;
    vec_minus(p1c, center);
    vec_minus(p2c, center);
    vec_minus(p3c, center);
    double radius1 = vec_norm(p1c);
    double radius2 = vec_norm(p2c);
    double radius3 = vec_norm(p3c);
    return max(max(radius1, radius2), radius3);
}

double great_circ_dist(vector<double>& p1, vector<double>& p2, double radius) { // finds great circle distance between two points
    double s = dot_prod(p1, p2) / (vec_norm(p1) * vec_norm(p2));
    double theta = acos(min(max(s, -1.0), 1.0));
    return theta * radius;
}

vector<double> BVE_gfunc(vector<double>& x, vector<double>& y) { // interaction function for barotropic vorticity equations
    double denom = 1 - dot_prod(x, y);
    vector<double> cross_prod = cross_product(x, y);
    scalar_mult(cross_prod, 1 / denom);
    return cross_prod;
}

int check_in_vec(vector<vector<double>>& x, vector<double>& y) { // checks if length 3 vector y is in vector of vectors x
    for (int i = 0; i < x.size(); i++) {
        if ((x[i][0] == y[0]) and (x[i][1] == y[1]) and (x[i][2] == y[2])) return i; // index where y is in x
    }
    return -1; // -1 if y not in x
}

int check_in_vec2(vector<double>& x, vector<double>& y, int max_points) { // checks if length 3 vector y is in state vector x
    // cout << "check in vec 2" << endl;
    // cout << x.size() << " " << 5 * max_points << endl;
    // cout << "y size " << y.size() << endl;
    // cout << x[5 * (max_points - 1) + 2] << " " << y[2] << endl;
    for (int i = 0; i < max_points; i++) {
        // cout << 5 * i + 2 << endl;
        // if ((x[5 * i] == y[0]) and (x[5 * i + 1] == y[1]) and (x[5 * i + 2] == y[2])) return i; // index where y is in x
        if ((abs(x[5 * i] - y[0]) < 1.0 / max_points) and (abs(x[5 * i + 1] - y[1]) < 1.0 / max_points) and (abs(x[5 * i + 2] - y[2]) < 1.0 / max_points)) return i;
    }
    return -1; // -1 if y not in x
}

int check_in_vec3(vector<vector<int>>& parent_verts, int v1, int v2) {
    for (int i = 0; i < parent_verts.size(); i++) {
        if ((parent_verts[i][0] == v1) and (parent_verts[i][1] == v2)) return i;
    }
    return -1;
}

void tri_interp(int iv1, int iv2, int iv3, vector<double>& v1, vector<double>& v2, vector<double>& v3, vector<double>& curr_state, vector<double>& target_points, vector<double>& curr_target, int i, double omega) {
    vector<double> bary;
    double absv1, absv2, absv3, absv;
    bary = norm_barycoords(v1, v2, v3, curr_target);
    // interpolate relative vorticity directly; results are decent
    // target_points[5 * i + 3] = bary[0] * curr_state[5 * iv1 + 3] + bary[1] * curr_state[5 * iv2 + 3] + bary[2] * curr_state[5 * iv3 + 3];
    // interpolate absolute vorticity; very similar results to interpolating relative vorticity directly, still diffusive
    absv1 = curr_state[5 * iv1 + 3] + 2 * omega * curr_state[5 * iv1 + 2];
    absv2 = curr_state[5 * iv2 + 3] + 2 * omega * curr_state[5 * iv2 + 2];
    absv3 = curr_state[5 * iv3 + 3] + 2 * omega * curr_state[5 * iv3 + 2];
    absv = bary[0] * absv1 + bary[1] * absv2 + bary[2] * absv3;
    target_points[5 * i + 3] = absv - 2 * omega * curr_target[2];
    // interpolate passive tracer
    target_points[5 * i + 4] = bary[0] * curr_state[5 * iv1 + 4] + bary[1] * curr_state[5 * iv2 + 4] + bary[2] * curr_state[5 * iv3 + 4];
}

vector<int> big_tri_verts(vector<vector<int>>& triangles, vector<vector<int>>& vert_tris, int iv1, int iv2, int iv3) {
    vector<int> points, iv1tris, iv2tris, iv3tris, v1tris, v2tris, v3tris;
    iv1tris = vert_tris[iv1];
    iv2tris = vert_tris[iv2];
    iv3tris = vert_tris[iv3];
    return points;
}

// int find_overlap(vector<int>& array1, vector<int>& array2, vector<int>& array3) { // finds common element of all 3 arrays
//     int common;
//     for (int i = 0; i < array1.size(); i++) {
//         common = array1[i];
//         for (int j = 0; j < array2.size(); j++) {
//             if (common == array2[j]) {
//                 for (int k = 0; k < array3.size(); k++) {
//                     if (common == array3[k]) return common; // return the common element
//                 }
//             }
//         }
//     }
//     return -1; // something is wrong
// }

void replace(vector<int>& vals, int find, int replacement) { // replace find in vals with replacement
    for (int i = 0; i < vals.size(); i++) {
        if (vals[i] == find) {
            vals[i] = replacement;
            break;
        }
    }
}

void read_points(int point_count, int tri_count, vector<double>& curr_state, vector<vector<int>>& vert_tris, vector<vector<int>>& triangles) {
    ifstream file1("../" + to_string(point_count) + "points_rh4.csv"); // ifstream = input file stream
    ifstream file2("../" + to_string(point_count) + "tris.csv");
    ifstream file3("../" + to_string(point_count) + "vert_tris.csv");
    ifstream file4("../" + to_string(point_count) + "vert_tri_count.csv");
    string line, word;
    int tri_counts;

    for (int i = 0; i < point_count; i++) {
        getline(file1, line);
        stringstream str1(line);
        for (int j = 0; j < 4; j++) { // read in initial condition of each point
            getline(str1, word, ',');
            curr_state[5 * i + j] = stod(word);
        }

        getline(file4, line);
        stringstream str4(line);
        getline(str4, word, ',');
        tri_counts = stod(word); // number of triangles each point borders

        vert_tris[i] = vector<int> (tri_counts);
        getline(file3, line);
        stringstream str3(line);
        for (int j = 0; j < tri_counts; j++) { // reads in each points adjacent triangles
            getline(str3, word, ',');
            vert_tris[i][j] = stod(word);
        }
    }

    for (int i = 0; i < tri_count; i++) { // reads in triangle vertex information
        getline(file2, line);
        stringstream str2(line);
        for (int j = 0; j < 3; j++) {
            getline(str2, word, ',');
            triangles[i][j] = stod(word);
        }
    }

    file1.close(); // close all the files we read from
    file2.close();
    file3.close();
    file4.close();
}

void area_init(vector<double>& curr_state, vector<double>& area, vector<vector<int>>& triangles, int tri_count) {
  int iv1, iv2, iv3;
  double curr_area;
  vector<double> v1, v2, v3;
  for (int i = 0; i < tri_count; i++) {
      // cout << i << endl;
      iv1 = triangles[i][0];
      iv2 = triangles[i][1];
      iv3 = triangles[i][2];
      // cout << iv1 << " " << iv2 << " " << iv3 << endl;
      v1 = slice(curr_state, 5 * iv1, 1, 3);
      v2 = slice(curr_state, 5 * iv2, 1, 3);
      v3 = slice(curr_state, 5 * iv3, 1, 3);
      curr_area = sphere_tri_area(v1, v2, v3, 1);
      area[iv1] += curr_area / 3.0;
      area[iv2] += curr_area / 3.0;
      area[iv3] += curr_area / 3.0;
      // cout << iv1 << endl;
  }
}

#endif
