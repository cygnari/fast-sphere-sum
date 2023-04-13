#ifndef F_H
#define F_H

#include <cmath>
#include <array>
#include <vector>
#include <Accelerate/Accelerate.h>
#include <cassert>
#include <algorithm>

extern "C" { // lapack
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

void scalar_mult(vector<double>& x, double scalar) { // multiplies x by scalar in place, modifies x
    for (int i = 0; i < x.size(); i++) x[i] *= scalar;
}

void vec_add(vector<double>& x, vector<double>& y) { // adds y to x in place, modifies x
    assert(x.size() == y.size());
    for (int i = 0; i < x.size(); i++) x[i] += y[i];
}

void vec_minus(vector<double>& x, vector<double>& y) { // subtracts y from x in place, modifies x
    assert(x.size() == y.size());
    for (int i = 0; i < x.size(); i++) x[i] -= y[i];
}

void scalar_mult2d(vector<vector<double>>& x, double scalar) {
    for (int i = 0; i < x.size(); i++) {
        for (int j = 0; j < x[i].size(); j++) x[i][j] *= scalar;
    }
}

void vec_add2d(vector<vector<double>>& x, vector<vector<double>>& y) {
    assert(x.size() == y.size());
    for (int i = 0; i < x.size(); i++) {
        assert(x[i].size() == y[i].size());
        for (int j = 0; j < x[i].size(); j++) x[i][j] += y[i][j];
    }
}

void vec_minus2d(vector<vector<double>>& x, vector<vector<double>>& y) {
    assert(x.size() == y.size());
    for (int i = 0; i < x.size(); i++) {
        assert(x[i].size() == y[i].size());
        for (int j = 0; j < x[i].size(); j++) x[i][j] -= y[i][j];
    }
}

vector<double> slice(vector<double>& x, int start, int stride, int length) {
    // returns a slice of length from x, starting at start, every stride elements
    vector<double> out(length);
    for (int i = 0; i < length; i++) out[i] = x[start + i * stride];
    return out;
}

void vector_copy(vector<double>& x, vector<double>& y, int start, int length) {
    // copies elements from y into x
    for (int i = 0; i < length; i++) x[start + i] = y[i];
}

void vector_copy2(vector<double>& x, vector<double> y, int start, int length) {
    // copies elements from y into x
    for (int i = 0; i < length; i++) x[start + i] = y[i];
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
    assert (p1.size() == 3);
    vector<double> sphere(3);
    double x = p1[0], y = p1[1], z = p1[2];
    sphere[0] = sqrt(x * x + y * y + z * z); // radius
    sphere[1] = atan2(sqrt(x * x + y * y), z); // colatitude
    sphere[2] = atan2(y, x); // longitude
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

vector<double> project_to_sphere_2(vector<double> p1, double radius) {
    // projects a point to the surface of a sphere of radius, returns the new coordinates
    double norm = vec_norm(p1);
    vector<double> newcords(p1.size());
    for (int i = 0; i < p1.size(); i++) newcords[i] = radius * p1[i] / norm;
    return newcords;
}

double great_circ_dist(vector<double>& p1, vector<double>& p2, double radius) {
    // finds great circle distance between two points
    double s = dot_prod(p1, p2) / (vec_norm(p1) * vec_norm(p2));
    double theta = acos(min(max(s, -1.0), 1.0));
    return theta * radius;
}

double sphere_tri_area(vector<double>& p1, vector<double>& p2, vector<double>& p3, double radius) {
    // finds the area of a spherical triangle
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

vector<double> barycoords(vector<double>& p1, vector<double>& p2, vector<double>& p3, vector<double>& p) {
    // finds triangle barycentric coordinates of point p
    assert (p1.size() == 3);
    assert (p2.size() == 3);
    assert (p3.size() == 3);
    assert (p.size() == 3);
    vector<double> coords(p.begin(), p.end());
    int dim = 3, nrhs = 1, info;
    vector<double> mat {p1[0], p1[1], p1[2], p2[0], p2[1], p2[2], p3[0], p3[1], p3[2]};
    vector<int> ipiv(3);
    dgesv_(&dim, &nrhs, &*mat.begin(), &dim, &*ipiv.begin(), &*coords.begin(), &dim, &info);
    return coords;
}

vector<double> normalized_barycoords(vector<double>& p1, vector<double>& p2, vector<double>& p3, vector<double>& p) {
    // returns normalized barycentric coordinates so b1 + b2 + b3 = 1
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

double circum_poly(double a, double b, double c) { // polynomial for circumcenter
    return a * a + b * b - c * c;
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

vector<double> circum_center(vector<double>& p1, vector<double>& p2, vector<double>& p3, double radius) {
    // finds triangle circumcenter x y z coords
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

vector<double> BVE_gfunc(vector<double>& x, vector<double>& y) { // interaction function for barotropic vorticity equations
    double denom = 1.0 - dot_prod(x, y);
    vector<double> cross_prod = cross_product(x, y);
    scalar_mult(cross_prod, 1.0 / denom);
    return cross_prod;
}

void replace(vector<int>& vals, int find, int replacement) { // replace find in vals with replacement
    for (int i = 0; i < vals.size(); i++) {
        if (vals[i] == find) {
            vals[i] = replacement;
            break;
        }
    }
}

void area_init(vector<double>& curr_state, vector<double>& area, vector<vector<int>>& triangles, int tri_count) {
    int iv1, iv2, iv3;
    double curr_area;
    vector<double> v1, v2, v3;
    for (int i = 0; i < tri_count; i++) {
        iv1 = triangles[i][0];
        iv2 = triangles[i][1];
        iv3 = triangles[i][2];
        v1 = slice(curr_state, 5 * iv1, 1, 3);
        v2 = slice(curr_state, 5 * iv2, 1, 3);
        v3 = slice(curr_state, 5 * iv3, 1, 3);
        curr_area = sphere_tri_area(v1, v2, v3, 1);
        area[iv1] += curr_area / 3.0;
        area[iv2] += curr_area / 3.0;
        area[iv3] += curr_area / 3.0;
    }
}

int check_point_exist(vector<vector<int>>& parent_points, int point_count, int iv1, int iv2) {
    for (int i = 0; i < point_count; i++) {
        if ((parent_points[i][0] == iv1) and (parent_points[i][1] == iv2)) return i;
    }
    return -1;
}

#endif
