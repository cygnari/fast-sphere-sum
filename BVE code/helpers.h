// helper functions

#ifndef F_H
#define F_H

#include <iostream>
#include <cmath>
#include <array>
#include <vector>
// #include <Accelerate/Accelerate.h>
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
    latlon[0] = M_PI / 2 - atan2(sqrt(x * x + y * y), z);
    latlon[1] = atan2(y, x);
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

void icos_init(vector<vector<double>>& verts, vector<vector<vector<double>>>& tri_info, vector<vector<vector<int>>>& tri_verts, double radius, int levels) { // initializes icosahedron for fast summation
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
}

void points_assign(vector<vector<vector<int>>>& tri_verts, vector<vector<double>>& verts, vector<double>& points, vector<vector<vector<int>>>& tri_points, vector<vector<int>>& point_locs, int levels, int point_count) { // finds which icosahedron triangle each dynamics point is in
    vector<double> v1, v2, v3, point;
    int iv1, iv2, iv3, lb, ub;
    point_locs.clear();
    tri_points.clear();
    point_locs.push_back(vector<int> (point_count, 0));
    tri_points[0] = vector<vector<int>> (20);
    for (int i = 0; i < point_count; i++) { // finds the placement of each point in the initial icosahedron
        point = slice(points, 5 * i, 1, 3);
        for (int j = 0; j < 20; j++) {
            iv1 = tri_verts[0][j][0];
            iv2 = tri_verts[0][j][1];
            iv3 = tri_verts[0][j][2];
            v1 = verts[iv1];
            v2 = verts[iv2];
            v3 = verts[iv3];
            if (check_in_tri(v1, v2, v3, point)) {
                point_locs[0][i] = j;
                tri_points[0][j].push_back(i);
                break;
            }
        }
    }
    for (int i = 1; i < levels; i++) { // finds point loc in refined faces
        point_locs.push_back(vector<int> (point_count, 0));
        tri_points[i] = vector<vector<int>> (20 * pow(4, i));
        for (int j = 0; j < point_count; j++) {
            point = slice(points, 5 * j, 1, 3);
            lb = 4 * point_locs[i-1][j]; // utilize tree structure to minimize searching
            ub = lb + 4;
            for (int k = lb; k < ub; k++) {
                iv1 = tri_verts[i][k][0];
                iv2 = tri_verts[i][k][1];
                iv3 = tri_verts[i][k][2];
                v1 = verts[iv1];
                v2 = verts[iv2];
                v3 = verts[iv3];
                if (check_in_tri(v1, v2, v3, point)) {
                    point_locs[i][j] = k;
                    tri_points[i][k].push_back(j);
                    break;
                }
            }
        }
    }
}

double interp_eval(vector<double>& alphas, double s, double t, int degree) { // evaluate interpolation polynomial with coefficients alpha and barycentric point (s, t)
    // return alphas[0] + alphas[1] * s + alphas[2] * t + alphas[3] * s * t + alphas[4] * s * s + alphas[5] * t * t;
    double accum = 0;
    int index;
    for (int i = 0; i < degree + 1; i++) {
        for (int j = 0; j < i + 1; j++) {
            index = i * (i + 1) / 2 + j;
            accum += pow(s, i - j) * pow(t, j) * alphas[index = i * (i + 1) / 2 + j];
        }
    }
    return accum;
}

// void __attribute__((optnone)) fekete_init(vector<vector<double>>& points, int degree)  { // initializes fekete matrix, local
void __attribute__((optimize(0))) fekete_init(vector<vector<double>>& points, int degree)  { // initializes fekete matrix, on GL
    double delta_x = 1.0 / degree;
    int index;
    double a, b, c, d;
    for (int i = 0; i < degree + 1; i++) {
        a = 1 - i * delta_x;
        for (int j = 0; j < i + 1; j++) {
            index = i * (i + 1) / 2 + j;
            c = j * delta_x;
            b = 1 - a - b;
            a = 0.5 * (1 + sin(M_PI / 2 * (2 * a - 1)));
            b = 0.5 * (1 + sin(M_PI / 2 * (2 * b - 1)));
            c = 0.5 * (1 + sin(M_PI / 2 * (2 * c - 1)));
            points[index][0] = a / (a + b + c);
            points[index][1] = b / (a + b + c);
            points[index][2] = c / (a + b + c);
        }
    }
}

void interp_mat_init(vector<double>& mat, vector<vector<double>>& points, int degree, int point_count) { // sets up matrix to interpolate with fekete points
    int index, place;
    double a, b;
    for (int i = 0; i < degree + 1; i++) {
        for (int j = 0; j < i + 1; j++) {
            index = i * (i + 1) / 2 + j;
            for (int k = 0; k < point_count; k++) {
                a = points[k][0];
                b = points[k][1];
                place = index * point_count + k;
                mat[place] = pow(a, i - j) * pow(b, j);
            }
        }
    }
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

void regrid_points(vector<double>& curr_state, vector<double>& target_points, vector<vector<int>>& triangles, vector<vector<int>>& vert_tris, int point_count, int tri_count, double omega, int lb, int ub, int ID) { // remesh back to original particle locations
    vector<double> curr_target, curr_pos, v1, v2, v3; // need to fix vorticity remeshing
    vector<int> poss_tris;
    // double curr_vor;
    int test_count, iv1, iv2, iv3;
    int success = 0, failure = 0, bad = 0;
    bool found;
    for (int i = lb; i < ub; i++) {
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
    // cout << "success: " << success << " failure: " << failure << " bad: " << bad << endl;
    // cout << "success: " << success << " failure: " << failure << " bad: " << bad << " ID " << ID << endl;
}

void regrid_point2(vector<double>& curr_state, vector<double>& target_points, vector<vector<int>>& triangles, vector<vector<int>>& vert_tris, int point_count, int tri_count, double omega, int ID, vector<int> particle_thread) { // remesh back to original particle locations
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

vector<int> amr(vector<double>& curr_state, vector<vector<int>>& triangles, vector<vector<int>>& vert_tris, vector<double>& areas, vector<vector<int>>& parent_verts, int tri_count, int point_count, int max_points) { // adaptive mesh refinement
    double circulation_tol = 0.0025 / 4.0; // for 2562 particles, reduce for more starting particles
    double vorticity_tol = 0.2 / 2.0; // for 2562 particles, reduce for more starting particles
    int iv1, iv2, iv3, iv12, iv23, iv31,  itriv1v12v31, itriv2v12v23, itriv3v31v23, itriv12v23v31, old_tri_count;
    double vor1, vor2, vor3, max_val, min_val, vor1n, vor2n, vor3n, circulation, tri_area;
    vector<double> v1, v2, v3, v12, v23, v31, p1, p2, p3;
    vector<int> return_values {tri_count, point_count}; // return new tri_count and point_count;
    vector<int> iv1tris, iv2tris, iv3tris;
    if (point_count >= max_points) return return_values;
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
    return_values = {tri_count, point_count};
    return return_values;
}

#endif
