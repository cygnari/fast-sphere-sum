// helper functions

#ifndef F_omp
#define F_omp

#include <iostream>
#include <cmath>
#include <array>
#include <vector>
// #include <Accelerate/Accelerate.h>
#include <cassert>

using namespace std;

extern "C" {
    extern int dgesv_(int*,int*,double*,int*,int*,double*,int*,int*);
    extern int dgetrf_(int*,int*,double*,int*,int*,int*);
    extern int dgetrs_(char*,int*,int*,double*,int*,int*,double*,int*,int*);
}

double dot_prod(vector<double>& x, vector<double>& y) {
    double sum = 0;
    assert (x.size() == y.size());
    for (int i = 0; i < x.size(); i++) sum += x[i] * y[i];
    return sum;
}

double vec_norm(vector<double>& x) {
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

void vec_add_omp(vector<double>& x, vector<double>& y, int lb, int ub) { // adds y to x in place, modifies x
  assert(x.size() == y.size());
  for (int i = lb; i < ub; i++) x[i] += y[i];
}

void scalar_mult(vector<double>& x, double scalar) { // multiplies x by scalar in place, modifies x
    for (int i = 0; i < x.size(); i++) x[i] *= scalar;
}

void scalar_mult_omp(vector<double>& x, double scalar, int lb, int ub) { // multiplies x by scalar in place, modifies x
    for (int i = lb; i < ub; i++) x[i] *= scalar;
    // for (int i = lb; i < ub; i++) {
    //     x[i] *= scalar;
    //     cout << "i " << i << endl;
    // }
}

void vec_add_scalar_mult_omp(vector<double>& x, vector<double>& y, double scalar, int lb, int ub) { // x = x + scalar * y
    for (int i = lb; i < ub; i++) x[i] += y[i] * scalar;
}

void vec_scalar_mult_add_omp(vector<double>& x, vector<double>& y, double scalar, int lb, int ub) { // x = scalar * x + y
    for (int i = lb; i < ub; i++) x[i] = y[i] + scalar * x[i];
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
    // vector<double> newcords(p1.size());
    for (int i = 0; i < p1.size(); i++) p1[i] *= radius / norm;
    // return newcords;
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
    vector<double> p1n = p1; // don't modify p1, modify p1n
    vector<double> p2n = p2;
    vector<double> p3n = p3;
    vec_minus(p2n, p3); // p2n = p2 - p3
    vec_minus(p3n, p1); // p3n = p3 - p1
    vec_minus(p1n, p2); // p1n = p1- p2
    double a = acos(1 - 0.5 * vec_norm(p2n));
    double b = acos(1 - 0.5 * vec_norm(p3n));
    double c = acos(1 - 0.5 * vec_norm(p1n));
    double s = (a + b + c) / 2;
    double z = tan(s / 2) * tan((s - a) / 2) * tan((s - b) / 2) * tan((s - c) / 2);
    double area = 4 * radius * radius * atan(sqrt(z));
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

vector<double> three_three_solve(vector<vector<double>> Amat, vector<double> b) {
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
    // vector<vector<double>> amatrix {{p1[0], p2[0], p3[0]}, {p1[1], p2[1], p3[1]}, {p1[2], p2[2], p3[2]}};
    vector<int> ipiv(3);
    int info;
    dgesv_(&dim, &nrhs, &*mat.begin(), &dim, &*ipiv.begin(), &*coords.begin(), &dim, &info);
    // return three_three_solve(amatrix, p);
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

double circum_poly(double a, double b, double c) {
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

vector<double> BVE_gfunc(vector<double>& x, vector<double>& y) {
    double denom = 1 - dot_prod(x, y);
    vector<double> cross_prod = cross_product(x, y);
    scalar_mult(cross_prod, 1 / denom);
    return cross_prod;
}

int check_in_vec(vector<vector<double>>& x, vector<double>& y) {
    for (int i = 0; i < x.size(); i++) {
        if ((x[i][0] == y[0]) and (x[i][1] == y[1]) and (x[i][2] == y[2])) return i; // index where y is in x
    }
    return -1; // -1 if y not in x
}

void icos_init(vector<vector<double>>& verts, vector<vector<vector<double>>>& tri_info, vector<vector<vector<int>>>& tri_verts, double radius, int levels) {
    double phi = (1 + sqrt(5)) / 2;
    verts.push_back(project_to_sphere_2(vector<double> {0, 1, phi}, radius));
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
    tri_verts[0][1].insert(tri_verts[0][1].end(), {1, 2, 11});
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
    for (int i = 0; i < 20; i++) {
        for (int j = 0; j < 3; j++) tri_verts[0][i][j] -= 1;
        int iv1 = tri_verts[0][i][0];
        int iv2 = tri_verts[0][i][1];
        int iv3 = tri_verts[0][i][2];
        vector<double> v1 = verts[iv1];
        vector<double> v2 = verts[iv2];
        vector<double> v3 = verts[iv3];
        vector<double> center = tri_center(v1, v2, v3, radius);
        tri_info[0][i].insert(tri_info[0][i].end(), center.begin(), center.end()); // index 0 1 2 is the triangle center
        tri_info[0][i].push_back(tri_radius(v1, v2, v3, center)); // index 3 is the triangle radius
    }
    for (int i = 0; i < levels; i++) {
        tri_info.push_back(vector<vector<double>> (20 * pow(4, i + 1), vector<double> (0)));
        tri_verts.push_back(vector<vector<int>> (20 * pow(4, i + 1), vector<int> (0)));
        for (int j = 0; j < 20 * pow(4, i); j++) {
            int iv1 = tri_verts[i][j][0];
            int iv2 = tri_verts[i][j][1];
            int iv3 = tri_verts[i][j][2];
            vector<double> v1 = verts[iv1];
            vector<double> v2 = verts[iv2];
            vector<double> v3 = verts[iv3];
            vector<double> v12 = v1;
            vector<double> v23 = v2;
            vector<double> v31 = v3;
            vec_add(v12, v2);
            vec_add(v23, v3);
            vec_add(v31, v1);
            scalar_mult(v12, 0.5);
            scalar_mult(v23, 0.5);
            scalar_mult(v31, 0.5);
            project_to_sphere(v12, radius);
            project_to_sphere(v23, radius);
            project_to_sphere(v31, radius);
            int iv12 = check_in_vec(verts, v12);
            int iv13 = check_in_vec(verts, v31);
            int iv23 = check_in_vec(verts, v23);
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
            tri_verts[i+1][4*j].insert(tri_verts[i+1][4*j].end(), {iv1, iv13, iv12});
            tri_verts[i+1][4*j+1].insert(tri_verts[i+1][4*j+1].end(),{iv3, iv23, iv13});
            tri_verts[i+1][4*j+2].insert(tri_verts[i+1][4*j+2].end(),{iv2, iv12, iv23});
            tri_verts[i+1][4*j+3].insert(tri_verts[i+1][4*j+3].end(),{iv12, iv13, iv23});

            vector<double> center = tri_center(v1, v12, v31, radius);
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

void points_assign(vector<vector<vector<int>>>& tri_verts, vector<vector<double>>& verts, vector<double>& points, vector<vector<vector<int>>>& tri_points, vector<vector<int>>& point_locs, int levels, int point_count, int lb_p, int ub_p, omp_lock_t writelock) {
    vector<double> v1;
    vector<double> v2;
    vector<double> v3;
    vector<double> point;
    int iv1, iv2, iv3;
    point_locs.clear();
    tri_points.clear();
    point_locs.push_back(vector<int> (point_count, 0));
    // tri_points.push_back(vector<vector<int>> (20, vector<int> (0)));
    tri_points[0] = vector<vector<int>> (20);
    for (int i = lb_p; i < ub_p; i++) {
        point = slice(points, 4 * i, 1, 3);
        // iv1 = tri_verts[0][point_locs[0][i]][0];
        // iv2 = tri_verts[0][point_locs[0][i]][1];
        // iv3 = tri_verts[0][point_locs[0][i]][2];
        // v1 = verts[iv1];
        // v2 = verts[iv2];
        // v3 = verts[iv3];
        // if (check_in_tri(v1, v2, v3, point)) {
        //     continue;
        // } else {
        //     for (int j = 0; j < 20; j++) {
        //         iv1 = tri_verts[0][j][0];
        //         iv2 = tri_verts[0][j][1];
        //         iv3 = tri_verts[0][j][2];
        //         v1 = verts[iv1];
        //         v2 = verts[iv2];
        //         v3 = verts[iv3];
        //         if (check_in_tri(v1, v2, v3, point)) {
        //             point_locs[0][i] = j;
        //             tri_points[0][j].push_back(i);
        //             cout << "point i " << i << " tri j " << j << endl;
        //             break;
        //         }
        //     }
        // }
        for (int j = 0; j < 20; j++) {
            iv1 = tri_verts[0][j][0];
            iv2 = tri_verts[0][j][1];
            iv3 = tri_verts[0][j][2];
            v1 = verts[iv1];
            v2 = verts[iv2];
            v3 = verts[iv3];
            if (check_in_tri(v1, v2, v3, point)) {
                point_locs[0][i] = j;
// #pragma omp critical
                // omp_set_lock(&writelock);
                tri_points[0][j].push_back(i);
                // omp_unset_lock(&writelock);
                // cout << "point i " << i << " tri j " << j << endl;
                break;
            }
        }
    }
    for (int i = 1; i < levels; i++) {
        point_locs.push_back(vector<int> (point_count, 0));
        // tri_points.push_back(vector<vector<int>> (20 * pow(4, i), vector<int> (0)));
        tri_points[i] = vector<vector<int>> (20 * pow(4, i));
        for (int j = lb_p; j < ub_p; j++) {
            point = slice(points, 4 * j, 1, 3);
            int lb = 4 * point_locs[i-1][j];
            int ub = lb + 4;
            for (int k = lb; k < ub; k++) {
                iv1 = tri_verts[i][k][0];
                iv2 = tri_verts[i][k][1];
                iv3 = tri_verts[i][k][2];
                v1 = verts[iv1];
                v2 = verts[iv2];
                v3 = verts[iv3];
                if (check_in_tri(v1, v2, v3, point)) {
                    point_locs[i][j] = k;
// #pragma omp critical
                    // omp_set_lock(&writelock);
                    tri_points[i][k].push_back(j);
                    // omp_unset_lock(&writelock);
                    break;
                }
            }
        }
    }
}

double interp_eval(vector<double>& alphas, double s, double t) {
    return alphas[0] + alphas[1] * s + alphas[2] * t + alphas[3] * s * t + alphas[4] * s * s + alphas[5] * t * t;
}

void zero_out(vector<double>& target, int lb, int ub) { // zeros out target vector
    for (int i = lb; i < ub; i++) target[i] = 0;
}

#endif
