#ifndef interp_H
#define interp_H

#include <iostream>
#include <cmath>
#include <vector>
#include <queue>
#include <chrono>
#include <Accelerate/Accelerate.h>

struct interp_struct {
    int degree;
    int point_count;
    int info;
    vector<double> matrix;
    vector<vector<double>> interp_points;
    int ipiv[];
};

void __attribute__((optnone)) fekete_init(vector<vector<double>>& points, int degree)  { // initializes fekete matrix, local
// void __attribute__((optimize(0))) fekete_init(vector<vector<double>>& points, int degree)  { // initializes fekete matrix, on GL
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

interp_struct interp_init(int degree) {
    interp_struct interp1;
    interp1.degree = degree;
    interp1.point_count = (degree + 1) * (degree + 2) / 2;
    interp1.interp_points = vector<vector<double>> (interp1.point_count, vector<double> (3, 0));
    interp1.matrix = vector<double> (interp1.point_count * interp1.point_count, 0);
    fekete_init(interp1.interp_points, degree);
    interp_mat_init(interp1.matrix, interp1.interp_points, degree, interp1.point_count);
    dgetrf_(&interp1.point_count, &interp1.point_count, &*interp1.matrix.begin(), &interp1.point_count, interp1.ipiv, &interp1.info);
    if (interp1.info > 0) {
        cout << "dgetrf: " << interp1.info << endl;
    }
    return interp1;
}

#endif
