#ifndef interp_H
#define interp_H

#include <iostream>
#include <cmath>
#include <vector>
#include <queue>
#include <chrono>
#include <Accelerate/Accelerate.h>
#include "struct_list.h"

void __attribute__((optnone)) fekete_init(vector<vector<double>>& points, int degree)  { // initializes fekete matrix, local
// void __attribute__((optimize(0))) fekete_init(vector<vector<double>>& points, int degree)  { // initializes fekete matrix, on GL
    double delta_x = 1.0 / degree;
    int index;
    double a, b, c;
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

void interp_init(run_information& config1) {
    // interp_struct interp1;
    // config1.degree = degree;
    config1.interp_count = (config1.degree + 1) * (config1.degree + 2) / 2;
    config1.interp_points = vector<vector<double>> (config1.interp_count, vector<double> (3, 0));
    config1.interp_matrix = vector<double> (config1.interp_count * config1.interp_count, 0);
    config1.ipiv = vector<int> (config1.interp_count, 0);
    fekete_init(config1.interp_points, config1.degree);
    interp_mat_init(config1.interp_matrix, config1.interp_points, config1.degree, config1.interp_count);
    dgetrf_(&config1.interp_count, &config1.interp_count, &*config1.interp_matrix.begin(), &config1.interp_count, &*config1.ipiv.begin(), &config1.info);
    if (config1.info > 0) {
        cout << "dgetrf: " << config1.info << endl;
    }
    // return interp1;
}

#endif
