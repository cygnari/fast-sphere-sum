#ifndef gfunc_H
#define gfunc_H

#include "general_utils.hpp"

vector<double> BVE_gfunc(vector<double>& x, vector<double>& y) { // interaction function for barotropic vorticity equations
    double denom = 1.0 - dot_prod(x, y);
    vector<double> cross_prod = cross_product(x, y);
    scalar_mult(cross_prod, 1.0 / denom);
    scalar_mult(cross_prod, -1.0 / (4.0 * M_PI));
    return cross_prod;
}

double stream_gfunc(vector<double>& x, vector<double>& y) {
    double interior = 1.0 - dot_prod(x, y);
    return -1.0 / (4 * M_PI) * log(interior);
}

#endif
