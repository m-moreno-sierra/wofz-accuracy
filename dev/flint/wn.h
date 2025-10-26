#ifndef CERF_FLINT_WN_H
#define CERF_FLINT_WN_H

#include <complex>
#include <vector>

struct Coeff {
    std::complex<double> fn;
    double abs;
    double rel_rounding_err;
};

std::complex<double> wofz(double x, double y);

std::vector<Coeff> w_n_vector(double x, double y, int N=53);

#endif // CERF_FLINT_WN_H
