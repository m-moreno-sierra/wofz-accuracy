//  ************************************************************************************************
//
//  Project:   Code generating code for Taylor expansion on square tiles
//
//  Reference: Joachim Wuttke,
//             Code generation for computing an analytical function with near machine precision
//             on square tiles, with application to the Faddeeva function
//             (in preparation, available upon request)
//
//  File:      wn.h
//
//  Purpose:   Declares applications-specific high-precision functions wofz, w_n_vector.
//
//  License:   MIT License (file LICENSE in toplevel directory)
//  Copyright: Forschungszentrum Jülich GmbH 2025
//  Author:    Joachim Wuttke <j.wuttke@fz-juelich.de>
//
//  ************************************************************************************************

#ifndef CERF_FLINT_WN_H
#define CERF_FLINT_WN_H

#include <complex>
#include <vector>

//! The reference function in standard double precision.

namespace Ref {

//! A complex Taylor coefficient, its absolute value and its relative rounding error.
struct Coeff {
    std::complex<double> fn;
    double abs;
    double rel_rounding_err;
};

std::complex<double> fref(double x, double y);

std::vector<Coeff> fn_vector(std::complex<double> z, int N = 53);

// Returns list of indices of lattice points within domain, for given inverse lattice constant
// inv_b.
std::vector<std::pair<int, int>> domain_grid(double inv_b);

} // namespace Ref

#endif // CERF_FLINT_WN_H
