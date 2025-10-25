//  ************************************************************************************************
//
//  Code generating code for piecewise Chebyshev approximation
//
//  Reference: Joachim Wuttke and Alexander Kleinsorge,
//             Algorithm 1XXX: Code generation for piecewise Chebyshev approximation
//
//  File:      f_imwofx.cpp
//
//  Arbitrary-precision implementation of function Im w(x), for use from PPAPP.
//
//  License:   GNU General Public License, version 3 or higher (see src/LICENSE)
//  Copyright: Forschungszentrum Jülich GmbH 2025
//  Author:    Joachim Wuttke <j.wuttke@fz-juelich.de>
//
//  ************************************************************************************************

#include <flint/arb_hypgeom.h>
#include "f.h"

//! The total intermediate domain [a,b) to be covered by piecewise polynomial approximation:

const std::pair<double, double> my_domain = { 0.5, 12. };

//! Evaluates the overflow-compensated imaginary error function
//!   f(x) = exp(-x^2) * erfi(x) = Im w(x),
//! using 'prec' binary digits.
//! Sets F to f(x) with x given by X.

void my_arb_f(arb_t& F, const arb_t& X, slong prec)
{
    arb_t T1; arb_init(T1);
    arb_t T2; arb_init(T2);

    arb_hypgeom_erfi(T1, X, prec);
    arb_sqr(T2, X, prec);
    arb_neg(T2, T2);
    arb_exp(T2, T2, prec); // exp(-x^2)
    arb_mul(F, T1, T2, prec);

    arb_clear(T1);
    arb_clear(T2);
}

//! A few test cases, just to secure against gross errors.

const std::vector<std::array<double,3>> my_testcases = {
    {.5, .478925, 1e-5},
    {1, .607158, 1e-5},
    {12, .0471808, 1e-5}
};
