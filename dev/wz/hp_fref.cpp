//  ************************************************************************************************
//
//  Project:   Code generating code for Taylor expansion on square tiles
//
//  Reference: Joachim Wuttke,
//             Code generation for computing an analytical function with near machine precision
//             on square tiles, with application to the Faddeeva function
//             (in preparation, available upon request)
//
//  File:      hp_fref.cpp
//
//  Purpose:   Provides functions in namespace HP that compute the target function f_ref(z)
//             and its Taylor coefficients with arbitrary precision.
//
//  License:   MIT License (file LICENSE in toplevel directory)
//  Copyright: Forschungszentrum Jülich GmbH 2025
//  Author:    Joachim Wuttke <j.wuttke@fz-juelich.de>
//
//  ************************************************************************************************

#include <cassert>
#include <flint/arb.h>
#include <flint/acb_hypgeom.h>
#include "hp_fref.h"

//! Computes 'F' := f_ref('Z') with a precision of 'prec' decimal digits.
//! Implementation for f_ref = w, Faddeeva function.

void HP::acb_fref(acb_t F, const acb_t Z, slong prec)
{
    acb_t T0; acb_init(T0);
    acb_t T1; acb_init(T1);
    acb_t T2; acb_init(T2);
    acb_t I;  acb_init(I);  acb_onei(I);

    acb_mul(T0, Z, I, prec);
    acb_neg(T0, T0);
    acb_hypgeom_erfc(T1, T0, prec); // erfc(-iz)
    acb_sqr(T0, T0, prec);
    acb_exp(T2, T0, prec); // exp(-z^2)=exp((-iz)^2)
    acb_mul(F, T1, T2, prec);

    acb_clear(T0);
    acb_clear(T1);
    acb_clear(T2);
    acb_clear(I);
}

//! Computes list of Taylor coefficient 'V', such that V[n] = (1/n!)*d^n f_ref(z)/dz^n
//! with a precision of 'prec' decimal digits.
//! Array V must be initialized with required size.
//! Implementation for f_ref = w, Faddeeva function, using the recursion DLMF 7.10.3.

void HP::acb_fcoeffs(std::vector<acb_t>& V, const acb_t Z, slong prec)
{
    acb_t T1; acb_init(T1);
    arb_t X; arb_init(X); arb_zero(X);
    arb_t Y; arb_init(Y);
    arb_const_pi(Y, prec);
    arb_sqrt(Y, Y, prec);
    arb_inv(Y, Y, prec);
    arb_neg(Y, Y); // -1/sqrt(pi)
    acb_set_arb_arb(T1, X, Y);

    assert(V.size() > 0);
    acb_fref(V[0], Z, prec);

    for(size_t n = 1; n<V.size(); ++n) {
	acb_mul(V[n], V[n-1], Z, prec);
	acb_add(V[n], V[n], n==1 ? T1 : V[n-2], prec);
	acb_mul_si(V[n], V[n], -2, prec);
	acb_div_ui(V[n], V[n], n, prec);
    }

    arb_clear(X);
    arb_clear(Y);
    acb_clear(T1);
}
