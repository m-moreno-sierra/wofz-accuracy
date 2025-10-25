//  ************************************************************************************************
//
//  High-precision computation of the half width of the Voigt function, for use from PPAPP.
//
//  File:      devtool/ppapp/f_hwhm_voigt.cpp
//
//  See also:  https://jugit.fz-juelich.de/mlz/ppapp
//
//  License:   GNU General Public License, version 3 or higher (see src/LICENSE)
//  Copyright: Forschungszentrum Jülich GmbH 2025
//  Authors:   Eva Mukherjee, Joachim Wuttke <j.wuttke@fz-juelich.de>
//
//  ************************************************************************************************

#include "f.h"
#include <array>
#include <flint/acb.h>
#include <flint/acb_hypgeom.h>
#include <flint/arb.h>
#include <iostream>
#include <vector>

//! The total intermediate domain [a,b) to be covered by piecewise polynomial approximation:

const std::pair<double, double> my_domain = { 0.25, 16. };

// Computing the Faddeeva Function
void acb_faddeeva_w(acb_t res, const acb_t z, slong prec) {
    acb_t iz, minus_iz, erfc_val, exp_val;
    acb_init(iz);
    acb_init(minus_iz);
    acb_init(erfc_val);
    acb_init(exp_val);
    // iz = -Im(z) + i Re(z)
    arb_t re, im;
    arb_init(re);
    arb_init(im);
    acb_get_real(re, z);
    acb_get_imag(im, z);
    arb_neg(acb_realref(iz), im); // real part = -Im(z)
    arb_set(acb_imagref(iz), re); // imag part = Re(z)

    arb_clear(re);
    arb_clear(im);

    acb_neg(minus_iz, iz);                      // minus_iz = -iz
    acb_hypgeom_erfc(erfc_val, minus_iz, prec); // erfc_val = erfc(minus_iz)
    acb_t z_squared, neg_z_squared;             // exp_val = exp(-z^2)
    acb_init(z_squared);
    acb_init(neg_z_squared);
    acb_mul(z_squared, z, z, prec);        // z^2
    acb_neg(neg_z_squared, z_squared);     // -z^2
    acb_exp(exp_val, neg_z_squared, prec); // exp(-z^2)

    acb_clear(z_squared);
    acb_clear(neg_z_squared);
    acb_mul(res, exp_val, erfc_val, prec); // res = exp(-z^2) * erfc(-i * z)
    acb_clear(iz);
    acb_clear(minus_iz);
    acb_clear(erfc_val);
    acb_clear(exp_val);
}


void my_arb_voigt(arb_t &F, const arb_t &q, const arb_t &sigma,
                  const arb_t &gamma, slong prec) {
    arb_t T1, T2;
    arb_init(T1);
    arb_init(T2);
    acb_t z_complex_arg, minus_iz_complex_arg, w_z_result;
    acb_init(z_complex_arg);
    acb_init(minus_iz_complex_arg);
    acb_init(w_z_result);

    // Calculating the complex argument z = q / (sqrt(2) * sigma) + i * gamma /
    // (sqrt(2) * sigma)
    arb_set_d(T1, 2.0);
    arb_sqrt(T1, T1, prec);       // T1 = sqrt(2)
    arb_mul(T2, T1, sigma, prec); // T2 = sqrt(2) * sigma
    arb_div(T1, q, T2, prec);     // T1 = q / (sqrt(2)*sigma) (real part)
    arb_div(T2, gamma, T2, prec); // T2 = gamma / (sqrt(2)*sigma) (imag part)
    acb_set_arb_arb(z_complex_arg, T1, T2); // z = T1 + i*T2

    // Calculating -i * z for the Faddeeva function (w(z) = erfcx(-iz))
    // If z = x + iy, then -iz = -i(x+iy) = -ix - i^2y = y - ix
    // So, the real part of (-iz) is Im(z), and the imaginary part of (-iz) is
    // -Re(z).
    arb_t real_part, imag_part;
    arb_init(real_part);
    arb_init(imag_part);
    acb_get_imag(imag_part, z_complex_arg);
    acb_get_real(real_part, z_complex_arg);
    arb_neg(real_part, real_part);
    acb_set_arb_arb(minus_iz_complex_arg, imag_part, real_part);
    // Faddeeva function
    acb_faddeeva_w(w_z_result, z_complex_arg, prec);
    acb_get_real(T1, w_z_result);

    // Calculating 1 / (sqrt(2*pi)*sigma)
    arb_const_pi(T2, prec);
    arb_mul_2exp_si(T2, T2, 1); // 2*pi
    arb_sqrt(T2, T2, prec);     // sqrt(2*pi)
    arb_mul(T2, T2, sigma, prec);
    arb_inv(T2, T2, prec);
    // Computing the Voigt Function
    arb_mul(F, T2, T1, prec);

    arb_clear(T1);
    arb_clear(T2);
    arb_clear(real_part);
    arb_clear(imag_part);
    acb_clear(z_complex_arg);
    acb_clear(minus_iz_complex_arg);
    acb_clear(w_z_result);
}

void my_arb_f(arb_t &F, const arb_t &X, slong prec) {
    arb_t sigma, gamma;
    arb_init(sigma);
    arb_init(gamma);

    // The input 'X' represents the ratio gamma/sigma.
    arb_set_d(sigma, 1.0);
    arb_set(gamma, X);

    // Calculate the target value: V(0; sigma, gamma) / 2.
    arb_t target_V;
    arb_init(target_V);
    arb_t q_zero;
    arb_init(q_zero);
    arb_set_d(q_zero, 0.0);
    my_arb_voigt(target_V, q_zero, sigma, gamma, prec);
    arb_mul_2exp_si(target_V, target_V, -1);
    arb_clear(q_zero);

    // Bisection Method
    // Finding the root of the function f(q) = V(q; sigma, gamma) - target value
    arb_t a, b, c, val_c, diff, tol;
    arb_init(a);
    arb_init(b);
    arb_init(c);
    arb_init(val_c);
    arb_init(diff);
    arb_init(tol);
    arb_zero(a);
    arb_set_d(b, 1.0);

    // Finding b such that V(b) < target_V
    for (int iter = 0; iter < 1000; iter++) {
        my_arb_voigt(val_c, b, sigma, gamma, prec);
        if (arb_lt(val_c, target_V))
            goto found_b;
        arb_mul_2exp_si(b, b, 1);
    }
    std::cerr << "Loop to find b failed" << std::endl;
    std::exit(1);
found_b:

    // Bisection loop
    slong max_iter = 20 * prec;
    for (slong i = 0; i < max_iter; i++) {
        arb_sub(diff, b, a, prec);
        arb_one(tol);
        arb_mul_2exp_si(tol, tol, -prec);
        if (arb_le(diff, tol))
            goto found_x;
        arb_add(c, a, b, prec);
        arb_mul_2exp_si(c, c, -1);
        my_arb_voigt(val_c, c, sigma, gamma, prec);
        if (arb_gt(val_c, target_V)) {
            arb_set(a, c); // V(c) too high => mid is too small
        } else {
            arb_set(b, c); // V(c) too low => mid is too large
        }
    }
found_x:

    arb_sub(diff, b, a, prec);
    arb_one(tol);
    arb_mul_2exp_si(tol, tol, -prec);
    if (arb_gt(diff, tol)) {
        std::cerr << "Warning: bisection may not have converged for X=";
        arb_fprintd(stderr, X, 10);
        std::cerr << "\n";
        std::exit(1);
    }
    // The result F is the HWHM which is midpoint of [a, b]
    arb_add(F, a, b, prec);
    arb_mul_2exp_si(F, F, -1);

    arb_clear(sigma);
    arb_clear(gamma);
    arb_clear(q_zero);
    arb_clear(target_V);
    arb_clear(a);
    arb_clear(b);
    arb_clear(c);
    arb_clear(val_c);
    arb_clear(diff);
    arb_clear(tol);
}

//! Function for HWHM calculation.
void hwhm_direct_calculation(arb_t &HWHM, const arb_t &sigma, const arb_t &gamma, slong prec) {
    // Calculating V(0; sigma, gamma)
    arb_t V_zero;
    arb_init(V_zero);
    arb_t q_zero;
    arb_init(q_zero);
    arb_set_d(q_zero, 0.0);
    my_arb_voigt(V_zero, q_zero, sigma, gamma, prec);
    arb_clear(q_zero);

    // Calculating the target value: 0.5 * V(0; sigma, gamma)
    arb_t target_V;
    arb_init(target_V);
    arb_mul_2exp_si(target_V, V_zero, -1);

    // Bisection Method
    // Finding the root of the function f(q) = V(q; sigma, gamma) - target value
    arb_t a, b, c, val_c, diff, tol;
    arb_init(a);
    arb_init(b);
    arb_init(c);
    arb_init(val_c);
    arb_init(diff);
    arb_init(tol);
    arb_zero(a);
    arb_set_d(b, 1.0);

    // Finding b such that V(b) < target_V
    for (int iter = 0; iter < 10; iter++) {
        my_arb_voigt(val_c, b, sigma, gamma, prec);
        if (arb_lt(val_c, target_V))
            break;
        arb_mul_2exp_si(b, b, 1);
    }

    // Bisection loop
    slong max_iter = 2 * prec;
    for (slong i = 0; i < max_iter; i++) {
        arb_sub(diff, b, a, prec);
        arb_one(tol);
        arb_mul_2exp_si(tol, tol, -prec);
        if (arb_le(diff, tol))
            break;
        arb_add(c, a, b, prec);
        arb_mul_2exp_si(c, c, -1);
        my_arb_voigt(val_c, c, sigma, gamma, prec);
        if (arb_gt(val_c, target_V)) {
            arb_set(a, c); // V(c) too high => mid is too small
        } else {
            arb_set(b, c); // V(c) too low => mid is too large
        }
    }

    // The result HWHM is the midpoint of [a, b]
    arb_add(HWHM, a, b, prec);
    arb_mul_2exp_si(HWHM, HWHM, -1);

    arb_clear(V_zero);
    arb_clear(target_V);
    arb_clear(a);
    arb_clear(b);
    arb_clear(c);
    arb_clear(val_c);
    arb_clear(diff);
    arb_clear(tol);
}


extern const std::vector<std::array<double, 3>> my_testcases = {
    {1.0, 1.80057, 1e-5}, {2.0, 2.93434, 1e-5}, {0.5, 1.28521, 1e-5}};
