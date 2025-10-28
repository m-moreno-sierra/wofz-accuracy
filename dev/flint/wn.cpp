//  ************************************************************************************************
//
//  Project:   Code generating code for Taylor expansion on square tiles
//
//  Reference: Joachim Wuttke,
//             Code generation for computing an analytical function with near machine precision
//             on square tiles, with application to the Faddeeva function
//             (in preparation, available upon request)
//
//  File:      wn.cpp
//
//  Purpose:   Provides applications-specific high-precision functions wofz, w_n_vector.
//
//  License:   GNU General Public License, version 3 or higher (see src/LICENSE)
//  Copyright: Forschungszentrum Jülich GmbH 2025
//  Author:    Joachim Wuttke <j.wuttke@fz-juelich.de>
//
//  ************************************************************************************************

#include <assert.h>
#include <iostream>
#include <stdio.h>
#include <flint/arb.h>
#include <flint/acb.h>
#include <flint/acb_hypgeom.h>
#include "wn.h"

void my_acb_wofz(acb_t W, const acb_t Z, slong prec)
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
    acb_mul(W, T1, T2, prec);

    acb_clear(T0);
    acb_clear(T1);
    acb_clear(T2);
    acb_clear(I);
}

void my_acb_w_coeffs(std::vector<acb_t>& V, const acb_t Z, slong prec)
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
    my_acb_wofz(V[0], Z, prec);

    for(int n = 1; n<V.size(); ++n) {
	const acb_t& T2 = n==1 ? T1 : V[n-2];
	acb_mul(V[n], V[n-1], Z, prec);
	acb_add(V[n], V[n], n==1 ? T1 : V[n-2], prec);
	acb_mul_si(V[n], V[n], -2, prec);
	acb_div_ui(V[n], V[n], n, prec);
    }

    arb_clear(X);
    arb_clear(Y);
    acb_clear(T1);
}

std::complex<double> wofz(double x, double y)
{
    acb_t Z;  acb_init(Z);
    arb_t E;  arb_init(E);
    acb_t W;  acb_init(W);

    arb_set_d(acb_realref(Z), x);
    arb_set_d(acb_imagref(Z), y);
    int prec = 103;
    my_acb_wofz(W, Z, prec);

    acb_get_real(E, W);
    double wx = arf_get_d(arb_midref(E), ARF_RND_NEAR);
    int xbits = arb_rel_accuracy_bits(E);
    acb_get_imag(E, W);
    double wy = arf_get_d(arb_midref(E), ARF_RND_NEAR);
    int ybits = arb_rel_accuracy_bits(E);

    assert(wx==0 || xbits==-1 || xbits>53);
    assert(wy==0 || xbits==-1 || xbits>53);

    acb_clear(Z);
    arb_clear(E);
    acb_clear(W);

    return {wx, wy};
}

std::vector<Coeff> w_n_vector(double x, double y, int N)
{
    acb_t Z;  acb_init(Z);
    acb_t F;  acb_init(F);
    arb_t E;  arb_init(E);

    std::vector<Coeff> ret(N);
    std::vector<acb_t> V(N);
    for (acb_t& v : V)
	acb_init(v);

    arb_set_d(acb_realref(Z), x);
    arb_set_d(acb_imagref(Z), y);
    const int prec = 343;
    my_acb_w_coeffs(V, Z, prec);
    for (int n = 0; n<N; ++n) {
	// Extract real and imag parts of V[n] and their accuracies.
	acb_get_real(E, V[n]);
	const double wx = arf_get_d(arb_midref(E), ARF_RND_NEAR);
	const int xbits = arb_rel_accuracy_bits(E);
	acb_get_imag(E, V[n]);
	const double wy = arf_get_d(arb_midref(E), ARF_RND_NEAR);
	const int ybits = arb_rel_accuracy_bits(E);
	assert(wx==0 || xbits==-1 || xbits>53);
	assert(wy==0 || xbits==-1 || xbits>53);

	// Absolute value of w_n:
	const double awn = hypot(wx,wy);

	// Relative rounding error:
	arb_set_d(acb_realref(F), wx);
	arb_set_d(acb_imagref(F), wy);
	acb_sub(F, F, V[n], 53);
	acb_abs(E, F, 53);
	const double rre = arf_get_d(arb_midref(E), ARF_RND_NEAR) / awn / pow(2, -53);

	//printf("w_%i(z) = ", n);
	//acb_printd(V[n], 12); printf(" [%i,%i bits]\n", xbits, ybits);
	ret[n] = {{wx, wy}, awn, rre};
    }

    acb_clear(Z);
    acb_clear(F);
    arb_clear(E);
    for (acb_t& v : V)
	acb_clear(v);

    return ret;
}
