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
//  License:   MIT License (file LICENSE in toplevel directory)
//  Copyright: Forschungszentrum Jülich GmbH 2025
//  Author:    Joachim Wuttke <j.wuttke@fz-juelich.de>
//
//  ************************************************************************************************

#include "ref.h"
#include <cassert>
#include <flint/arb.h>
#include "hp_fref.h"

std::complex<double> Ref::fref(double x, double y)
{
    acb_t Z;  acb_init(Z);
    arb_t E;  arb_init(E);
    acb_t W;  acb_init(W);

    arb_set_d(acb_realref(Z), x);
    arb_set_d(acb_imagref(Z), y);
    int prec = 103;
    HP::acb_fref(W, Z, prec);

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

std::vector<Ref::Coeff> Ref::fn_vector(double x, double y, int N)
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
    HP::acb_fcoeffs(V, Z, prec);

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

	ret[n] = {{wx, wy}, awn, rre};
    }

    acb_clear(Z);
    acb_clear(F);
    arb_clear(E);
    for (acb_t& v : V)
	acb_clear(v);

    return ret;
}

// Returns list of indices of lattice points within domain, for given inverse lattice constant inv_b.
std::vector<std::pair<int,int>> Ref::domain_grid(double inv_b)
{
    std::vector<std::pair<int,int>> ret;
    const double R = 7;
    for (int ix = 0;; ++ix) {
	double x = ix / inv_b;
	if (x>R)
	    break;
	for (int iy = 0;; ++iy) {
	    double y = iy / inv_b;
	    if (x*x + y*y >= R*R)
		break;
	    ret.emplace_back(std::pair<int,int>{ix, iy});
	}
    }
    return ret;
}
