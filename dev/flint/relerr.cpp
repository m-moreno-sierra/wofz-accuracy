#include <assert.h>
#include <iostream>
#include <stdio.h>
#include <flint/arb.h>
#include <flint/acb.h>
#include <flint/acb_hypgeom.h>
#include <cerf.h>
#include "relerr.h"
#include "wn.h"

std::complex<double> old_w_of_z(std::complex<double> z);

void my_acb_wofz(acb_t W, const acb_t Z, slong prec);

double relerr_wofz(bool old, double x, double y)
{
    acb_t Z;  acb_init(Z);
    arb_t E;  arb_init(E);
    acb_t W;  acb_init(W);
    acb_t W2; acb_init(W2);
    arb_t A1; arb_init(A1);
    arb_t A2; arb_init(A2);

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

    const std::complex<double> w2 = old ?
	old_w_of_z({x, y}) :
	w_of_z(std::complex<double>{x, y});
    arb_set_d(acb_realref(W2), w2.real());
    arb_set_d(acb_imagref(W2), w2.imag());

    acb_sub(W2, W, W2, prec);
    acb_abs(A2, W2, prec);
    acb_abs(A1, W, prec);
    arb_div(A1, A2, A1, prec);

    double ret = arf_get_d(arb_midref(A1), ARF_RND_NEAR);

    acb_clear(Z);
    arb_clear(E);
    acb_clear(W);
    acb_clear(W2);
    arb_clear(A1);
    arb_clear(A2);

    return ret;
}
