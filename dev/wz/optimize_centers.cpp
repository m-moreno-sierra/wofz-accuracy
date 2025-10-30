//  ************************************************************************************************
//
//  Project:   Code generating code for Taylor expansion on square tiles
//
//  Reference: Joachim Wuttke,
//             Code generation for computing an analytical function with near machine precision
//             on square tiles, with application to the Faddeeva function
//             (in preparation, available upon request)
//
//  File:      optimize_centers.cpp
//
//  Purpose:   Program optimize_centers is Algorithm 2 [Acenters] of the reference paper.
//             It determines near-optimum expansion centers, and writes a data file
//             that is an essential input to the subsequent set cover computation.
//
//  License:   MIT License (file LICENSE in toplevel directory)
//  Copyright: Forschungszentrum Jülich GmbH 2025
//  Author:    Joachim Wuttke <j.wuttke@fz-juelich.de>
//
//  ************************************************************************************************

#include <cassert>
#include <chrono>
#include <cstdio>
#include <cstdlib>
#include <format>
#include <iostream>
#include <limits>
#include <set>
#include <omp.h>
#include "terms.h"

#define SQR(x) ((x)*(x))

namespace {

//! Returns list of squared polyomino circumcircle diameters for given parities pi_x, pi_y.
//! This is Algorithm 1 of the reference paper.
std::vector<int> sorted_diameters(int pi_x, int pi_y, int kappamax)
{
    std::set<int> set;
    set.insert(0);
    for (int p=1; p<sqrt(kappamax); ++p) {
	for (int q=1;; ++q) {
	    int kappa = SQR(2*p-pi_x) + SQR(2*q-pi_y);
	    if (kappa > kappamax)
		break;
	    set.insert(kappa);
	}
    }
    return std::vector<int>(set.begin(), set.end());
}

int iSref(bool pi_x, bool pi_y)
{
    if (pi_x && pi_y)
	return 0;
    else if ((!pi_x) && (!pi_y))
	return 2;
    else
	return 1;
}

double double_neighbor(double x, int n)
{
    double ret = x;
    for (int i=0; i<std::abs(n); ++i)
	ret = nextafter(ret, n>0 ? std::numeric_limits<double>::infinity() :
			-std::numeric_limits<double>::infinity());
    return ret;
}

std::vector<std::vector<double>> wOnGrid(double inv_b, int kappamax)
{
    std::vector<std::vector<double>> ret;
    const double R=7;
    const double jrmax = R*(inv_b/2) + sqrt(kappamax)/2;

    for (int jx = 0; jx <= jrmax; ++jx) {
	std::vector<double> wm;
	for (int jy = 0; jx*jx+jy*jy <= jrmax*jrmax; ++jy)
	    wm.emplace_back(abs(Ref::fref(jx/(inv_b/2), jy/(inv_b/2))));
	ret.emplace_back(wm);
    }
    return ret;
}

double wmin(int ix, int iy, int kappa, const std::vector<std::vector<double>>& W)
{
    double ret = std::numeric_limits<double>::infinity();
    for (int dix = ix%2; dix < sqrt(kappa); dix += 2) {
        int diy = iy%2 + int((sqrt(kappa-dix*dix)-iy%2)/2);
        ret = std::min(ret, W.at((ix+dix)/2).at((iy+diy)/2));
    }
    return ret;
}

//! Full characterization of one Taylor expansion center.
//! Such centers are computed with respect to a given b-grid,
//! i.e. a square lattice with lattice constant b, which is
//! one half of the lattice constant a of the given square tiling.
struct OneCenter {
    int ix;                 //!< x index of nearby b-grid point.
    int iy;                 //!< y index of nearby b-grid point.
    int kappa;                 //!< squared covered polyomino circumcircle diameter in units of b^2.
    std::complex<double> z; //!< exact coordinates of expansion center.
};

//! Determines maximum kappa (with rho<=delta) for given expansion center c.
//! Return values:
//! - itau, index in Si such that kappa_max = Si[itau].
//! - maxerr, maximum error within kappa, in units of delta.
std::tuple<int, double> max_kappa(int N, double delta, double inv_b, const OneCenter& c,
                               const std::vector<int>& Si, const std::vector<std::vector<double>>& WW)
{
    const int nS = Si.size();
    const int nSlb = int(log2(nS));
    int itau = 0;
    double maxerr = std::numeric_limits<double>::infinity();
    std::vector<Ref::Coeff> WN = Ref::fn_vector(c.z);
    for (int n = nSlb; n>=0; --n) {
        int itmp = itau + (1 << n);
        if (itmp >= nS)
            continue;
        const int kappa = Si[itmp];
        const double tau = sqrt(kappa) / inv_b;
        const double te = Terms::truncation_error(c.z, N, tau, WN);
        if (std::isinf(te))
            continue;
        const double re = Terms::rounding_error(c.z, N, tau, WN);
        const double wmi = ::wmin(c.ix, c.iy, kappa, WW);
        double err = (te+re) / wmi;
        if (err <= delta) {
            itau = itmp;
            maxerr = err;
        }
    }
    if (itau >= nS-1)
        throw std::runtime_error("increase kappamax!");
    assert(maxerr <= delta || itau == 0);
    return {itau, maxerr};
}

} // namespace


//  Purpose:   Program optimize_centers is Algorithm 2 [Acenters] of the reference paper.
//             It determines near-optimum expansion centers, and writes a data file
//             that is an essential input to the subsequent set cover computation.
int main(int argc, char *argv[])
{
    if (argc != 5) {
        fprintf(stderr, "For each grid point, determine optimized expansion center and range.\n");
        fprintf(stderr, "Usage: %s N_Taylor delta inverse_a M_offsets\n", argv[0]);
        fprintf(stderr, "where\n");
        fprintf(stderr, "  N_Taylor  = order of expansion\n");
        fprintf(stderr, "  delta     = maximum error in units of epsilon\n");
        fprintf(stderr, "  inverse_a = 1/a where a is lattice constant of square tiling\n");
        fprintf(stderr, "  M_offsets = M : try neighbor grid points from -M to M as expansion centers\n");
        fprintf(stderr, "Results are written to stdout.\n");
        fprintf(stderr, "Redirect them into a file for use in subsequent set-cover computation\n");
        return 1;
    }
    char *endptr;
    const int N = strtol(argv[1], &endptr, 10);
    assert(*endptr == '\0');
    const double delta = strtod(argv[2], &endptr);
    assert(*endptr == '\0');
    // Grid of expansion centers has half the lattice constant of the give square tiling:
    const double inv_b = 2 * strtod(argv[3], &endptr);
    assert(*endptr == '\0');
    const int M = strtol(argv[4], &endptr, 10);
    assert(*endptr == '\0');

    const auto now = std::chrono::system_clock::now();
    printf("# Created by %s at %s\n", argv[0], std::format("{:%d-%m-%Y %H:%M:%OS}", now).c_str());
    printf("# Given parameters:\n");
    printf("#   N_Taylor = %i (order of expansion)\n", N);
    printf("#   delta = %g (maximum error in units of epsilon)\n", delta);
    printf("#   M_recenter = %i (potential expansion centers from range -M...M)\n", M);
    printf("#   1/a = %g (inverse lattice constant of square tiling)\n", inv_b / 2);
    printf("# Output format:\n");
    printf("#   Block header line entry:\n");
    printf("#     x of grid point\n");
    printf("#   Block lines entries:\n");
    printf("#     y of grid point\n");
    printf("#     kappa, squared covered polyomino circumcircle diameter in units of (a/2)^2\n");
    printf("#     x' of expansion center (close to x)\n");
    printf("#     y' of expansion center (close to y)\n");

    const int kappamax = 1023;

    // kappa sequences start with 0. Otherwise as in paper.
    const std::vector<std::vector<int>> S {
	::sorted_diameters(0, 0, kappamax),
	::sorted_diameters(0, 1, kappamax),
	::sorted_diameters(1, 1, kappamax) };
    assert(S[0][4]==40);
    assert(S[1][4]==25);
    assert(S[2][4]==26);

    const std::vector<std::vector<double>> WW = ::wOnGrid(inv_b, kappamax);

    std::vector<std::pair<int,int>> I = Ref::domain_grid(inv_b); // List of grid points within domain.

    std::vector<::OneCenter> VR(I.size());

    #pragma omp parallel for
    for (size_t i=0; i<VR.size(); ++i) {
	const auto [ix, iy] = I[i];
	const double x = ix / inv_b;
	const double y = iy / inv_b;
	const std::vector<int>& Si = S[::iSref(ix%2, iy%2)];

        OneCenter c = {ix, iy, 0, {x, y}}; // will be optimized in-place

	// Determine maximum kappa (with rho<=delta) for lattice point z.
        auto [itau, maxerr] = ::max_kappa(N, delta, inv_b, c, Si, WW);

	// Determine nearby zd with minimum rho at kappa determined above.
	const int kappa = Si[itau];
	const double tau = sqrt(kappa) / inv_b;
	const double wmi = ::wmin(ix, iy, kappa, WW);
	for (int idx=-M; idx<=M; ++idx) {
	    if (ix==0 && idx!=0)
		continue;
	    const double xd = ::double_neighbor(x, idx);
	    for (int idy=-M; idy<=M; ++idy) {
		if (idx==0 and idy==0)
		    continue;
		if (iy==0 && idy!=0)
		    continue;
		const double yd = ::double_neighbor(y, idy);
		const std::complex<double> zd{xd, yd};
		const std::vector<Ref::Coeff> WN = Ref::fn_vector(zd);
		const double te = Terms::truncation_error(zd, N, tau, WN);
		if (std::isinf(te))
		    continue;
		const double re = Terms::rounding_error(zd, N, tau, WN);
		const double err = (te+re) / wmi;
		if (err < maxerr) {
		    maxerr = err;
		    c.z = zd;
		}
	    }
	}

	// Determine maximum kappa (with rho<=delta) for z determined above.
        auto [itau2, dummy] = ::max_kappa(N, delta, inv_b, c, Si, WW);
        c.kappa = Si[itau2];
	VR[i] = c;
    }

    // Print list of results, divided in blocks with same x.
    int ixold = -1;
    for (const auto& r : VR) {
	const auto [ix, iy, kappa, z] = r;
	if (ix != ixold) {
	    if (ix>0)
		std::cout << "\n";
	    std::cout << std::format("{:f}", ix / inv_b) << "\n";
	}
	std::cout << std::format("{:f} {:d} 0x{:a} 0x{:a}\n", iy / inv_b, kappa, z.real(), z.imag());
	ixold = ix;
    }
    std::cout << "\n";
    return 0;
}
