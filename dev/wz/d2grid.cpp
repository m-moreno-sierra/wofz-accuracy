//  ************************************************************************************************
//
//  Project:   Code generating code for Taylor expansion on square tiles
//
//  Reference: Joachim Wuttke,
//             Code generation for computing an analytical function with near machine precision
//             on square tiles, with application to the Faddeeva function
//             (in preparation, available upon request)
//
//  File:      d2grid.cpp
//
//  Purpose:   Provides program d2grid.
//
//  License:   GNU General Public License, version 3 or higher (see src/LICENSE)
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
#include "terms.h"
#include <omp.h>

#define SQR(x) ((x)*(x))

std::vector<int> sorted_diameters(int sx, int sy, int d2max)
{
    std::set<int> set;
    set.insert(0);
    for (int p=1; p<sqrt(d2max); ++p) {
	for (int q=1;; ++q) {
	    int d2 = SQR(2*p-sx) + SQR(2*q-sy);
	    if (d2 > d2max)
		break;
	    set.insert(d2);
	}
    }
    return std::vector<int>(set.begin(), set.end());
}

int iSref(bool sx, bool sy)
{
    if (sx && sy)
	return 0;
    else if ((!sx) && (!sy))
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

int main(int argc, char *argv[])
{
    if (argc != 5) {
        fprintf(stderr, "Usage: %s N_Taylor delta inverse_a M_offsets\n", argv[0]);
        return 1;
    }
    char *endptr;
    const int N = strtol(argv[1], &endptr, 10);
    assert(*endptr == '\0');
    const double delta = strtod(argv[2], &endptr);
    assert(*endptr == '\0');
    const double inv_a = strtod(argv[3], &endptr);
    assert(*endptr == '\0');
    const int M = strtol(argv[4], &endptr, 10);
    assert(*endptr == '\0');

    const auto now = std::chrono::system_clock::now();
    printf("# Created by %s at %s\n", argv[0], std::format("{:%d-%m-%Y %H:%M:%OS}", now).c_str());
    printf("# N_Taylor = %i\n", N);
    printf("# delta = %g\n", delta);
    printf("# M_recenter = %i\n", M);
    printf("# 1/a = %g\n", inv_a);

    const double R = 7;
    const double inv_b = 2 * inv_a;
    const int d2max = 1023;

    // D2 sequences start with 0. Otherwise as in paper.
    const std::vector<std::vector<int>> S {
	sorted_diameters(0, 0, d2max),
	sorted_diameters(0, 1, d2max),
	sorted_diameters(1, 1, d2max) };
    assert(S[0][4]==40);
    assert(S[1][4]==25);
    assert(S[2][4]==26);

    const std::vector<std::vector<double>> WW = wOnGrid(inv_b, d2max);

    std::vector<std::pair<int,int>> I;
    for (int ix = 0;; ++ix) {
	double x = ix / inv_b;
	if (x>R)
	    break;
	std::string out = std::format("{:f}\n", x);
	for (int iy = 0;; ++iy) {
	    double y = iy / inv_b;
	    if (x*x + y*y >= R*R)
		break;
	    I.emplace_back(std::pair<int,int>{ix, iy});
	}
    }

    struct Result {
	int ix;
	int iy;
	int d2;
	double x;
	double y;
    };
    std::vector<Result> VR(I.size());

    #pragma omp parallel for
    for (int i=0; i<VR.size(); ++i) {
	const auto [ix, iy] = I[i];
	double x = ix / inv_b;
	double y = iy / inv_b;
	std::complex<double> z{x, y};
	const std::vector<int>& Si = S[iSref(ix%2, iy%2)];
	const int nS = Si.size();
	const int nSlb = int(log2(nS));

	// Determine maximum d2 (with rho<=delta) for lattice point z.
	int itau = 0;
	int imin = 0;
	double minerr = std::numeric_limits<double>::infinity();
	std::vector<Coeff> WN = w_n_vector(x, y);
	for (int n = nSlb; n>=0; --n) {
	    int itmp = itau + (1 << n);
	    if (itmp >= nS)
		continue;
	    const int d2 = Si[itmp];
	    const double tau = sqrt(d2) / (2 * inv_a);
	    const double te = truncation_error(z, N, tau, WN);
	    if (std::isinf(te))
		continue;
	    const double re = rounding_error(z, N, tau, WN);
	    const double wmi = wmin(ix, iy, d2, WW);
	    double err = (te+re) / wmi;
	    if (err <= delta) {
		itau = itmp;
		minerr = err;
	    }
	}
	if (itau >= nS-1)
	    throw std::runtime_error("increase d2max!");
	assert(minerr <= delta || itau == 0);

	// Determine nearby zd with minimum rho at d2 determined above.
	const int d2 = Si[itau];
	const double tau = sqrt(d2) / (2 * inv_a);
	const double wmi = wmin(ix, iy, d2, WW);
	for (int idx=-M; idx<=M; ++idx) {
	    if (ix==0 && idx!=0)
		continue;
	    double xd = double_neighbor(x, idx);
	    for (int idy=-M; idy<=M; ++idy) {
		if (idx==0 and idy==0)
		    continue;
		if (iy==0 && idy!=0)
		    continue;
		double yd = double_neighbor(y, idy);
		const std::complex<double> zd{xd, yd};
		WN = w_n_vector(xd, yd);
		const double te = truncation_error(zd, N, tau, WN);
		if (std::isinf(te))
		    continue;
		const double re = rounding_error(zd, N, tau, WN);
		double err = (te+re) / wmi;
		if (err < minerr) {
		    minerr = err;
		    z = zd;
		}
	    }
	}
	x = z.real();
	y = z.imag();

	// Determine maximum d2 (with rho<=delta) for z determined above.
	itau = 0;
	WN = w_n_vector(x, y);
	for (int n = nSlb; n>=0; --n) {
	    int itmp = itau + (1 << n);
	    if (itmp >= nS)
		continue;
	    const int d2 = Si[itmp];
	    const double tau = sqrt(d2) / (2 * inv_a);
	    const double te = truncation_error(z, N, tau, WN);
	    if (std::isinf(te))
		continue;
	    const double re = rounding_error(z, N, tau, WN);
	    const double wmi = wmin(ix, iy, d2, WW);
	    const double err = (te+re) / wmi;
	    if (err <= delta) {
		itau = itmp;
	    }
	}
	if (itau >= nS-1)
	    throw std::runtime_error("increase d2max!");

	VR[i] = Result{ix, iy, Si[itau], x, y};
    }

    int ixold = -1;
    for (const auto& r : VR) {
	const auto [ix, iy, d2, x, y] = r;
	if (ix != ixold) {
	    if (ix>0)
		std::cout << "\n";
	    std::cout << std::format("{:f}", ix / inv_b) << "\n";
	}
	std::cout << std::format("{:f} {:d} 0x{:a} 0x{:a}\n", iy / inv_b, d2, x, y);
	ixold = ix;
    }
    std::cout << "\n";
    return 0;
}
