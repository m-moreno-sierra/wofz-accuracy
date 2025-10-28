//  ************************************************************************************************
//
//  Project:   Code generating code for Taylor expansion on square tiles
//
//  Reference: Joachim Wuttke,
//             Code generation for computing an analytical function with near machine precision
//             on square tiles, with application to the Faddeeva function
//             (in preparation, available upon request)
//
//  File:      worst_per_tile.cpp
//
//  Purpose:   Provides program worst_per_tile.
//
//  License:   GNU General Public License, version 3 or higher (see src/LICENSE)
//  Copyright: Forschungszentrum Jülich GmbH 2025
//  Author:    Joachim Wuttke <j.wuttke@fz-juelich.de>
//
//  ************************************************************************************************

#include <cassert>
#include <iostream>
#include <fstream>
#include <cstdio>
#include <cstdlib>
#include <chrono>
#include <format>
#include <random>
#include <cerf.h>
#include <omp.h>
#include "terms.h"

int main(int argc, char *argv[])
{
    if (argc != 3) {
        fprintf(stderr, "Usage: %s inverse_a n_per_tile\n", argv[0]);
        return 1;
    }
    char *endptr;
    double inv_a = strtod(argv[1], &endptr);
    assert(*endptr == '\0');
    int npt = strtol(argv[2], &endptr, 10);
    assert(*endptr == '\0');

    const int ntot = int(10*inv_a);
    struct Result {
	double x; // nominal x of tile
	double y; // nominal y of tile
	double d; // worst relative error
	double worst_x;
	double worst_y;
    };
    std::vector<Result> result;
    result.reserve(ntot*ntot);
    for (int jx=0; jx<ntot; ++jx) {
	const double x = jx/inv_a;
	for (int jy=0; jy<ntot; ++jy) {
	    const double y = jy/inv_a;
	    result.push_back(Result{x, y, 0});
	}
    }

    std::random_device rdev;

    #pragma omp parallel for
    for (int j = 0; j < size(result); ++j) {
	Result& res = result[j];

	std::mt19937 reng(j);
	std::uniform_real_distribution<double> rdis(0, 1/inv_a);

	res.d = 0;
	for (int i=0; i<npt; ++i) {
	    double x = res.x + rdis(reng);
	    double y = res.y + rdis(reng);
	    const std::complex<double> wref = wofz(x, y);
	    const std::complex<double> wlib = w_of_z(std::complex<double>{x, y});
	    const double relerr = std::abs(wlib-wref) / std::abs(wref);
	    if (relerr > res.d) {
		res.d = relerr;
		res.worst_x = x;
		res.worst_y = y;
	    }
	}
    }

    auto now = std::time(nullptr);
    {
	std::string fname = "/tmp/worst_per_tile.tab";
        std::ofstream f(fname);
        if (!f)
	    throw std::runtime_error("Failed to open " + fname);
        f << "# Created by " << argv[0] << " inv_a=" << inv_a << "n/tile=" << npt
          << " on " << std::put_time(std::localtime(&now), "%H:%M:%S") << "\n";
	for (int jx=0; jx<ntot; ++jx) {
	    const double x = jx/inv_a;
	    f << x << "\n";
	    for (int jy=0; jy<ntot; ++jy) {
		Result& res = result[jx*ntot + jy];
		f << std::format("{:9f} {:12e} {:24.16} {:24.16}",
				 res.y, res.d, res.worst_x, res.worst_y) << "\n";
	    }
	    f << "\n";
	}
	std::cout << "wrote " << fname << std::endl;
    }
    return 0;
}
