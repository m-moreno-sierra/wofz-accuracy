//  ************************************************************************************************
//
//  Project:   Code generating code for Taylor expansion on square tiles
//
//  Reference: Joachim Wuttke,
//             Code generation for computing an analytical function with near machine precision
//             on square tiles, with application to the Faddeeva function
//             (in preparation, available upon request)
//
//  File:      statistics.cpp
//
//  Purpose:   Provides program statistics.
//
//  License:   MIT License (file LICENSE in toplevel directory)
//  Copyright: Forschungszentrum Jülich GmbH 2025
//  Author:    Joachim Wuttke <j.wuttke@fz-juelich.de>
//
//  ************************************************************************************************

#include "relerr.h"
#include <cassert>
#include <chrono>
#include <cstdio>
#include <cstdlib>
#include <format>
#include <fstream>
#include <iostream>
#include <omp.h>
#include <random>

int main(int argc, char* argv[])
{
    if (argc != 5) {
        fprintf(stderr, "Usage: %s #region range #channels #draws\n", argv[0]);
        return 1;
    }
    char* endptr;
    int region = strtol(argv[1], &endptr, 10);
    assert(*endptr == '\0');
    double range = strtod(argv[2], &endptr);
    assert(*endptr == '\0');
    int nch = strtol(argv[3], &endptr, 10);
    assert(*endptr == '\0');
    double ndraws = strtod(argv[4], &endptr);
    assert(*endptr == '\0');

    const double R = 7;
    const double dx = range * pow(2, -53) / nch;

    std::vector<int> histo(nch, 0);


    const int nthreads = 40;

#pragma omp parallel for
    for (int j = 0; j < nthreads; ++j) {

        std::mt19937 reng(j);
        std::uniform_real_distribution<double> rdis(0, 1);

        for (int idraws = 0; idraws < ndraws / nthreads;) {
            double x = rdis(reng);
            double y = rdis(reng);
            if (region == 1) {
                x *= .23;
                y *= .23;
                const double z2 = x * x + y * y;
                if (z2 >= 0.053)
                    continue;
            } else if (region == 2 || region == 4) {
                x *= 7;
                y *= 7;
                const double z2 = x * x + y * y;
                if (z2 >= R * R || z2 <= 0.053)
                    continue;
            } else if (region == 3 || region == 5) {
                x *= 14;
                y *= 14;
                const double z2 = x * x + y * y;
                if (z2 <= R * R)
                    continue;
            } else
                assert(0);

            const double relerr = relerr_wofz(x, y);

            const int k = int(relerr / dx);
            ++idraws;
            ++histo[k >= nch ? nch - 1 : k];
        }
    }

    auto now = std::time(nullptr);
    std::cout << "# Created by " << argv[0] << "region=" << region << " range=" << range
              << " #ch=" << nch << "#draws=" << ndraws << " on "
              << std::put_time(std::localtime(&now), "%H:%M:%S") << "\n";
    for (int k = 0; k < nch; ++k)
        std::cout << std::format("{:12f} {:12d}", (k + 0.5) * range / nch, histo[k]) << "\n";
    if (histo[nch - 1] != 0)
        std::cout << "# Warning: nonzero overflow channel\n";
    return 0;
}
