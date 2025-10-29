//  ************************************************************************************************
//
//  Project:   Code generating code for Taylor expansion on square tiles
//
//  Reference: Joachim Wuttke,
//             Code generation for computing an analytical function with near machine precision
//             on square tiles, with application to the Faddeeva function
//             (in preparation, available upon request)
//
//  File:      runat.cpp
//
//  Purpose:   Provides program runat.
//
//  License:   MIT License (file LICENSE in toplevel directory)
//  Copyright: Forschungszentrum Jülich GmbH 2025
//  Author:    Joachim Wuttke <j.wuttke@fz-juelich.de>
//
//  ************************************************************************************************

#include <cassert>
#include <iostream>
#include <cstdio>
#include <cstdlib>
#include <chrono>
#include <format>
#include <future>
#include <limits>
#include "terms.h"

int main(int argc, char *argv[])
{
    if (argc != 6) {
        fprintf(stderr, "Usage: %s N_Taylor d2 inverse_a x y\n", argv[0]);
        return 1;
    }
    char *endptr;
    int N = strtol(argv[1], &endptr, 10);
    assert(*endptr == '\0');
    int d2 = strtol(argv[2], &endptr, 10);
    assert(*endptr == '\0');
    double inv_a = strtod(argv[3], &endptr);
    assert(*endptr == '\0');
    double x = strtod(argv[4], &endptr);
    assert(*endptr == '\0');
    double y = strtod(argv[5], &endptr);
    assert(*endptr == '\0');

    const auto now = std::chrono::system_clock::now();
    printf("# Created by %s at %s\n", argv[0], std::format("{:%d-%m-%Y %H:%M:%OS}", now).c_str());
    printf("# N_Taylor = %i\n", N);
    printf("# d2 = %i\n", d2);
    printf("# 1/a = %g\n", inv_a);
    printf("# z = %g + i %g\n", x, y);

    const double R = 7;
    const double inv_b = 2 * inv_a;
    const double tau = sqrt(d2)/(2*inv_a);
    const std::complex<double> z{x, y};
    const std::vector<Coeff> WN = w_n_vector(x, y);

    const double te = Terms::truncation_error(z, N, tau, WN);
    const double re = Terms::rounding_error(z, N, tau, WN);
    std::cout << std::format("te={:f} re={:f}", te, re) << std::endl;

    return 0;
}
