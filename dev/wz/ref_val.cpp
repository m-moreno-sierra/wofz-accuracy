//  ************************************************************************************************
//
//  Project:   Code generating code for Taylor expansion on square tiles
//
//  Reference: Joachim Wuttke,
//             Code generation for computing an analytical function with near machine precision
//             on square tiles, with application to the Faddeeva function
//             (in preparation, available upon request)
//
//  File:      ref_val.cpp
//
//  Purpose:   Program ref_val prints high-precision reference value for given z=x+iy.
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
#include "wn.h"

//! Program ref_val prints high-precision reference function value f(z)
//! for z=x+iy given by command-line arguments x, y.

int main(int argc, char *argv[])
{
    if (argc != 3) {
        fprintf(stderr, "Print a reference value f(z) for given z=x+iy.\n", argv[0]);
        fprintf(stderr, "Usage: x y\n", argv[0]);
        return 1;
    }
    char *endptr;
    double x = strtod(argv[1], &endptr);
    assert(*endptr == '\0');
    double y = strtod(argv[2], &endptr);
    assert(*endptr == '\0');

    const std::complex<double> fz = wofz(x, y);

    std::cout << std::format("f(z) = {:21.16e} + i {:21.16e}", fz.real(), fz.imag()) << std::endl;

    return 0;
}
