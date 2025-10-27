#include <cassert>
#include <iostream>
#include <cstdio>
#include <cstdlib>
#include <chrono>
#include <format>
#include <future>
#include <limits>
#include "wn.h"

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
