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
    if (argc != 5) {
        fprintf(stderr, "Usage: %s d2 inverse_a ix iy\n", argv[0]);
        return 1;
    }
    char *endptr;
    int d2 = strtol(argv[1], &endptr, 10);
    assert(*endptr == '\0');
    double inv_a = strtod(argv[2], &endptr);
    assert(*endptr == '\0');
    int ix = strtol(argv[3], &endptr, 10);
    assert(*endptr == '\0');
    int iy = strtol(argv[4], &endptr, 10);
    assert(*endptr == '\0');

    const auto now = std::chrono::system_clock::now();
    printf("# Created by %s at %s\n", argv[0], std::format("{:%d-%m-%Y %H:%M:%OS}", now).c_str());
    printf("# d2 = %i\n", d2);
    printf("# 1/a = %g\n", inv_a);
    printf("# z = %g + i %g\n", ix*2*inv_a, iy*2*inv_a);

    const double wm = wmin(ix, iy, 2*inv_a, d2, 600);

    std::cout << std::format("wm={:f}", wm) << std::endl;

    return 0;
}
