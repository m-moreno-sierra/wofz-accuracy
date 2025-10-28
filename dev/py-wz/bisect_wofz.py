#!/bin/env python

"""
Computes accuracy of external implementation of function w_of_z
by comparison with an internal high-precision function.

The comparison is done for z=r2z(r,s) with r,s from some grids,
and for additional points obtained by bisection around points r
where the employed algorithm or number of terms has changed.

The grids are hard-coded here, and are changed ad hoc during development.

Output is either a table with one line for each r,
or a report on the worst case (the r with largest inaccuracy of f(x)).
"""

from mpmath import *
import os, sys
dir_path = os.path.dirname(os.path.realpath(__file__))
sys.path.insert(0, dir_path+'/../shared')
import hp_funcs as hp
import runtool as rt

mp.dps = 48
mp.pretty = True

rt.print_line = lambda locus, a, n, rr, f, f2: print('%2i %3i %3i  %22.16e  %22.16e %22.16e %22.16e  %8e %8e %8e' %
          (locus, a, n, rr, re(f), im(f), abs(f), re(f-f2)/abs(f2), im(f-f2)/abs(f2), abs((f-f2)/f2)))


if __name__ == '__main__':
    if len(sys.argv)<2:
        print(f"Usage: {sys.argv[0]} <scan_mode>")
        print(f"   where <scan_mode>   is any of b (bisect) s (scan)")
        sys.exit(-1)
    if sys.argv[1] == 'b':
        do_scan = rt.scan_and_bisect
    elif sys.argv[1] == 's':
        do_scan = rt.scan_wo_bisect
    else:
        raise Exception("Invalid arg1, expected b|s")

    rt.external_program = "run/run_wofz"
    rt.output_mode = 't'

#    S = rt.loggrid(111, .2, 8)
#    R = rt.loggrid(111, .2, 8)
    R = rt.lingrid(173, 0, 10.25)
    S = rt.lingrid(173, 0, 10.25)

    for s in S:
        print(s)
        def ext(r):
            a = rt.external_function2d(float(s), float(r))
            return (a[0].imag, a[1], a[2], a[3])
        rt.f_ext = ext
        rt.f_hp = lambda r : hp.wofz(mpc(float(s), float(r)))
        do_scan(R)
        print()
