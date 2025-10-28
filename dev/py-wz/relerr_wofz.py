#!/bin/env python

"""
Computes relative error of cerf{|w(z)|} on a grid, and prints worst case.
"""

from mpmath import *
import sys
sys.path.insert(0, '../shared')
import hp_funcs as hp
import runtool as rt

mp.dps = 48
mp.pretty = True

if __name__ == '__main__':
    rt.external_program = "run/run_wofz"
    rt.range_mode = 'p'
    rt.output_mode = 't'

#    S = rt.loggrid(111, .2, 8)
#    R = rt.loggrid(111, .2, 8)
    R = rt.lingrid(111, 0, 7.75)
    S = rt.lingrid(111, 0, 7.75)

    A = [7, 7.5, 8, 9, 10, 12, 15, 20, 23, 23.5, 30, 50, 100, 150, 300, 2000, 50000, 6e7, 7e7, 9e7]

    for ia in range(1, len(A)):
        wc_re = 0
        wc_z = None
        for az in rt.lingrid(1371, A[ia-1], A[ia]):
            for ph in rt.lingrid(1373, 0, pi/2):
                x = az*cos(ph)
                y = az*sin(ph)
                z = mpc(x, y)
                zback, v, d1, d2 = rt.external_function2d(x, y)
                assert(abs((zback-z)/z) < 2.3e-16)
                vref = hp.wofz(z)
                relerr = abs((v-vref)/vref)
                if relerr > wc_re:
                    wc_re = relerr
                    wc_z = z
                    print("%2i %6g %6g %12g %12g %12g" %
                          (ia, A[ia-1], A[ia], wc_z.real, wc_z.imag, wc_re/2**-53))
