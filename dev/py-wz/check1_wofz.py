#!/bin/env python

"""
Check w(z) for one z.
"""

from mpmath import *
import os, sys
dir_path = os.path.dirname(os.path.realpath(__file__))
sys.path.insert(0, dir_path+'/../shared')
import hp_funcs as hp
import runtool as rt

mp.dps = 48
mp.pretty = True
rt.external_program = "run/run_wofz"

def run_cerf(z0):
    return rt.external_function2d(z0.real, z0.imag)

def print_deviation(txt, wc, wr):
    print("%s: im %8g re %8g abs %8g" %
          (txt, abs((wc-wr).real/wr.real), abs((wc-wr).imag/wr.imag), abs((wc-wr)/wr)))

if __name__ == '__main__':
    if len(sys.argv)<3:
        print(f"Usage: {sys.argv[0]} <x> <y>")
    z = mpc(sys.argv[1], sys.argv[2])
    z, wc, a, n = run_cerf(z)
    wr = hp.wofz(z, True)

    dr = 0
    if wr.real != 0:
        dr = abs((wc-wr).real/wr.real)
    di = 0
    if wr.imag != 0:
        di = abs((wc-wr).imag/wr.imag)

    print("z:      %22.15e %22.15e" % (z.real, z.imag))
    print("mpmath: %22.15e %22.15e" % (wr.real, wr.imag))
    print("cerf:   %22.15e %22.15e algo %i nterms %i rdiff %7.2e %7.2e %7.2e" %
          (wc.real, wc.imag, a, n, dr, di, abs((wc-wr)/wr)))
