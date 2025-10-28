#!/bin/env python

"""
Generates a map relerr(x,y), where relerr is the worst out of n random draws in square at x,y.
"""

from mpmath import *
import os, sys
dir_path = os.path.dirname(os.path.realpath(__file__))
sys.path.insert(0, dir_path+'/../shared')
import hp_funcs as hp
import runtool as rt
import random

mp.dps = 48
mp.pretty = True

def run_cerf(z0):
    return rt.external_function2d(z0.real, z0.imag)

if __name__ == '__main__':
    if len(sys.argv)<3:
        print(f"Usage: {sys.argv[0]} <n> <inverse_a>")
        print(f"   where <n> is number of random draws per square")
        sys.exit(-1)
    NDraw = int(sys.argv[1])
    inverse_a = int(sys.argv[2])
    Rmap = 10

    rt.external_program = "run/run_wofz"
    rt.output_mode = 't'

    for jx in range(inverse_a*Rmap):
        xj = jx/inverse_a
        print(xj)
        for jy in range(inverse_a*Rmap):
            yj = jy/inverse_a

            worst = 0
            for i in range(NDraw):
                dx = random.uniform(0, 1/inverse_a)
                dy = random.uniform(0, 1/inverse_a)

                z = mpc(xj+dx, yj+dy)
                z, wc, a, n = run_cerf(z)
                wr = hp.wofz(z, True)
                re = abs((wr-wc)/wr)
                if re > worst:
                    worst = re

            print(yj, worst)
        print()
