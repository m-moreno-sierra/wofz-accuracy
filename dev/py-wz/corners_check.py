#!/bin/env python

"""
Check w(z) around corners of square tiling.
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

def check_at(x,y):
    z = mpc(float(x), float(y))
    z, wc, g, t = run_cerf(z)
    wr = hp.wofz(z, False)

    dr = 0
    if wr.real != 0:
        dr = abs((wc-wr).real/wr.real)
    di = 0
    if wr.imag != 0:
        di = abs((wc-wr).imag/wr.imag)
    da = abs((wc-wr)/wr)

    return x, y, dr, di, da, g, t

if __name__ == '__main__':
    Dr = 0
    Di = 0
    Da = 0
    zr = None
    zi = None
    za = None
    h = 1+1e-14
    n = 0
    inverse_a = 7
    R = 7
    for ix in range(inverse_a*R):
        xc= ix/8
        for iy in range(inverse_a*R):
            yc= iy/8
            if (xc==0 and yc==0) or hypot(xc,yc)>=R:
                continue
            for x,y in ((xc/h,yc/h), (xc*h,yc/h), (xc/h,yc*h), (xc*h,yc*h)):
                if hypot(x,y)<.23 or hypot(x,y)>=R:
                    continue
                _x, _y, dr, di, da, g, t = check_at(x, y)
                n += 1
                assert(_x == x)
                assert(_y == y)
                if g != 900:
                    continue
                if dr > Dr:
                    Dr = dr
                    zr = mpc(x, y)
                    print("worst rel re  of %8.4e at %21.16f %21.16f subdomain %2i" % (Dr, zr.real, zr.imag, t))
                if di > Di:
                    Di = di
                    zi = mpc(x, y)
                    print("worst rel im  of %8.4e at %21.16f %21.16f subdomain %2i" % (Di, zi.real, zi.imag, t))
                if da > Da:
                    Da = da
                    za = mpc(x, y)
                    print("worst rel abs of %8.4e at %21.16f %21.16f subdomain %2i" % (Da, za.real, za.imag, t))
    print(f"Did {n} tests.")
