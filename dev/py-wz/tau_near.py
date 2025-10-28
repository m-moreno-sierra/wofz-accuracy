#!/bin/env python

"""
Computes tau for given c, N, delta, and also for points near c.
"""

from mpmath import *
import sys
import reltot as rt
import derive_w as dw
import math

if __name__ == '__main__':
    if len(sys.argv)!=3:
        raise Exception(f'Usage: {sys.argv[0]} zre zim')

    xc = float(sys.argv[1])
    yc = float(sys.argv[2])

    N = 20
    delta = 3

    for jx in range(-3,4):
        for jy in range(-3,4):
            x = math.nextafter(xc, math.copysign(math.inf, jx), steps=abs(jx))
            y = math.nextafter(yc, math.copysign(math.inf, jy), steps=abs(jy))
            z = mpc(x, y)
            T = dw.forward(z, 60, False)
            tau = findroot(lambda tau: rt.rho_of_cNt(z, N, max(0,tau), T) - delta, 0.5, solver='newton')
            print("%i %i %12.5g" % (jx, jy, tau))
