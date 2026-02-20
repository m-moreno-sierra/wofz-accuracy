#!/bin/env python

"""
Print absolute error and wmin, to test monotonicity.
"""

from mpmath import *
import os, sys
dir_path = os.path.dirname(os.path.realpath(__file__))
sys.path.insert(0, dir_path+'/../shared')
import hp_funcs as hp
import derive_w as dw
import functool as fut
import reltot as rt
import bisect, math

def sorted_diameters(N, sx, sy):
    D2 = set()
    for j in range(1+sy, N, 2):
        for i in range(1+sx, N, 2):
            d2 = j**2 + i**2
            if d2 > N:
                continue
            D2.add(d2)
    return sorted(D2)

def wmin(ix, iy, Nb, d2):
    # print(f"WMIN ix={ix} iy={iy} d2={d2}")
    wmin = inf
    for dix in range(ix%2, int(sqrt(d2)), 2):
        diy = iy%2 + int((sqrt(d2-dix**2)-iy%2)/2)
        wmin = min(wmin, abs(hp.wofz(mpc((ix+dix)/Nb, (iy+diy)/Nb))))
        #print(f"... dix={dix} diy={diy} wmin={'%12g' % wmin}")
    #print(f"RETURN wmin={'%12g' % wmin}")
    return wmin

if __name__ == '__main__':
    if len(sys.argv)!=4:
        print(f"Usage: {sys.argv[0]} N_Taylor delta ix")
        sys.exit(1)
    N = int(sys.argv[1])
    delta = float(sys.argv[2])
    ix = int(sys.argv[3])

    d2max = 600
    S = [[sorted_diameters(d2max, sx, sy) for sy in [0, 1]] for sx in [0, 1]]

    Nb = 16

    nT = int((1.2*Nb)**2) + 1

    x0 = ix / Nb
    for iy in range(113):
        y0 = iy / Nb
        print(y0)
        if x0**2+y0**2>49:
            break
        Sxy = S[ix%2][iy%2]

        nGenS = int(math.log(len(Sxy), 2)) # number of binary generations
        imin = 0
        x,y = x0,y0

        z = mpc(x, y)
        T = dw.forward(z, 60, False)
        assert(rt.rho_of_cNt(z, N, 0, T) <= delta)

        for d2 in Sxy:
            tau = sqrt(d2) / Nb
            te = rt.truncation_error(z, N, tau, T)
            re = rt.rounding_error(z, N, tau, T)
            wm = wmin(ix, iy, Nb, d2)
            print("%4i %24.16e %24.16e %24.16e %24.16e" % (d2, te, re, wm, (te+re)/wm))

        print()
