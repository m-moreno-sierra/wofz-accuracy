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
    if len(sys.argv)!=5:
        print(f"Usage: {sys.argv[0]} N_Taylor d2 x y")
        sys.exit(1)
    N = int(sys.argv[1])
    d2 = int(sys.argv[2])
    x = float(sys.argv[3])
    y = float(sys.argv[4])
    inv_a = 10

    d2max = 600
    assert(d2<=d2max)
    S = [[sorted_diameters(d2max, sx, sy) for sy in [0, 1]] for sx in [0, 1]]

    z = mpc(x, y)
    T = dw.forward(z, 60, False)
    tau = sqrt(d2) / (2* inv_a)
    te = rt.truncation_error(z, N, tau, T)
    re = rt.rounding_error(z, N, tau, T)
    wm = wmin(int(x*2*inv_a), int(y*2*inv_a), int(2*inv_a), d2)
    print(f"x = {x}")
    print(f"y = {y}")
    print(f"d2 = {d2}")
    print(f"tau = {tau}")
    print(f"te = {te}")
    print(f"re = {re}")
    print(f"wm = {wm}")
    print(f"rh = {(te+re)/wm}")
