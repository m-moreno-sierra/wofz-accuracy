#!/bin/env python

"""
Computes tau for grid points c_ij and given N, delta.
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

def wmin(ix, iy, inverse_b, d2):
    # print(f"WMIN ix={ix} iy={iy} d2={d2}")
    wmin = inf
    for dix in range(ix%2, int(sqrt(d2)), 2):
        diy = iy%2 + int((sqrt(d2-dix**2)-iy%2)/2)
        wmin = min(wmin, abs(hp.wofz(mpc((ix+dix)/inverse_b, (iy+diy)/inverse_b))))
        #print(f"... dix={dix} diy={diy} wmin={'%12g' % wmin}")
    #print(f"RETURN wmin={'%12g' % wmin}")
    return wmin

if __name__ == '__main__':
    if len(sys.argv)!=5:
        print(f"Usage: {sys.argv[0]} N_Taylor delta inverse_a n_offsets")
        sys.exit(1)
    N = int(sys.argv[1])
    delta = float(sys.argv[2])
    inverse_a = float(sys.argv[3])
    nr = int(sys.argv[4])

    fut.print_provenience()
    print(f"# N_Taylor = {N}")
    print(f"# delta = {delta}")
    print(f"# M_recenter = {nr}")
    print(f"# 1/a = {inverse_a}")

    d2max = 600
    S = [[sorted_diameters(d2max, sx, sy) for sy in [0, 1]] for sx in [0, 1]]

    inverse_b = 2*inverse_a

    for ix in range(int(7*inverse_b)):
        x0 = ix / inverse_b
        print(x0)

        for iy in range(int(7*inverse_b)):
            y0 = iy / inverse_b
            if x0**2+y0**2>=49:
                continue
            Sxy = S[ix%2][iy%2]

            nGenS = int(math.log(len(Sxy), 2)) # number of binary generations
            imin = 0
            x,y = x0,y0

            z = mpc(x, y)
            T = dw.forward(z, 60, False)
            assert(rt.rho_of_cNt(z, N, 0, T) <= delta)

            itau = 0
            minerr = inf
            for n in reversed(range(nGenS)):
                inew = itau + 2**n
                if inew<imin or inew>=len(Sxy):
                    continue
                d2 = Sxy[inew]
                tau = sqrt(d2) / inverse_b
                te = rt.truncation_error(z, N, tau, T)
                if isinf(te):
                    continue
                re = rt.rounding_error(z, N, tau, T)
                err = (te+re) / wmin(ix, iy, inverse_b, d2)
                if err <= delta:
                    itau = inew
                    minerr = err
            if itau >= len(Sxy)-1:
                raise Exception("increase d2max!")

            if nr > 0:
                imin = itau
                rx = range(-nr,nr+1)
                if ix==0:
                    rx = [0]
                ry = range(-nr,nr+1)
                if iy==0:
                    ry = [0]
                d2 = Sxy[itau]
                tau = sqrt(d2) / inverse_b
                wmi = wmin(ix, iy, inverse_b, d2)
                for jx in rx:
                    for jy in ry:
                        if jx==0 and jy==0:
                            continue
                        xj = math.nextafter(x0, math.copysign(math.inf, jx), steps=abs(jx))
                        yj = math.nextafter(y0, math.copysign(math.inf, jy), steps=abs(jy))
                        zj = mpc(xj, yj)
                        T = dw.forward(zj, 60, False)
                        err = rt.abserr(zj, N, tau, T) / wmi
                        if err < minerr:
                            minerr = err
                            z = zj

                T = dw.forward(z, 60, False)
                assert(rt.rho_of_cNt(z, N, 0, T) <= delta)
                itau = 0
                for n in reversed(range(nGenS)):
                    inew = itau + 2**n
                    if inew<imin:
                        itau = inew
                        continue
                    if inew>=len(Sxy):
                        continue
                    d2 = Sxy[inew]
                    tau = sqrt(d2) / inverse_b
                    te = rt.truncation_error(z, N, tau, T)
                    if isinf(te):
                        continue
                    re = rt.rounding_error(z, N, tau, T)
                    err = (te+re) / wmin(ix, iy, inverse_b, d2)
                    if err <= delta:
                        itau = inew
                if itau >= len(Sxy)-1:
                    raise Exception("increase d2max!")

            if isfinite(minerr):
                print("%10.5g %3i %s %s" %
                      (y0, Sxy[itau], fut.double2hexstring(z.real), fut.double2hexstring(z.imag)))
            else:
                print("%10.5g   0 0x0 0x0" % (y0))
        print()
