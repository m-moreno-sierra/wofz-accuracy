#!/bin/env python

"""
Find approximative solution of the set cover problem.
Reads file that provides d2(x,y).
Writes two C source files that are included by w_of_z.c.
"""

import os, sys
dir_path = os.path.dirname(os.path.realpath(__file__))
sys.path.insert(0, dir_path+'/../shared')
from mpmath import *
from math import sqrt
import hp_funcs as hp
import functool as fut
import enumerate_polyominoes as ep
import datetime, re

mp.dps = 48
mp.pretty = True

def sorted_diameters(N, sx, sy):
    """
    Returns sorted list of possible values of squared circumcircle diameter
    in units of a, for given parities of center.
    """
    D2 = set()
    for j in range(1+sy, N, 2):
        for i in range(1+sx, N, 2):
            d2 = j**2 + i**2
            if d2 > N:
                continue
            D2.add(d2)
    return sorted(D2)

def polyomino_pattern(d2, sx, sy):
    rows = []
    nj2 = int(sqrt(d2)) + 1
    for k in range(2-sx, nj2, 2):
        dx = sqrt(d2 - k**2)/2 - (sy)/2
        row = 2 * int(dx) + (sy)
        if row > 0:
            if k > sx:
                rows.insert(0, row)
            rows.append(row)
    return rows

def read_d2_file(fname):
    """
    Reads file that provides d2(x,y).
    Returns 1/a, N, D2 (squared diameters for b-lattice points), and Z (expansion center cx, cy).
    """
    with open(fname, 'r') as f:
        t = ''
        f.readline() # header
        assert(m := re.match(r'# N_Taylor = (\d+)', f.readline()))
        N = int(m.group(1))
        assert(m := re.match(r'# delta = (\d+)', f.readline()))
        assert(m := re.match(r'# M_recenter = (\d+)', f.readline()))
        assert(m := re.match(r'# 1/a = (\d+)', f.readline()))
        inverse_a = float(m.group(1))
        for line in f:
            assert(not line[0] =='#')
            t += line
    D2 = []
    Z  = []
    for block in t.split('\n\n'):
        if block=='':
            break
        a = block.split('\n')
        x = float(a[0])
        if x == 0:
            ix = 0
        else:
            ix += 1
            if round(2*inverse_a*x) != ix:
                raise Exception(f'Unexpected x entry')
        B = []
        zz = []
        for iy in range(len(a)-1):
            l = a[1+iy]
            ww = l.split()
            if len(ww) != 4:
                raise Exception(f'Unexpected data line')
            y = float(ww[0])
            if ix==0 and iy==0:
                pass
            else:
                if round(2*inverse_a*y) != iy:
                    raise Exception(f'Unexpected y entry')
            B.append(int(ww[1]))
            cx = float.fromhex(ww[2])
            cy = float.fromhex(ww[3])
            zz.append((cx, cy))
        if B:
            D2.append(B)
            Z.append(zz)
    return inverse_a, N, D2, Z

def covered_squares(D2, P, ix, iy, inverse_a):
    """
    Returns list of coordinates of square tiles that may be covered by Taylor expansion around b-lattice point ix,iy.
    """
    d2 = D2[ix][iy]
    if d2==0:
        raise Exception('Cannot add expansion as d2=0')
    pat = P[ix%2][iy%2][d2]

    ret = []
    mx = ix//2
    my = iy//2
    lx = len(pat)
    for nx in range(lx):
        jx = mx + nx - lx//2
        ly = pat[nx]
        for ny in range(ly):
            jy = ny + my - ly//2
            if jx>=0 and jy>=0 and (jx/inverse_a)**2+(jy/inverse_a)**2<49:
                ret.append((jx,jy))
    return ret

def add_expansion(F, C, Q, ix, iy):
    n = len(C)
    C.append( (ix, iy) )

    qs = Q[ix][iy]

    count = 0
    for jx, jy in qs:
        nold = F[jx][jy]
        if nold == -1:
            F[jx][jy] = n
            count += 1
        elif nold >= 0:
            ixold, iyold = C[nold]
            if ixold!=0 and iyold!=0:
                ijx, ijy = 2*jx+1, 2*jy+1
                if (ijx-ix)**2 + (ijy-iy)**2 < (ijx-ixold)**2 + (ijy-iyold)**2:
                    F[jx][jy] = n

    return count

def n_naked(qs, F):
    count = 0
    for jx, jy in qs:
        if F[jx][jy] == -1:
            count += 1
    return count

if __name__ == '__main__':
    if len(sys.argv) != 2:
        raise Exception(f'Usage: {sys.argv[0]} <file with x blocks with y tau d2 lines>')

    # Load d2(x,y) from file.
    fname = sys.argv[1]
    inverse_a, Ntay, D2, Z = read_d2_file(fname)
    d2max = max([max(line) for line in D2])
    tmax = int(sqrt(d2max)) + 1

    Nax  = int(inverse_a*7+1)
    C = [] # expansion centers

    S = [[sorted_diameters(d2max, sx, sy) for sy in [0, 1]] for sx in [0, 1]]
    P = [[{d2:polyomino_pattern(d2, sx, sy) for d2 in S[sx][sy]} for sy in [0, 1]] for sx in [0, 1]]
    Q = [[covered_squares(D2, P, ix, iy, inverse_a) for iy in range(len(D2[ix]))] for ix in range(len(D2))]

    # Set up field of base squares.
    # Value -1: not yet covered by polyomino.
    # Value -2: inner reserved area
    # Value -3: outer reserved area
    F = [[-1 for jy in range(Nax)] for jx in range(Nax)]
    nF = Nax**2
    for jx in range(Nax):
        for jy in range(Nax):
            if (jx/inverse_a)**2+(jy/inverse_a)**2 >= 49:
                F[jx][jy] = -3
                nF -= 1
            elif ((jx+1)/inverse_a)**2+((jy+1)/inverse_a)**2 < 0.053:
                F[jx][jy] = -2
                nF -= 1

    # Expand around points on y axis.
    while True:
        # Search first square not yet covered
        for jy0 in range(1, Nax):
            if F[0][jy0] == -1:
                break
        if F[0][jy0] == -3: # not in domain
            break
        best_n = 0
        iy = None
        # print("DEBUG", 2*jy0+1, 2*Nax, 2*jy0+2*tmax+1)
        for iy1 in range(2*jy0+1, min(2*Nax, 2*jy0+2*tmax+1)):
            if F[0][iy1//2] < -1 or iy1 >= len(D2[0]):
                break
            qs = Q[0][iy1]
            if (0, jy0) in qs and len(qs) > best_n:
                best_n = n_naked(qs, F)
                iy = iy1
        if iy is None:
            break
        nF -= add_expansion(F, C, Q, 0, iy)
        # print(f"iy={iy} -> nF={nF}")

    # Expand around points on x axis.
    while True:
        # Search first square not yet covered
        for jx0 in range(1, Nax):
            if F[jx0][0] == -1:
                break
        if F[jx0][0] == -3: # not in domain
            break
        best_n = 0
        ix = None
        # print("DEBUG", 2*jx0+1, 2*Nax, 2*jx0+2*tmax+1)
        for ix1 in range(2*jx0+1, min(2*Nax, 2*jx0+2*tmax+1)):
            if F[ix1//2][0] < -1 or ix1 >= len(D2):
                break
            qs = Q[ix1][0]
            if (jx0, 0) in qs and len(qs) > best_n:
                best_n = n_naked(qs, F)
                ix = ix1
        if ix is None:
            break
        nF -= add_expansion(F, C, Q, ix, 0)
        # print(f"ix={ix} -> nF={nF}")

    # Sort remaining squares by distance from (-1,0).
    JJ = []
    for jx in range(1, Nax):
        for jy in range(1, Nax):
            if F[jx][jy] == -1:
                JJ.append((jx, jy))
    JJ = sorted(JJ, key=lambda jj: hypot(jj[0]+20*inverse_a, jj[1]))

    # Create expansion points inside the quadrant.
    while nF>0:
        # Search first square not yet covered
        for jx0, jy0 in JJ:
            if F[jx0][jy0] == -1:
                break
        best_n = 0
        ix = None
        for ix1 in range(len(D2)):
            for iy1 in range(len(D2[ix1])):
                if F[ix1//2][iy1//2] < -1 or D2[ix1][iy1] == 0:
                    continue
                qs = Q[ix1][iy1]
                if (jx0, jy0) in qs and len(qs) > best_n:
                    nn = n_naked(qs, F)
                    if nn > best_n:
                        ix, iy = ix1, iy1
                        best_n = nn
        if ix is None:
            raise Exception("Could not cover all tiles")
        nF -= add_expansion(F, C, Q, ix, iy)
        # print(f"i={ix},{iy} -> nF={nF}")

    assert(nF==0)
    print(f"covered domain with {len(C)} subdomains")

    fname1 = "/tmp/w_taylor_centers.tab"
    with open(fname1, "w") as f:
        print("# Created by %s on %s" % (" ".join(sys.argv), datetime.datetime.now().time()),
              file=f)
        print("# %i expansion centers approximately at b-lattice points:" % len(C), file=f)
        for n in range(len(C)):
            ix, iy = C[n]
            print("%3i %3i" % (ix, iy), file=f)
    print(f"wrote {fname1}")

    fname2 = "/tmp/w_taylor_cover.c"
    with open(fname2, "w") as f:
        print("// Created by %s on %s" % (" ".join(sys.argv), datetime.datetime.now().time()),
              file=f)
        typ = "short"
        if len(C)<128:
            typ = "signed char"
        print("static const double inverseA = %f;" % inverse_a, file=f)
        print("static const int nXcover = %i;" % Nax, file=f)
        print("alignas(64) static const %s Cover[%i] = {" % (typ, Nax**2), file=f)
        for jx in range(Nax):
            for jy in range(Nax):
                print("%2i," % F[jx][jy], end="", file=f)
            print("", file=f)
        print("};", file=f)
    print(f"wrote {fname2}")

    fname3 = "/tmp/w_taylor_coeffs.c"
    with open(fname3, "w") as f:
        print("// Created by %s on %s" % (" ".join(sys.argv), datetime.datetime.now().time()),
              file=f)
        print("static const int NTay = %i;" % Ntay, file=f)
        print("alignas(64) static const double TaylorCoeffs[2*%i*%i] = {" % (Ntay+1,len(C)), file=f)
        for c in C:
            ix, iy = c
            x, y = Z[ix][iy]
            z = mpc(x, y)
            W = hp.wofz_taylor(z, Ntay)
            print(" %s, %s," % (fut.double2hexstring(x), fut.double2hexstring(y)), file=f);
            for w in W:
                print(" %s, %s," % (fut.double2hexstring(w.real), fut.double2hexstring(w.imag)),
                      end="", file=f)
            print(" // x=%8g y=%8g" % (x, y), file=f)
        print("};", file=f)
    print(f"wrote {fname3}")
