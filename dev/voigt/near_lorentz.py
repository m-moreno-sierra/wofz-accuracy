#!/usr/bin/env python3

# File near_gauss.py
#   Compute Taylor coefficients of hwhm_voigt for sigma/gamma<<1.
#
# Reference:
#   J Wuttke, Power series for the half width of the Voigt function
#   (in preparation)
#
# Copyright:
#   (C) 2025 Forschungszentrum Jülich GmbH
#
# Licence:
#   Permission is hereby granted, free of charge, to any person obtaining
#   a copy of this software and associated documentation files (the
#   "Software"), to deal in the Software without restriction, including
#   without limitation the rights to use, copy, modify, merge, publish,
#   distribute, sublicense, and/or sell copies of the Software, and to
#   permit persons to whom the Software is furnished to do so, subject to
#   the following conditions:
#
#   The above copyright notice and this permission notice shall be
#   included in all copies or substantial portions of the Software.
#
#   THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
#   EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
#   MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
#   NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE
#   LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION
#   OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION
#   WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
#
# Author:
#   Joachim Wuttke, Forschungszentrum Jülich, 2025
#
# Revision history:
#   oct25 Initial version published in libcerf/dev/voigt.

import sys
sys.path.insert(0, '../shared')
from fractions import Fraction
import functool as ft
import math
from decimal import Decimal, getcontext

getcontext().prec = 100  # adjust as needed

M=100

# --------------------
# Auxiliary functions
# --------------------

psi = [2, -2, 1, 0, -2, 2, -1, 0]

def b(n, k):
    m = 2*n + k
    return Fraction(math.comb(m, k) * psi[m%8], 4**(m//4))

# --------------------
# Recurrence for g_n
# --------------------
def compute_g():
    a = [1]
    for n in range(1, M+1):
        a.append(a[-1] * (2*n-1))

    G = [[], [None, Fraction(3,2)]]

    for n in range(2,M+1):
        G.append([None for n in range(n+1)])
        for k in range(2,n+1):
            sum = 0
            for j in range(1, n-k+2):
                sum += G[j][1] * G[n-j][k-1]
            G[n][k] = sum
        Sum = 0
        for m in range(n+1):
            kstart = 1
            if m==0:
                kstart += 1
            sum = 0
            for k in range(kstart, n-m+1):
                sum += b(m,k) * G[n-m][k]
            Sum += a[m] * sum
        G[n][1] = (a[n]*(b(n,0)-2*(-1)**n) + Sum) / 2
    return [1] + [GG[1] for GG in G[1:]]

# --------------------
# Compute (1+g)^2 coefficients
# --------------------
def compute_1plusg_sq(g):
    N=len(g)-1
    res=[Fraction(0)]*(N+1)
    res[0]=Fraction(1)
    for n in range(1,N+1):
        s=2*g[n]
        for k in range(1,n):
            s += g[k]*g[n-k]
        res[n]=s/2**n
    return res

# --------------------
# Output helpers
# --------------------
# Check whether i can be written as 2^j, and return j
def lb(i):
    """Return j such that i = 2**j, or None if i is not a power of 2."""
    if i <= 0:
        raise ValueError(f"{i} is not a positive integer, cannot be a power of 2.")
    if (i & (i - 1)) != 0:
        raise ValueError(f"{i} is not a power of 2.")
    # i is power of 2, compute j
    j = i.bit_length() - 1
    return j

def frac2float(f):
    return float(Decimal(f.numerator) / Decimal(f.denominator))

# --------------------
# Main driver
# --------------------
if __name__=="__main__":
    g=compute_g()
    h=compute_1plusg_sq(g)

    print("Coefficients [(s/2)^n](1+g)^2, for cross-check against WangZZ22")
    for n in range(0,M+1):
        print(f"{n:2d}  {h[n].numerator}/2^{lb(h[n].denominator)} = {frac2float(h[n]): .10f}")

    print("Coefficients g_n=[s^n]g")
    for n in range(0,M+1):
        print(f"g_<{n}> &= {g[n].numerator :32d}/2^{lb(g[n].denominator) :2d} &={frac2float(g[n]) :40.20f}\\\\")
    print("Ditto as floats")
    for n in range(0,M+1):
        print(f"{frac2float(g[n]): 40.24g} {n:2d}")
    print("Ditto as hexstring")
    for n in range(0,M+1):
        print(f"{ft.double2hexstring(frac2float(g[n])):24s} {n:2d}")

    print("Ratio |g_n/g_{n-1}|")
    for n in range(1,M+1):
        print(f"{n:2d}  {g[n]/g[n-1]: 16.10g}")

    s = 1./256
    print(f"Check expansion for s={s}")
    sum = 1
    sfac = s
    oldterm = 1
    for n in range(1,M+1):
        term = frac2float(g[n])*sfac
        if abs(term) > abs(oldterm):
            break
        sum += term
        oldterm = term
        sfac *= s
        print(f"{n:2d}  {term: .24e}  {sum: .24f}")
