#!/bin/env python

# File near_gauss.py
#   Compute Taylor coefficients of hwhm_voigt for gamma/sigma<<1.
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

from mpmath import *
import sys
sys.path.insert(0, '../shared')
import derive_w as dw
import hp_funcs as hp
import functool as ft
import runtool as rt
from functools import lru_cache

mp.pretty = True

def voigt(q_by_sigma,gamma_by_sigma):
    '''Voigt function, for use in tests'''
    wq = hp.wofz(mpc(q_by_sigma/sqrt(2), gamma_by_sigma/sqrt(2))).real
    w0 = hp.wofz(mpc(0, gamma_by_sigma/sqrt(2))).real
    return wq / w0

def Hn(dps):
    mp.dps = dps
    M = 32
    h0 = sqrt(log(2))
    W = dw.forward(h0, M, True, 1e-41)
    A = [(-1)**n / (2 * gamma(n / 2 + 1)) for n in range(M)]
    I = mpc(0, 1)
    h1 = (W[1].imag + A[1]) / W[1].real + I

    H = [h0, h1]

    @lru_cache(maxsize=None)
    def F(n, N):
        if n == 2:
            return sum(H[m] * H[N - m] for m in range(1, N))
        return sum(H[m] * F(n - 1, N - m) for m in range(1, N - n + 2))

    for N in range(2, M):
        hN = A[N]
        for n in range(2, N + 1):
            hN -= (W[n] * F(n, N)).real
        hN /= W[1].real
        H.append(hN)

    return [h.real for h in H]

if __name__ == '__main__':
    # print("sqrt(2*ln2) = " + ft.double2hexstring(sqrt(2*log(2))))
    H  = Hn(41)
    H2 = Hn(44)

    print(f"Coefficients for use in libcerf:")
    for N in range(len(H)):
        relerr = (H[N]-H2[N])/H[N]
        g = H[N] / sqrt(2)**(N-1)
        print("%2i %s %s %13.2e" % (N-1, nstr(g,32,min_fixed = +inf, max_fixed = -inf),
                                    ft.double2hexstring(g), relerr))
        if relerr > 1e-16:
            raise ValueError("Insufficient precision")

    print(f"Coefficients for comparison with WaZZ22:")
    for N in range(len(H)):
        relerr = (H[N]-H2[N])/H[N]
        print("%2i  %s  %13.2e" % (N, nstr(H[N],32,min_fixed = +inf, max_fixed = -inf), relerr))
        if relerr > 1e-32:
            raise ValueError("Insufficient precision")

    print("Tests:")
    Nj = 30
#    for j in range(Nj):
#        t = (j+1.)/Nj*.8
    for t in [0.125, 0.25, .5]:
        h = H[0]
        for n in range(1,len(H)):
            h += H[n]*(t/sqrt(2))**n
            if (abs(voigt(sqrt(2)*h,t)-0.5) < 1e-17):
                break
        print("%9.5g %2i %21.16f %9g" % (t, n, voigt(sqrt(2)*h,t), voigt(sqrt(2)*h,t)-0.5))
