#!/bin/env python

"""
For given N, compute smallest z for which the remainder R_N
of the asymptotic expansion of w(z) falls below zeta*eps*|w(z)|.
Protokollbuch sc1, 158.
"""

from mpmath import *
import sys
sys.path.insert(0, '../shared')
import hp_funcs as hp

mp.dps = 64
mp.pretty = True

eps = 2**(-53)
zeta = 1./25
I = mpc(0,1)

def wmin(x):
    z = mpc(x, 0)
    w = hp.wofz(z)
    return abs(w)

def run(N):
    x = ( gamma(N+1./2)/sqrt(pi)/zeta/eps )**(1./(2*N)) # first order of asymptotic expansion
    # iterate using exact wmin(x) instead of the above estimate
    for iter in range(100):
        x = (gamma(N+1./2)/pi/zeta/eps/wmin(x))**(1./(2*N+1))
    print(N, x)

if __name__ == '__main__':
    for N in range(1,50):
        run(N)
