#!/bin/env python

"""
Print relative-error estimate for the asymptotic expansion of w(z).
"""

from mpmath import *
import os, sys
dir_path = os.path.dirname(os.path.realpath(__file__))
sys.path.insert(0, dir_path+'/../shared')
import runtool as rt
import hp_funcs as hp

mp.dps = 64
mp.pretty = True

eps = 2**(-53)
I = mpc(0,1)
la = sqrt(5)

def wmin(x):
    z = mpc(x, 0)
    w = hp.wofz(z)
    return abs(w)

def err_asymex(z, N):
    NN = 30
    assert(N<NN)
    An = [gamma(mpf(n)+1/2)/pi for n in range(NN)]
    Vn = [abs(a - float(a))/abs(a)/eps for a in An]
    Bn = [Vn[n] + 6 + n*(3+la) for n in range(NN)]
    Bn[0] -= 1
    Bn[N-1] -= 1

    re = 0
    for n in range(N):
        re += Bn[n] * An[n] / z**(2*n+1)
    wmi = wmin(z)
    ae = An[N]  / z**(2*N+1) / eps

    print("z    = %9g" % (z))
    print("wmin = %9g" % (wmi))
    print("re   = %9g" % (re))
    print("ae   = %9g" % (ae))
    print("tre  = %9g" % ((re+ae) / wmi))

if __name__ == '__main__':
    if len(sys.argv)!=3:
        print(f"Usage: {sys.argv[0]} |z| N")
        sys.exit(1)
    z = float(sys.argv[1])
    N = int(sys.argv[2])

    err_asymex(z, N)
