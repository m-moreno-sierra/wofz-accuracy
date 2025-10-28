#!/bin/env python

"""
Print expansion coefficients of w(z) for the asymptotic expansion and for the Maclaurin series.
Used to generate code for direct inclusion in w_of_z.c.
"""

from mpmath import *

mp.dps = 64
mp.pretty = True

eps = 2**(-53)
I = mpc(0,1)

if __name__ == '__main__':
    print("# Asymptotic expansion")
    for n in reversed(range(23)):
        c = gamma(mpf(n)+1/2)/pi
        print(" + %22.16e) / (z*z) // n=%i" % (c, n))

    print("# Maclaurin series")
    for n in reversed(range(1,27)):
        c = 1 / gamma(mpf(n)/2+1)
        if n%4==2 or n%4==3:
            c = -c
        if n%2==0: # real
            print(" + C(%+22.16e, 0) ) * z // z^%i" % (c, n))
        else:
            print(" + C(0, %+22.16e) ) * z // z^%i" % (c, n))
    print(" + 1;")
