#!/bin/env python

"""
Print reference value of w(z) for one z.
"""

from mpmath import *
import sys
sys.path.insert(0, '../shared')
import hp_funcs as hp

mp.dps = 48
mp.pretty = True

if __name__ == '__main__':
    if len(sys.argv)<3:
        print(f"Usage: {sys.argv[0]} <x> <y>")
    z = mpc(sys.argv[1], sys.argv[2])
    w = hp.wofz(z)
    print("%22.15e %22.15e" % (w.real, w.imag))
