#!/bin/env python

"""
Contour of |w(z)|=const.
"""

from mpmath import *
import sys
sys.path.insert(0, '../shared')
import hp_funcs as hp
import rdp

mp.dps = 48
mp.pretty = True

def per_v(v):
    ret = []

    f = lambda z: abs(hp.wofz(z))
    ymin = -5
    x = findroot(lambda x: f(mpc(x,ymin))-v, [-20, -.01], solver='bisect')
    z = mpc(x, ymin)

    ret.append((z.real, z.imag))

    z0, phi0 = z, None
    step = .02
    while z.imag>=ymin:
        a = lambda phi: z0 + step*exp(mpc(0,phi))
        if phi0:
            phi = findroot(lambda phi: f(a(phi)) - v, phi0, solver='newton')
        else:
            phi = findroot(lambda phi: f(a(phi)) - v, [0, pi], solver='bisect')
        z = a(phi)
        if z.imag<ymin:
            x = findroot(lambda x: f(mpc(x,ymin)) - v, [0, 2*z.real], solver='bisect')
            z = mpc(x,ymin)
            ret.append((z.real, z.imag))
            break
        ret.append((z.real, z.imag))
        z0, phi0 = z, phi

    return ret

if __name__ == '__main__':
    for v in [.14,.165,.2,.25,.32,.4,.55,.75,1.,1.4,2.,3.,6.,15.,40.,100,250]:
        raw = per_v(v)
        out = rdp.rdp(raw, .01)

        print(v)
        for o in out:
            print("%10.5f %10.5f" % (o[0], o[1]))
        print()
