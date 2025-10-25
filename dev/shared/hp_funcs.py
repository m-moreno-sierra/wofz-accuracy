# File hp_funcs.py:
#   Python module providing high-precision functions related to w(z),
#   for use in code generation for libcerf.
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
#   Joachim Wuttke, Forschungszentrum Jülich, 2024
#
# Revision history:
#   oct25 Initial version published in libcerf/dev.

from mpmath import *

def cagree(u, v):
    """
    Returns true if complex numbers agree with each other within epsilon.
    """
    m = (u+v)/2
    d = u-v
    t = 1e-17
    return abs(d.real) < t * (abs(m.real) + t*abs(m)) and abs(d.imag) < t * (abs(m.imag) + t*abs(m))

def erfcx(x, doublecheck=False):
    result = exp(x**2)*erfc(x)
    if doublecheck:
        # Check mpmath-computed reference value against mpmath-based brute-force integration
        if x<6:
            r2 = exp(x**2) - 2/sqrt(pi) * quad(lambda t : exp(x**2-t**2), [0, x])
        else:
            r2 = 2/sqrt(pi) * quad(lambda t : exp(x**2-t**2), [x, mpf('inf')])
        if abs(result-r2)/result > 1e-17:
            raise Exception(f"mpmath inaccurate")
    return result

def imwx(x, doublecheck=False):
    fz = exp(-x**2)*erfc(mpc(0, -x))
    result = fz.imag
    if doublecheck:
        # Check mpmath-computed reference value against mpmath-based brute-force integration
        r2 = 2/sqrt(pi) * exp(-x**2) * quad(lambda t : exp(t**2), [0, x])
        if abs(result-r2)/result > 1e-17:
            raise Exception(f"mpmath inaccurate")
    return result

def wofz(z, doublecheck=False):
    j = mpc('0', '1')
    r1 = exp(-z**2)*erfc(-j*z)
    if (not doublecheck) or z.imag==0:
        return r1
    # Check mpmath-computed reference value against mpmath-based brute-force integration
    r2 = mpc(0,1)/pi*quad(lambda t: exp(-t**2)/(z-t), [-inf, +inf])
    if cagree(r1, r2):
        return r1
    x = z.real
    r3 = mpc(0,1)/pi*(quad(lambda t: exp(-t**2)/(z-t), [-inf, x])
                      + quad(lambda t: exp(-t**2)/(z-t), [x, +inf]))
    if cagree(r1, r3):
        return r1
#            w = z
#            for k in reversed(range(2000)):
#                w = z - k/2/w
#            r4 = mpc(0,1)/sqrt(pi)/w
    raise Exception(f"mpmath inaccurate for z=%8g+i%8g: r1=%8g+i%8g, d2=%8g+i%8g, d3=%8g+i%8g" %
                    (z.real, z.imag, r1.real, r1.imag,
                     (r2-r1).real, (r2-r1).imag, (r3-r1).real, (r3-r1).imag))

def wofz_taylor(z, N):
    """
    Taylor coefficients of w(z), forward computed.
    """
    W = []
    W.append(wofz(z, True))
    W.append(-2*z*W[0] + mpc(0,2)/sqrt(pi))
    for k in range(2,N):
        W.append(-2*(z*W[k-1]+(k-1)*W[k-2]))
    R = []
    fac = 1
    for k in range(N):
        R.append(W[k] * fac)
        fac /= (k+1)
    return R
