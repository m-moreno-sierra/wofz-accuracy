import ctypes
import flint
import math
import os

PATH_TO_SO = "/home/miriam/libcerf/build/libc/libcerf.so.3"

class ComplexDouble(ctypes.Structure):
    _fields_ = [("real", ctypes.c_double), ("imag", ctypes.c_double)]

try:
    lib = ctypes.CDLL(PATH_TO_SO)
    w_func = lib.w_of_z
    w_func.restype = ComplexDouble
    w_func.argtypes = [ctypes.c_double, ctypes.c_double]
except Exception as e:
    print(f"Error cargando libcerf: {e}")
    lib = None

flint.ctx.prec = 96


def w_asymptotic_T23(z):
    i_unit = flint.acb(0, 1)
    series = 1 + 0.5/(z**2) + 0.75/(z**4) + 1.875/(z**6)
    w_theo = (i_unit / (flint.acb.pi().sqrt() * z)) * series
    return 2 * (-z**2).exp() + w_theo

def get_exact_zero(n):
    pi = flint.arb.pi()
    lamb = ((n - 0.125) * pi).sqrt()
    mu = (2 * (2 * pi).sqrt() * lamb).log()
    xn, yn = -lamb + 0.25*mu/lamb, lamb + 0.25*mu/lamb
    
    z = flint.acb(xn, yn)
    deriv_const = -2/flint.arb.pi().sqrt()
    for _ in range(15):
        f = z.erfc()
        f_prime = deriv_const * (-z**2).exp()
        z = z - f/f_prime
        
    re_part = z.real
    im_part = z.imag
    z_conj = flint.acb(re_part, -im_part)
    
    return flint.acb(0, 1) * z_conj


def main():
    if not lib: return

    print(f"{'n':<2} | {'Re(z)':<10} | {'Im(z)':<10} | {'Err libcerf':<12} | {'Err Asymptotic':<12}")
    print("-" * 65)
    
    for n in range(1, 11):
        z_exact = get_exact_zero(n)
        re_val = float(z_exact.real)
        im_val = float(z_exact.imag)
        
        res = w_func(ctypes.c_double(re_val), ctypes.c_double(im_val))
        err_lib = math.sqrt(res.real**2 + res.imag**2)
        
        w_theo = w_asymptotic_T23(z_exact)
        err_theo = float(w_theo.abs_upper())
        
        print(f"{n:<2} | {re_val:<10.4f} | {im_val:<10.4f} | {err_lib:.2e}   | {err_theo:.2e}")

if __name__ == "__main__":
    main()