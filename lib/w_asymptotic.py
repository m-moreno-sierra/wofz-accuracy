import flint
import sys
import math

def get_w_reference(z_flint):
    i = flint.acb(0, 1)
    return (-z_flint**2).exp() * (-(i * z_flint)).erfc()

def w_asymptotic_order(z_flint, N):
    i = flint.acb(0, 1)
    coeffs = [0.5, 0.75, 1.875]
    series = flint.acb(1)
    for m in range(1, N + 1):
        series += coeffs[m-1] / (z_flint**(2*m))
    
    w_alg = (i / (flint.arb.pi().sqrt() * z_flint)) * series
    
    if z_flint.imag < 0:
        reflection = 2 * (-z_flint**2).exp()
        return reflection + w_alg
    else:
        return w_alg

def main():
    flint.ctx.prec = 128
    
    zeros = [
        (1.991466842834, -1.354810128112),
        (2.691149024251, -2.177044906090),
        (3.235330868353, -2.784387613230),
        (3.697309702468, -3.287410789390),
        (4.106107284683, -3.725948719446),
        (4.476815692968, -4.119635227612),
        (4.818488291883, -4.479832797731),
        (5.137067271266, -4.813806682044),
        (5.436703910734, -5.126531545497),
        (5.720434851015, -5.421588576923)
    ]

    print(f"{'n':<2} & {'Re(zn)':<15} & {'Im(zn)':<15} & {'N=1':<6} & {'N=2':<6} & {'N=3':<6} \\\\")
    print("-" * 75)

    for n, (x, y) in enumerate(zeros, 1):
        z = flint.acb(x, y)
        results = []
        for N in [1, 2, 3]:
            w_asy = w_asymptotic_order(z, N)
            val = float(w_asy.abs_upper())
            digits = -math.log10(val) if val > 0 else 16.0
            results.append(f"{digits:>6.2f}")
        
        print(f"{n:<2} & {x:<15.12f} & {y:<15.12f} & {' & '.join(results)} \\\\")

if __name__ == "__main__":
    main()
