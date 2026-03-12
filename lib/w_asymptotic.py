import flint
import sys

def w_asymptotic(z_val):
   
    flint.ctx.prec = 96
    z = flint.acb(z_val.real, z_val.imag)
    i = flint.acb(0, 1)
    
    term0 = flint.acb(1)
    term1 = 0.5 / (z**2)
    term2 = 0.75 / (z**4)
    term3 = 1.875 / (z**6)
    
    series = term0 + term1 + term2 + term3
    
    w_alg = (i / (flint.arb.pi().sqrt() * z)) * series
    
    if z.imag < 0:
        reflection = 2 * (-z**2).exp()
        return reflection + w_alg
    else:
        return w_alg

def main():
    if len(sys.argv) < 3:
        print("Usage: python3 w_asymptotic.py <x> <y>")
        return

    x, y = float(sys.argv[1]), float(sys.argv[2])
    z = complex(x, y)
    
    result = w_asymptotic(z)
    print(f"Asymptotic w({x} + {y}i) = {result.str(24)}")

if __name__ == "__main__":
    main()