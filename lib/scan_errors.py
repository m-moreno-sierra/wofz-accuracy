import subprocess
import flint

def main():
    flint.ctx.prec = 96
    i = flint.acb(0, 1)
    epsilon = 2.22e-16

    print("--- y > 0 ---")
    for x, y in [(1.23, 4.56), (5.67, 8.91), (9.01, 0.12)]:
        out = subprocess.check_output(["../build/run/run_wofz", str(x), str(y)]).decode().split()
        w_c = flint.acb(out[0], out[1])
        z = flint.acb(x, y)
        w_ref = (-z**2).exp() * (-(i * z)).erfc()
        ratio = (abs(complex(w_c) - complex(w_ref.mid())) / abs(complex(w_ref.mid()))) / epsilon
        print(f"x={x:<5} | y={y:<5} | Error/eps: {ratio:.2f}")

    print("\n--- y < 0 ---")
    zeros = [
        (1.99146684, -1.35481013), 
        (3.61361182, -1.50201383),
        (4.83515863, -1.60216447)
    ]
    
    for x, y in zeros:
        out = subprocess.check_output(["../build/run/run_wofz", str(x), str(y)]).decode().split()
        w_c = flint.acb(out[0], out[1])
        z = flint.acb(x, y)
        w_ref = (-z**2).exp() * (-(i * z)).erfc()
        error = abs(complex(w_c) - complex(w_ref.mid())) / abs(complex(w_ref.mid()))
        print(f"x={x:<10} | y={y:<10} | Rel. Error: {error:.2e}")

if __name__ == "__main__":
    main()