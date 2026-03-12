import sys
import subprocess
import flint

def main():
    if len(sys.argv) < 3:
        print("Usage: python3 wofz_error.py <x> <y>")
        return

    x_val = sys.argv[1]
    y_val = sys.argv[2]

    result = subprocess.check_output(["../build/run/run_wofz", x_val, y_val]).decode().split()
    w_c = flint.acb(result[0], result[1])

    flint.ctx.prec = 96 #aprox 24 digits
    z = flint.acb(x_val, y_val)
    i = flint.acb(0, 1)
    iz = i * z
    w_ref = (-z**2).exp() * (-iz).erfc()

    error = abs(w_c - w_ref.mid()) / abs(w_ref.mid())

    print(f"Relative error: {error.str(10)}")

if __name__ == "__main__":
    main()