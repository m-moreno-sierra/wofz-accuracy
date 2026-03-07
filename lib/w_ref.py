import sys
import flint

def main():
    if len(sys.argv) < 3:
        print("Uso: python3 w_ref.py <x> <y>")
        return

    x_val = sys.argv[1]
    y_val = sys.argv[2]

    flint.ctx.prec = 96 

    z = flint.acb(x_val, y_val)

    i = flint.acb(0, 1)
    iz = i * z
    wz = (-z**2).exp() * (-iz).erfc()

    print(wz)

if __name__ == "__main__":
    main()
