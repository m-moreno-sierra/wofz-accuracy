import subprocess
import flint

flint.ctx.prec = 96
points = [("1.0", "0.5"), ("5.0", "0.001"), ("0.1", "10.0")]

print(f"{'Point z':<20} | {'Relative Error (TRE)':<20}")
print("-" * 45)

for x, y in points:
    try:
        raw_output = subprocess.check_output(["../build/run/run_wofz", x, y]).decode().split()
        w_c = flint.acb(raw_output[0], raw_output[1])

        z = flint.acb(float(x), float(y))
        i = flint.acb(0, 1)
        w_ref = (-z**2).exp() * (-(i * z)).erfc()

        diff = w_c.mid() - w_ref.mid()
        tre = abs(diff) / abs(w_ref.mid())

        print(f"{x} + {y}j {' ':<8} | {float(tre):.4e}") 
    except Exception as e:
        print(f"Error processing z = {x} + {y}j: {e}")
