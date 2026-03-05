import subprocess
import flint

result=subprocess.check_output(["../build/run/run_wofz", "1.0", "0.5"]).decode().split()
w_c=flint.acb(result[0], result[1])

flint.ctx.prec = 96 
z = flint.acb(1.0, 0.5)
i = flint.acb(0, 1)
iz = i * z
w_ref = (-z**2).exp() * (-iz).erfc()

error_abs = abs(w_c - w_ref.mid())
error_rel = error_abs/abs(w_ref.mid())

print(f"Relative error: {error_rel.str(10)}")

