import matplotlib.pyplot as plt
import numpy as np

n = np.arange(1, 11)
err_libcerf = [8.33e-17, 1.89e-16, 3.75e-16, 1.19e-16, 2.73e-16, 
               3.93e-16, 5.72e-16, 6.25e-16, 2.81e-16, 1.12e-16]
err_asymp = [1.25e-03, 5.10e-05, 7.82e-06, 2.08e-06, 7.47e-07, 
             3.24e-07, 1.60e-07, 8.72e-08, 5.10e-08, 3.16e-08]

plt.figure(figsize=(10, 6))

plt.semilogy(n, err_libcerf, 'o-', label='libcerf Error (Absolute)', color='#005493', linewidth=2) # Azul Jülich
plt.semilogy(n, err_asymp, 's--', label='Asymptotic Expansion (A&S 7.1.23)', color='#f07d00', linewidth=2) # Naranja

plt.grid(True, which="both", ls="-", alpha=0.5)
plt.xlabel('Zero Index (n)', fontsize=12)
plt.ylabel('Absolute Error $|w(z_n)|$', fontsize=12)
plt.title('Numerical Audit of the Faddeeva Function Zeros', fontsize=14, fontweight='bold')
plt.xticks(n)
plt.legend(fontsize=11)

plt.axhline(y=2.22e-16, color='gray', linestyle=':', alpha=0.7)
plt.text(1, 2.5e-16, 'Machine Epsilon ($\epsilon$)', color='gray', fontsize=10)


plt.tight_layout()

plt.savefig('audit_comparison.png', dpi=300)
