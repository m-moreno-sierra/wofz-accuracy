import matplotlib.pyplot as plt

re_z = [1.9915, 2.6911, 3.2353, 3.6973, 4.1061, 4.4768, 4.8185, 5.1371, 5.4367, 5.7204]
im_z = [-1.3548, -2.1770, -2.7844, -3.2874, -3.7259, -4.1196, -4.4798, -4.8138, -5.1265, -5.4216]

plt.figure(figsize=(8, 6))
plt.scatter(re_z, im_z, color='red', marker='x', s=100, label='Zeros $z_n$')

for i, txt in enumerate(range(1, 11)):
    plt.annotate(f"n={txt}", (re_z[i], im_z[i]), xytext=(5,5), textcoords='offset points')

plt.axhline(0, color='black', linewidth=1)
plt.axvline(0, color='black', linewidth=1)
plt.grid(True, linestyle='--', alpha=0.7)
plt.xlabel('$Re(z)$')
plt.ylabel('$Im(z)$')
plt.title('Distribution of $w(z)$ Zeros in the Lower Half-Plane')
plt.legend()

plt.savefig('zeros_map.png', dpi=300)
