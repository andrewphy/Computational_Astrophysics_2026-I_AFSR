import sys
sys.path.insert(0, '/home/pipe-desktop/pluto-4.4-patch3/PLUTO/Tools/pyPLUTO')
import pyPLUTO as pp
import numpy as np
import matplotlib.pyplot as plt

# Cargar datos de PLUTO
w_dir = '/home/pipe-desktop/pluto-4.4-patch3/PLUTO/Test_Problems/MHD/CP_Alfven/'
D = pp.pload.pload(0, w_dir=w_dir, datatype='vtk')

# Solucion analitica - onda viaja a vA=1, en t=1 recorre L=1
t = 1.0
eps = 0.1
phi = 2.0 * np.pi * (D.x1 - t)  # onda viajando hacia la derecha
By_analitica = eps * np.sin(phi)

# Error L1
error = np.abs(D.Bx2 - By_analitica)
error_L1 = np.mean(error)

# Grafica
fig, axes = plt.subplots(2, 1, figsize=(9, 7))

axes[0].plot(D.x1, D.Bx2, label='PLUTO By')
axes[0].plot(D.x1, By_analitica, '--', label='ANAL By')
axes[0].set_ylabel('By')
axes[0].legend()
axes[0].grid(True, alpha=0.3)

axes[1].plot(D.x1, error)
axes[1].set_ylabel('|error|')
axes[1].set_xlabel('x')
axes[1].set_title(f'Error L1 = {error_L1:.2e}')
axes[1].grid(True, alpha=0.3)

plt.tight_layout()
plt.savefig('alfven_error.png', dpi=150)
print(f"Error L1 = {error_L1:.2e}")
print("Guardado: alfven_error.png")