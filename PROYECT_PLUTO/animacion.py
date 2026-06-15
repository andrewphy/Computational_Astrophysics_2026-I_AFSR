import matplotlib
matplotlib.use("Agg")
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import matplotlib.animation as animation
import sys, os

JET_DIR = "/home/pipe-desktop/pluto-4.4-patch3/PLUTO/Test_Problems/MHD/Jet"
sys.path.insert(0, JET_DIR)
import pyPLUTO as pp

# carga todos los frames disponibles (0..8)
frames = []
for n in range(9):
    try:
        d = pp.Load(nout=n, path=JET_DIR, text=False)
        rho = d.d_vars["rho"] if "rho" in d.d_vars else d.rho
        R_full = np.concatenate([-d.x1[::-1], d.x1])
        rho_full = np.concatenate([rho[::-1, :], rho], axis=0)
        frames.append((d.timelist[n], R_full, d.x2, rho_full))
    except Exception as e:
        print(f"frame {n}: {e}")

# rango de color global
all_vals = np.concatenate([f[3].ravel() for f in frames])
vmin = all_vals[all_vals > 0].min()
vmax = all_vals.max()
norm = mcolors.LogNorm(vmin=vmin, vmax=vmax)

# figura
fig, ax = plt.subplots(figsize=(5, 9))
t0, R0, Z0, rho0 = frames[0]
im = ax.pcolormesh(R0, Z0, rho0.T, cmap="inferno", norm=norm, shading="auto")
plt.colorbar(im, ax=ax, label="log ρ")
title = ax.set_title(f"PLUTO HLLD   t = {t0:.2f}", fontsize=12)
ax.set_xlabel("R"); ax.set_ylabel("z")
plt.tight_layout()

def update(i):
    t, R, Z, rho = frames[i]
    im.set_array(rho.T.ravel())
    title.set_text(f"PLUTO HLLD   t = {t:.2f}")
    return im, title

ani = animation.FuncAnimation(fig, update, frames=len(frames), interval=500, blit=True)
fname = os.path.join(JET_DIR, "resultados", "pluto_rho.gif")
os.makedirs(os.path.dirname(fname), exist_ok=True)
ani.save(fname, writer="pillow", fps=2, dpi=130)
plt.close()
print(f"-> {fname}")
