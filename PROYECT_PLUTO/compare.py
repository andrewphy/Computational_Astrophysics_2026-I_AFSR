

import matplotlib
matplotlib.use("Agg")    
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import sys, os

# directorio donde estan los datos de PLUTO y los scripts de Python
JET_DIR = "/home/pipe-desktop/pluto-4.4-patch3/PLUTO/Test_Problems/MHD/Jet"
# directorio de salida para todas las figuras
OUT_DIR = os.path.join(JET_DIR, "resultados")
os.makedirs(OUT_DIR, exist_ok=True)   # crea el directorio si no existe

# agrega JET_DIR al path para poder importar mhd_jet_muscl
sys.path.insert(0, JET_DIR)

# frames de PLUTO a comparar (con dbl=1.0 en pluto.ini, frame N = tiempo N)
FRAMES = [1, 4, 8]

# =============================================================================
# 1. SIMULACION PYTHON
# =============================================================================
print("=" * 60)
print("  Simulacion Python  MUSCL + RK2")
print("=" * 60)

import mhd_jet_muscl as m   # importa el solver Python con todos sus parametros

# ejecuta la simulacion y guarda snapshots en t=1,4,8
# run() retorna una lista de tuplas (tiempo, array_Q)
snaps = m.run(outdir=os.path.join(OUT_DIR, "py_snaps"), save_times=(1.0, 4.0, 8.0))

def py_snap(t_target):
    """Devuelve el snapshot (t, Q) cuyo tiempo es mas cercano a t_target."""
    return min(snaps, key=lambda s: abs(s[0] - t_target))

# =============================================================================
# 2. CARGA DE DATOS PLUTO
# =============================================================================
print("\nCargando datos PLUTO ...")
import pyPLUTO as pp   # libreria oficial de lectura de datos PLUTO

def load_pluto(frame):
    """Carga el frame N de PLUTO desde el directorio JET_DIR."""
    return pp.Load(nout=frame, path=JET_DIR, text=False)

def pluto_var(d, *names):
    """
    Busca una variable en el objeto de datos PLUTO por varios nombres posibles.
    PLUTO puede guardar Bphi como 'Bx3' o 'bx3' segun la version.
    """
    for n in names:
        if hasattr(d, "d_vars") and n in d.d_vars:
            return d.d_vars[n]   # diccionario de variables del nuevo pyPLUTO
        if hasattr(d, n):
            return getattr(d, n)  # atributo directo del objeto (version vieja)
    raise KeyError("Variable no encontrada: " + str(names))

# =============================================================================
# 3. FUNCIONES AUXILIARES PARA LOS PLOTS
# =============================================================================

def mirror(arr):
    """
    Refleja un array 2D (NR, NZ) en la direccion R para obtener (-R_max..R_max, NZ).
    El jet es axisimetrico, entonces el hemisferio R<0 es imagen del R>0.
    arr[::-1, :] invierte el orden en R (de R_max a 0), concatenado con arr (de 0 a R_max).
    """
    return np.concatenate([arr[::-1, :], arr], axis=0)

def mirror_R(r):
    """
    Construye el eje R completo (negativo y positivo) a partir del semi-eje r > 0.
    -r[::-1] son los valores negativos en orden decreciente de |r|.
    """
    return np.concatenate([-r[::-1], r])

def log_norm(arrays):
    """
    Calcula una normalizacion LogNorm comun a varias arrays.
    vmin = minimo positivo (escala log no acepta ceros)
    vmax = maximo global
    Usar la misma norma en PLUTO y Python permite comparar colores directamente.
    """
    vmin = min(a[a > 0].min() for a in arrays if (a > 0).any())
    vmax = max(a.max() for a in arrays)
    return mcolors.LogNorm(vmin=vmin, vmax=vmax)

def sym_norm(arrays):
    """
    Calcula una normalizacion lineal simetrica centrada en 0 para varias arrays.
    Usada para campos que cambian de signo (Bphi, velocidades).
    """
    am = max(np.abs(a).max() for a in arrays)   # amplitud maxima
    return mcolors.Normalize(vmin=-am, vmax=am)

# =============================================================================
# 4. PLOTS 2D: rho y Bphi para cada tiempo
# =============================================================================
print("\nGenerando plots 2D ...")

for frame in FRAMES:
    # intenta cargar el frame de PLUTO (puede fallar si el archivo no existe)
    try:
        d = load_pluto(frame)
    except Exception as e:
        print("  skip frame " + str(frame) + ": " + str(e)); continue

    # tiempo real del frame: d.timelist es la lista de tiempos de todos los frames
    # d.timelist[frame] da el tiempo exacto del frame N (no necesariamente = N exacto)
    t_pluto = d.timelist[frame]

    # snapshot de Python mas cercano al tiempo del frame
    t_py, Q = py_snap(float(frame))

    # extrae las variables de interes
    rho_p  = pluto_var(d, "rho")       # densidad de PLUTO, shape (NR_pluto, NZ_pluto)
    rho_y  = Q[m.iRHO]                 # densidad de Python, shape (NR_py, NZ_py)
    bphi_p = pluto_var(d, "Bx3", "bx3")  # campo toroidal de PLUTO
    bphi_y = Q[m.iBF]                  # campo toroidal de Python

    # ejes R completos (reflejados) para PLUTO y Python
    R_p = mirror_R(d.x1)   # d.x1 son los centros de celda en R de PLUTO
    R_y = mirror_R(m.R)    # m.R son los centros de celda en R de Python

    # normalizacion comun: mismo rango de colores en PLUTO y Python
    norm_rho  = log_norm([mirror(rho_p),  mirror(rho_y)])
    norm_bphi = sym_norm([mirror(bphi_p), mirror(bphi_y)])

    # figura con 4 paneles: 2 filas (rho, Bphi) x 2 columnas (PLUTO, Python)
    fig, axes = plt.subplots(2, 2, figsize=(13, 14))
    fig.suptitle("Jet MHD  t=" + str(frame), fontsize=15, fontweight="bold")

    # lista de paneles: (eje, datos, R, Z, norma, colormap, etiqueta, titulo)
    pairs = [
        (axes[0,0], mirror(rho_p),  R_p, d.x2, norm_rho,  "inferno", "log rho",
         "PLUTO HLLD   t=" + str(round(t_pluto, 2))),
        (axes[0,1], mirror(rho_y),  R_y, m.Z,  norm_rho,  "inferno", "log rho",
         "Python MUSCL+RK2   t=" + str(round(t_py, 2))),
        (axes[1,0], mirror(bphi_p), R_p, d.x2, norm_bphi, "RdBu_r",  "Bphi",
         "PLUTO HLLD   t=" + str(round(t_pluto, 2))),
        (axes[1,1], mirror(bphi_y), R_y, m.Z,  norm_bphi, "RdBu_r",  "Bphi",
         "Python MUSCL+RK2   t=" + str(round(t_py, 2))),
    ]

    for ax, data, Rv, Zv, norm, cmap, label, title in pairs:
        # pcolormesh dibuja una malla coloreada: Rv (NR_full,), Zv (NZ,), data (NR_full, NZ)
        # data.T porque pcolormesh espera (NZ, NR_full) con los ejes como se definen
        im = ax.pcolormesh(Rv, Zv, data.T,
                           cmap=plt.colormaps[cmap], norm=norm, shading="auto")
        plt.colorbar(im, ax=ax, label=label)   # barra de color con unidades
        ax.set_title(title, fontsize=11)
        ax.set_xlabel("R"); ax.set_ylabel("z")

    plt.tight_layout()
    # nombre del archivo: compare_t01.png, compare_t04.png, compare_t08.png
    fname = os.path.join(OUT_DIR, "compare_t" + str(frame).zfill(2) + ".png")
    plt.savefig(fname, dpi=150, bbox_inches="tight"); plt.close()
    print("  -> " + fname)

# =============================================================================
# 5. PERFILES AXIALES rho(z) en R=0
# =============================================================================
print("\nGenerando perfiles axiales ...")

cols = ["tab:blue", "tab:orange", "tab:green"]   # un color por tiempo
fig, axes = plt.subplots(1, 3, figsize=(15, 5))
fig.suptitle("Perfil axial  rho(z)  en R=0  PLUTO vs Python",
             fontsize=13, fontweight="bold")

l2_list = []   # guardara los errores L2 para imprimirlos al final

for ax, frame, col in zip(axes, FRAMES, cols):
    rho_p = None   # se usara para verificar si PLUTO cargo correctamente

    # perfil de PLUTO: celda R=0 es el indice 0 en la primera dimension
    try:
        d     = load_pluto(frame)
        t_p   = d.timelist[frame]          # tiempo real del frame
        rho_p = pluto_var(d, "rho")[0, :]  # rho en R=0 para todos los z
        # semilogy: escala logaritmica en Y para ver varios ordenes de magnitud
        ax.semilogy(d.x2, rho_p, color=col, lw=2.5,
                    label="PLUTO  t=" + str(round(t_p, 2)))
    except Exception as e:
        print("  skip frame " + str(frame) + ": " + str(e))

    # perfil de Python
    t_py, Q = py_snap(float(frame))
    rho_y   = Q[m.iRHO][0, :]   # rho en R=0 (primera celda radial) para todos z
    ax.semilogy(m.Z, rho_y, color=col, lw=2, ls="--",
                label="Python MUSCL  t=" + str(round(t_py, 2)))

    ax.set_xlabel("z"); ax.set_ylabel("rho")
    ax.legend(fontsize=9); ax.grid(True, alpha=0.3)

    # calcula el error L2 entre Python y PLUTO
    if rho_p is not None:
        # interpola rho de Python a la grilla de PLUTO (mas fina) para comparar
        # np.interp(x_nuevo, x_conocido, y_conocido)
        rho_yi = np.interp(d.x2, m.Z, rho_y)
        # error L2 = sqrt(media de los errores al cuadrado)
        l2 = np.sqrt(np.mean((rho_yi - rho_p)**2))
        l2_list.append((frame, l2))
        ax.set_title("t=" + str(frame) + "  L2=" + str(round(l2, 3)), fontsize=11)
    else:
        ax.set_title("t=" + str(frame), fontsize=11)

plt.tight_layout()
fname = os.path.join(OUT_DIR, "perfiles_rho.png")
plt.savefig(fname, dpi=150, bbox_inches="tight"); plt.close()
print("  -> " + fname)

# =============================================================================
# 6. RESUMEN DE ERRORES
# =============================================================================
print("\n" + "=" * 60)
print("  Errores L2  rho(z) en R=0")
for frame, l2 in l2_list:
    print("  t=" + str(frame) + "  L2=" + str(round(l2, 4)))
print("  Plots en: " + OUT_DIR)
print("=" * 60)
