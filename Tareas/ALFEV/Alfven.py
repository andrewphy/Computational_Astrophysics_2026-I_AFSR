"""
Animacion - Onda de Alfven Circularmente Polarizada 1D
"""
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation

# ─── Parametros ────────────────────────────────────────────────────────────
N      = 128
L      = 1.0
dx     = L / N
x      = np.linspace(dx/2, L - dx/2, N)
rho0   = 1.0
P0     = 0.1
B0     = 1.0
eps    = 0.1
mu0    = 1.0
CFL    = 0.4
t_stop = 5.0

# ─── Condicion inicial ─────────────────────────────────────────────────────
phi0 = 2.0 * np.pi * x
rho = np.ones(N) * rho0
vx  = np.zeros(N)
vy  = eps * np.sin(phi0)
vz  = eps * np.cos(phi0)
Bx  = np.ones(N) * B0
By  = eps * np.sin(phi0)
Bz  = eps * np.cos(phi0)
P   = np.ones(N) * P0

def ddx(f):
    return (np.roll(f, -1) - np.roll(f, 1)) / (2.0 * dx)

def rhs(vy, vz, By, Bz):
    dvy_dt = (B0 / (mu0 * rho0)) * ddx(By)
    dvz_dt = (B0 / (mu0 * rho0)) * ddx(Bz)
    dBy_dt = B0 * ddx(vy)
    dBz_dt = B0 * ddx(vz)
    return dvy_dt, dvz_dt, dBy_dt, dBz_dt

def paso_dt():
    return CFL * dx / 1.0  # vA = 1

# ─── Simulacion guardando snapshots ────────────────────────────────────────
snapshots = []
snap_times = []
t = 0.0
step = 0
save_every = 8  # guardar cada N pasos

print("Corriendo simulacion...")
while t < t_stop:
    dt = min(paso_dt(), t_stop - t)

    d1 = rhs(vy, vz, By, Bz)
    vy_s = vy + dt * d1[0]
    vz_s = vz + dt * d1[1]
    By_s = By + dt * d1[2]
    Bz_s = Bz + dt * d1[3]

    d2 = rhs(vy_s, vz_s, By_s, Bz_s)
    vy = vy + 0.5 * dt * (d1[0] + d2[0])
    vz = vz + 0.5 * dt * (d1[1] + d2[1])
    By = By + 0.5 * dt * (d1[2] + d2[2])
    Bz = Bz + 0.5 * dt * (d1[3] + d2[3])

    t += dt
    step += 1

    if step % save_every == 0:
        snapshots.append((By.copy(), Bz.copy(), vy.copy()))
        snap_times.append(t)

print(f"Simulacion terminada. {len(snapshots)} frames guardados.")

# ─── Animacion ────────(pura IA pero se ve lindo)───────────────────────────────────────────────────
fig, axes = plt.subplots(2, 1, figsize=(10, 6))
fig.suptitle('Onda de Alfvén CP - 1D MHD Ideal (RK2)', fontsize=13)

ax1, ax2 = axes
ax1.set_xlim(0, L)
ax1.set_ylim(-eps*1.5, eps*1.5)
ax1.set_ylabel('By', fontsize=11)
ax1.grid(True, alpha=0.3)
ax1.plot(x, eps*np.sin(phi0), 'k--', alpha=0.3, label='t=0')

ax2.set_xlim(0, L)
ax2.set_ylim(-eps*1.5, eps*1.5)
ax2.set_ylabel('Bz', fontsize=11)
ax2.set_xlabel('x', fontsize=11)
ax2.grid(True, alpha=0.3)
ax2.plot(x, eps*np.cos(phi0), 'k--', alpha=0.3, label='t=0')

line_By, = ax1.plot([], [], 'b-', lw=2, label='By(t)')
line_Bz, = ax2.plot([], [], 'r-', lw=2, label='Bz(t)')
time_text = ax1.text(0.02, 0.88, '', transform=ax1.transAxes, fontsize=10,
                     bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))
ax1.legend(loc='upper right', fontsize=9)
ax2.legend(loc='upper right', fontsize=9)

def init():
    line_By.set_data([], [])
    line_Bz.set_data([], [])
    time_text.set_text('')
    return line_By, line_Bz, time_text

def animate(i):
    By_i, Bz_i, _ = snapshots[i]
    line_By.set_data(x, By_i)
    line_Bz.set_data(x, Bz_i)
    time_text.set_text(f't = {snap_times[i]:.2f}')
    return line_By, line_Bz, time_text

ani = animation.FuncAnimation(fig, animate, init_func=init,
                               frames=len(snapshots), interval=30, blit=True)

plt.tight_layout()
ani.save('alfven_animation.gif', writer='pillow', fps=30, dpi=100)
print("Guardado: alfven_animation.gif")