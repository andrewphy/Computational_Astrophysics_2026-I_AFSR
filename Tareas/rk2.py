import numpy as np
from matplotlib import pyplot as plt


# Función diferencial
def f(t, I):
    S = 0.3
    return -I + S

# Condiciones iniciales
t0 = 0
I0 = 1
h = 0.01

# Tiempo final
tf = 5

# Listas para guardar resultados
t_vals = [t0]
I_vals = [I0]

t = t0
I = I0
        
def euler(f, t0, y0, h, tf):
    t, y = t0, y0
    t_vals, y_vals = [t], [y]

    while t < tf:
        y = y + h*f(t, y)
        t += h
        t_vals.append(t)
        y_vals.append(y)

    return t_vals, y_vals

def rk2_midpoint(f, t0, y0, h, tf):
    t = t0
    y = y0

    t_vals = [t]
    y_vals = [y]

    while t < tf:
        k1 = f(t, y)
        k2 = f(t + h/2, y + (h/2)*k1)

        y = y + h * k2
        t += h

        t_vals.append(t)
        y_vals.append(y)

    return t_vals, y_vals

def rk4_runge(f, t0, y0, h, tf):
    t = t0
    y = y0

    t_vals = [t]
    y_vals = [y]

    while t < tf:
        k1 = f(t, y)
        k2 = f(t + h/2, y + (h/2)*k1)
        k3 = f(t + h/2, y + (h/2)*k2)
        k4 = f(t + h, y + h*k3)

        y = y + (h/6)*(k1 + 2*k2 + 2*k3 + k4)
        t += h

        t_vals.append(t)
        y_vals.append(y)

    return t_vals, y_vals

def rk4_kutta(f, t0, y0, h, tf):
    t = t0
    y = y0

    t_vals = [t]
    y_vals = [y]

    while t < tf:
        k1 = f(t, y)
        k2 = f(t + h/3, y + (h/3)*k1)
        k3 = f(t + 2*h/3, y - (h/3)*k1 + h*k2)
        k4 = f(t + h, y + h*k1 - h*k2 + h*k3)

        y = y + (h/8)*(k1 + 3*k2 + 3*k3 + k4)
        t += h

        t_vals.append(t)
        y_vals.append(y)

    return t_vals, y_vals



trkm_vals, Irkm_vals = rk2_midpoint(f, t0, I0, h, tf)
teuler_vals, Ieuler_vals = euler(f, t0, I0, h, tf)
trk4r_vals, Irk4r_vals = rk4_runge(f, t0, I0, h, tf)
trk4k_vals, Irk4k_vals = rk4_kutta(f, t0, I0, h, tf)

plt.plot(trkm_vals, Irkm_vals,linewidth=1, label='RK2 Midpoint')
plt.plot(teuler_vals, Ieuler_vals, linewidth=1, label='Euler')
plt.plot(trk4r_vals, Irk4r_vals,  linewidth=1, label='RK4 Runge-Kutta')
plt.plot(trk4k_vals, Irk4k_vals,  linewidth=1, label='RK4 Kutta')
plt.axhline(0.5, color='red', linestyle='--', label='Concentración 0.5')
plt.xlabel('Tiempo (s)')
plt.ylabel('I')
plt.title(rf'Solve $\frac{{dI}}{{dt}}=-I+S$')
plt.legend()
plt.grid()
plt.show()

# --- Evaluar I(tau) resolviendo la ODE con RK4 ---
def I_of_tau(tau_target, h=0.01):
    S = 0.3

    def f_local(t, I):
        return -I + S

    t = 0.0
    I = 1.0

    while t < tau_target:
        h_step = min(h, tau_target - t)

        k1 = f_local(t, I)
        k2 = f_local(t + h_step/2, I + (h_step/2)*k1)
        k3 = f_local(t + h_step/2, I + (h_step/2)*k2)
        k4 = f_local(t + h_step, I + h_step*k3)

        I = I + (h_step/6)*(k1 + 2*k2 + 2*k3 + k4)
        t += h_step

    return I


# --- Bisección pura ---
def bisection(tol=1e-6, max_iter=100):
    a, b = 0.0, 5.0

    def g(tau):
        return I_of_tau(tau) - 0.5

    fa = g(a)
    fb = g(b)

    if fa * fb > 0:
        raise ValueError("No hay cambio de signo")

    for _ in range(max_iter):
       
        c = (a + b)/2
        fc = g(c)

        if abs(fc) < tol:
            return c

        if fa * fc < 0:
            b = c
            fb = fc
        else:
            a = c
            fa = fc

    return (a + b)/2

tau_half = bisection()
print(f"Tau_m (bisección pura) = {tau_half:.6f}")