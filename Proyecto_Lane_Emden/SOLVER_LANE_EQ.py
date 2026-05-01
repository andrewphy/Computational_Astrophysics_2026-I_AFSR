import numpy as np
from scipy.interpolate import CubicSpline
from scipy.interpolate import CubicSpline
import matplotlib.pyplot as plt

class LaneEmdenSolver:
    """
    Solver para la ecuación de Lane-Emden usando RK4
    """

    def __init__(self, n, h=0.001, xi_max=12):
        self.n = n  #INDICE POLITROPICO
        self.h = h  #PASO
        self.xi_max = xi_max #VALOR MAX RADIO (\chi)

        # Resultados
        self.xi = None 
        self.theta = None
        self.dtheta = None
        self.xi1 = None #theta(xi1)=0

    # ============================================================
    # ECUACIONES DEL SISTEMA
    # ============================================================

    def rhs(self, xi, y):
        theta, dtheta = y

        if xi == 0:
            return np.array([dtheta, 0])

        # evitar problemas numéricos , y theta negativa (no tiene sentido físico)
        if theta < 0:
            theta = 0

        d2theta = -theta**self.n - (2/xi)*dtheta

        return np.array([dtheta, d2theta])

    # ============================================================
    # PASO RK4
    # ============================================================

    def rk4_step(self, xi, y): #rk4 
        h = self.h

        k1 = self.rhs(xi, y)
        k2 = self.rhs(xi + h/2, y + h*k1/2)
        k3 = self.rhs(xi + h/2, y + h*k2/2)
        k4 = self.rhs(xi + h, y + h*k3)

        return y + (h/6)*(k1 + 2*k2 + 2*k3 + k4)

    # ============================================================
    # CONDICIONES INICIALES (TAYLOR)
    # ============================================================

    def initial_conditions(self):
        xi = 1e-3 #\chi inicial (cercano de cero, se usa taylor)

        theta = 1 - xi**2/6 + self.n*xi**4/120
        dtheta = -xi/3 + self.n*xi**3/30

        return xi, np.array([theta, dtheta])

    # ============================================================
    # SOLVER PRINCIPAL
    # ============================================================

    def solve(self):

        xi, y = self.initial_conditions()

        xi_vals = [0, xi]  #inicializamos con el punto inicial (cercano a cero), el primer punto es exactamente xi=0, theta=1, dtheta=0 por las condiciones de iniciales
        theta_vals = [1, y[0]]
        dtheta_vals = [0, y[1]]

        while xi < self.xi_max:
            y_new = self.rk4_step(xi, y)
            xi += self.h

            xi_vals.append(xi)
            theta_vals.append(y_new[0])
            dtheta_vals.append(y_new[1])

        # ======================================================
        # DETECCIÓN DE LA RAÍZ + CUBICSPLINE
        # ======================================================
            if y_new[0] <= 0:
                xi_sub = np.array(xi_vals[-5:]) # tomar últimos puntos (mejor usar varios)
                theta_sub = np.array(theta_vals[-5:])

                # spline cúbico
            
                cs = CubicSpline(xi_sub, theta_sub) #es una función que interpola los puntos xi_sub, theta_sub con un spline cúbico

                # refinamiento fino
                xi_fine = np.linspace(xi_sub[0], xi_sub[-1], 1000000)
                theta_fine = cs(xi_fine)

                # encontrar raíz
                idx = np.where(theta_fine <= 0)[0][0] #np where devuelve los índices donde se cumple la condición, [0][0] para obtener el primer índice donde theta_fine es menor o igual a cero
                self.xi1 = xi_fine[idx] #xi1 es el valor de xi donde theta se anula, es decir, el radio de la estrella polytropica

                break # salir del loop principal una vez encontrada la raíz

        # protección NaN
            if np.isnan(y_new[0]):
                print("NaN detectado, deteniendo")
                break

            y = y_new

        self.xi = np.array(xi_vals)
        self.theta = np.array(theta_vals)
        self.dtheta = np.array(dtheta_vals)

    # ============================================================
    # VERIFICACIÓN (FINITE DIFFERENCES)
    # ============================================================

    def finite_difference(self):
        """
        Calcula derivada usando diferencias finitas
        """

        h = 0.001
        theta = self.theta

        dtheta_fd = []

        for i in range(1, len(theta)-1): #no se puede calcular la derivada en los extremos con esta fórmula, por eso se itera desde 1 hasta len(theta)-1
            deriv = (theta[i+1] - theta[i-1]) / (2*h)
            dtheta_fd.append(deriv)
        return np.array(dtheta_fd)

