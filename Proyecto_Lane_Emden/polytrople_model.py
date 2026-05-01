from dbm import error

import numpy as np
from scipy.interpolate import CubicSpline


class PolytropeModel: #Parametros físicos de la estrella polytropica a partir de la solución de Lane-Emden
    """
    Modelo físico basado en la solución de Lane-Emden
    """

    def __init__(self, solver):

        # Verificación básica
        if solver.xi is None or solver.xi1 is None:
            raise ValueError("Debes ejecutar solver.solve() primero")

        self.solver = solver
        self.n = solver.n
        self.xi1 = solver.xi1
        self.xi = solver.xi
        self.theta = solver.theta

        self.dtheta_xi1 = None

    # ============================================================
    # DERIVADA EN LA SUPERFICIE
    # ============================================================

    def compute_surface_gradient(self):
        """
        Calcula dtheta(chi_1) usando spline
        """

        xi_sub = self.xi[-5:]
        theta_sub = self.theta[-5:]

        cs = CubicSpline(xi_sub, theta_sub)

        self.dtheta_xi1 = cs(self.xi1, 1)

        return self.dtheta_xi1

    # ============================================================
    # RELACIÓN DE DENSIDADES
    # ============================================================

    def density_ratio(self):
        """
        pc / p
        """

        if self.dtheta_xi1 is None:
            self.compute_surface_gradient()

        return -self.xi1 / (3 * self.dtheta_xi1)
    
    def central_density_physical(xi1, dtheta_xi1):
        """
        Calcula densidad central física (g/cm^3)
        para un politropo usando M☉ y R☉
        """

        G = 6.674e-8
        M_sun = 1.989e33      # g
        R_sun = 6.957e10      # cm

        rho_c = (M_sun / (4 * np.pi * R_sun**3)) * (xi1 / (-dtheta_xi1))

        return rho_c

    # ============================================================
    # MASA ADIMENSIONAL
    # ============================================================

    def mass_parameter(self):
        """
        def in pdf
        """

        if self.dtheta_xi1 is None:
            self.compute_surface_gradient()

        return -self.xi1**2 * self.dtheta_xi1

    # ============================================================
    # CONSTANTE ω_n
    # ============================================================

    def omega_n(self):
        """
        def in pdf
        """

        if self.dtheta_xi1 is None:
            self.compute_surface_gradient()

        return -self.xi1**((self.n + 1)/(self.n - 1)) * self.dtheta_xi1

    # ============================================================
    # CONSTANTE N_n
    # ============================================================

    def N_n(self):
        """
        def in pdf
        """

        omega = self.omega_n()

        return (1 / (self.n + 1)) * (4 * np.pi / (omega**(self.n - 1)))**(1/self.n)

    # ============================================================
    # CONSTANTE W_n
    # ============================================================

    def W_n(self):
        """
        def in pdf
        """

        mp = self.mass_parameter()

        return (1 / (4 * np.pi * (self.n + 1))) * (mp**2)

    # ============================================================
    # PRESIÓN CENTRAL
    # ============================================================

    def central_pressure(self):
        """
        Pc en dyn/cm²
        """

        G = 6.674e-8
        M_sun = 1.989e33
        R_sun = 6.957e10

        Wn = self.W_n()
        Pc = Wn * (G * M_sun**2 / R_sun**4)

        return Pc / 1e17
    def derivate_check(self):
        """
        Verificación RK4 vs diferencias finitas
        """

        dtheta_fd = self.solver.finite_difference()
        dtheta_rk = self.solver.dtheta[1:-1]

        error = np.abs(dtheta_rk - dtheta_fd)

        return {
            "error_mean": np.mean(error),
            "error_max": np.max(error)
            }

    # ============================================================
    # RESUMEN 
    # ============================================================

    def summary(self):
        """
        Devuelve todos los parámetros en un diccionario
        """

        return {
            "n": self.n,
            "xi1": self.xi1,
            "dtheta_xi1": self.compute_surface_gradient(),
            "rho_ratio": self.density_ratio(),
            "mass_param": self.mass_parameter(),
            "omega_n": self.omega_n(),
            "N_n": self.N_n(),
            "W_n": self.W_n(),
            "Pc": self.central_pressure(),
            "error_mean": self.derivate_check()["error_mean"],
            "error_max": self.derivate_check()["error_max"]
        }