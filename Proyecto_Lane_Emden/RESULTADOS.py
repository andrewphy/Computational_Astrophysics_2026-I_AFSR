import matplotlib.pyplot as plt   
import numpy as np 
import pandas as pd
from calibration import start
from SOLVER_LANE_EQ import LaneEmdenSolver
from polytrople_model import PolytropeModel
from scipy.interpolate import interp1d
#PRUEBA CON n=1 raiz= pi 
calibration = start() 
n_values = [3.0, 3.15, 3.30, 3.45, 3.60, 3.75]
results=[]

plt.figure(figsize=(8,6))

for n in n_values:
    print(f'--Calculando n={n}')
    solver = LaneEmdenSolver(n)
    solver.solve()
    # grafica
    plt.plot(solver.xi,solver.theta,label=rf"$n={n},\ \xi_1\approx {solver.xi1:.3f}$")
    plt.scatter(solver.xi1, 0, zorder=3)
    # modelo
    model=PolytropeModel(solver)
    #tabla
    data=model.summary()
    results.append(data)
# ============================================================
# CONFIG PLOT
# ============================================================
plt.axhline(0, linestyle='--')# línea en θ = 0 (superficie)
plt.xlabel(r"$\xi$")
plt.ylabel(r"$\theta(\xi)$")
plt.title("Solución de la ecuación de Lane-Emden")
plt.legend()
plt.grid(alpha=0.4)
plt.minorticks_on()
plt.tight_layout()
#plt.savefig("lane_emden_plot.pdf", dpi=300)
plt.show()

# ============================================================
# TABLA
# ============================================================
df = pd.DataFrame(results)
df = df.round({"xi1": 4,"dtheta_xi1": 4,"rho_ratio": 4,"mass_param": 4,"omega_n": 4,"N_n": 4,"W_n": 4,"Pc": 4})
df = df[["n","xi1","dtheta_xi1","rho_ratio","mass_param","omega_n","N_n","W_n","Pc","error_mean","error_max"]]
print("\nTabla de resultados:\n")
print(df)
# guardar en CSV
#df.to_csv("/home/pipe-desktop/Computational_Astrophysics_2026-I_AFSR/PROYECT_LAN/tabla_lane_emden.csv", index=False)

# ============================================================
# Comaparacion con datos 
# ============================================================
    ## Datos del modelo solar 
    ##===========================================================
#url = "https://users-phys.au.dk/jcd/solar_models/cptrho.l5bi.d.15c"
solar = pd.read_csv(
    "/home/pipe-desktop/Computational_Astrophysics_2026-I_AFSR/PROYECT_LAN/data_solar_real.csv",
    sep=r"\s+",
    comment="#",
    names=["r", "c", "rho", "p", "gamma1", "T"])

#print(solar.head())
r_solar = solar["r"].values
r_solar=r_solar[::-1]#invertir para que vaya de 0 a 1
rho_solar = solar["rho"].values
rho_solar=rho_solar[::-1]
rho_centra=rho_solar[0] #densidad central
rho_solar_normalized = rho_solar / rho_centra #densidad normalizada


    ##===========================================================
    ## Datos de mi modelo n=3.0
    ##===========================================================
solver_n3 = LaneEmdenSolver(3.0)
solver_n3.solve()
poly_n3 = PolytropeModel(solver_n3)
xi = solver_n3.xi
theta = solver_n3.theta
xi1 = solver_n3.xi1
dxi1 = poly_n3.compute_surface_gradient()
rho_central_model = PolytropeModel.central_density_physical(xi1, dxi1)
# radio normalizado
r_model = xi / xi1

# densidad 
rho_model = rho_central_model * theta**3

interp = interp1d(r_model, rho_model, kind='cubic', fill_value="extrapolate")
rho_model_interp = interp(r_solar)
#print(rho_model_interp)
    ##===========================================================
    ## PLOT COMPARACIÓN
    ##===========================================================
plt.figure(figsize=(8,6))
plt.plot(r_solar, rho_solar, label="Modelo Solar (Datos)", color="orange")
plt.plot(r_solar,rho_model_interp, label="Modelo Polytropic n=3 (Interpolado)", color="green", linestyle='--')
plt.plot(r_model, rho_model, label="Modelo Polytropic n=3", color="blue", alpha=0.4)
plt.axvspan(0, 0.2, color='red', alpha=0.1)
plt.text(0.05, max(rho_solar), "NUCLEO", color="red")
plt.xlabel("Radio Normalizado (r/R)")
plt.ylabel(fr"Densidad ($\rho \, [g/cm^3]$)")
plt.title("Comparación de Densidad: Modelo Solar vs Polytropic n=3")
plt.legend()
plt.grid(alpha=0.4)
plt.minorticks_on()
plt.tight_layout()
#plt.savefig("comparacion_densidad_n3_datosrealesvsmodelo.pdf", dpi=300)
plt.show()  