import numpy as np
def start():#PRUEBA CON n=1 raiz= pi 
    from SOLVER_LANE_EQ import LaneEmdenSolver
    solver = LaneEmdenSolver(n=1.0)
    solver.solve()
    error = 100*abs(solver.xi1 - np.pi)/np.pi
    print("xi end =", solver.xi[-1])
    print("xi1 numérico =", solver.xi1)
    print("xi1 exacto   =", np.pi)
    print("error %      =", error)
