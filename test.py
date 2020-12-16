# ------------
# DEPENDENCIES
# ------------
from  wfl_opt_data import *         # data file; containing all relevant data; place in ./
import docplex.mp.model as cpx      # docplex requiring cplex installation
import matplotlib.pyplot as plt     # pyplot required for plotting
import numpy as np                  # numpy package

# ---------
# VARIABLES AND LISTS
# ---------
xAxis = 10           # size of x axis in m
yAxis = 10           # size of y axis in m
Nmin  = 0            # minimal number of turbines
Nmax  = 10           # maximal number of turbines
Dmin  = 0            # minimum distance between turbines
k       = 0.05       # wake decay coefficient
V       = 15         # m/s
D       = 1          # rotor diameter
v_wind  = 15
w_dir  = np.transpose(np.array([[1,1]]))

# -----------
# Set the optimization environment
OPenv = wfl_environmet(xAxis, yAxis, Dmin, k, V, D, v_wind, w_dir)
OPenv.find_initial_sol(Nmin, Nmax)

sol = [0]*OPenv.n
sol[0:1] = [1.0]*2
fig,ax = plt.subplots()
OPenv.plot_turbines(OPenv.initial_sol)
OPenv.plot_interference(OPenv.initial_sol)
plt.show()