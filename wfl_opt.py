# #######################################
# OPTIMIZING AN OFFSHORE WIND PARK LAYOUT
# #######################################
#
# TODO:
# - add kpi of the number of turbines that are installed
# - interference matrix has to be included
# - cost of foundation has to be calculated included
# - adapt objective to represent interference and foundation cost
# - prettify the plots
#
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
xAxis = 20           # size of x axis in m
yAxis = 20           # size of y axis in m
Nmin  = 0            # minimal number of turbines
Nmax  = 200           # maximal number of turbines
Dmin  = 5            # minimum distance between turbines
k       = 0.05       # wake decay coefficient
V       = 15         # m/s
D       = 1          # rotor diameter
v_wind  = 15
w_dirs  = np.array([[1,0], [0,1], [1,1], [-1,0], [-1,-1], [0,-1]])

# -----------
# Set the optimization environment
OPenv = wfl_environmet(xAxis, yAxis, Dmin, k, V, D, v_wind, w_dirs)


##############
# CPLEX SOLVER
##############
# ----
# INITIALISE MODEL
# ----
# docplex instance with name and force output after solving
OP = cpx.Model(name="Wind Farm Layout", log_output=True)

# -----
# DVARS
# -----
# create a number of binary decsision variables based on th set_V (from data file)
# "x_{0}".format defines the name
x_vars = {(i): OP.binary_var(name="x_{0}".format(i)) for i in OPenv.setV}

# -----------
# CONSTRAINTS
# -----------
# min max constraint of number of turbines
OP.add_constraint(OP.sum(x_vars[i] for i in OPenv.setV) <= Nmax)
OP.add_constraint(OP.sum(x_vars[i] for i in OPenv.setV) >= Nmin)

# geometric constraint
# requires the set_E[i] to be defined as the set of nodes that are in the distance
# of Dmin around the chosen node [i]
# this adds a constraint for from every node [i] to every node[j] for all i,j el of grid
for cnt in range(0,OPenv.n):
    for el in OPenv.setE[cnt]:
        OP.add_constraint(x_vars[cnt] + x_vars[el] <= 1)


# ---------
# OBJECTIVE
# ---------
def w(i, ifm, grid):
    # -----
    # calculate the interference caused by site i to all sites j
    #
    # i:    integer of the current node in the set_V
    # ifm:  array-like, precalculated interference matrix
    # ----- 
    ifm_wi = 0
    for j in OPenv.setV:
        if j != i:
            ifm_wi = ifm_wi + ifm[i][ grid[j][0] ][ grid[j][1] ]*x_vars[j]
    return ifm_wi


# objective funtion
# P() is the (linear) power function of the turbine depending on the wind speed
obj =   OP.sum(OPenv.P(15)*x_vars[i] - 
        w(i, OPenv.infer_matrix, OPenv.grid)*x_vars[i]
        for i in OPenv.setV)

# -----
# SENSE
# -----
# the sense can be minimize or maximize;
# we want our power to be maximized
OP.maximize(obj)

# ----
# SOLVE
# ----
# solve OP, output is forced in the initialization, otherwise use print()
OPenv.warmstart(OP, x_vars, Nmin, Nmax)

OP.solve()

print("Number of turbines built: ", sum(x_vars[i].solution_value for i in OPenv.setV))
OPenv.sol = [OP.solution.get_value(x_vars[element]) for element in OPenv.setV]

print(OP.solution.objective_value)
# ########
# PLOTTING
# ########
fig, ax = plt.subplots()            # create new figure
OPenv.plot_grid(numbers=False)       # plot_grid() in data file; plots grid nodes
OPenv.plot_turbines(OPenv.sol)
OPenv.plot_turbines(OPenv.initial_sol)
OPenv.plot_interference(OPenv.sol)

# show the plot, can be removed to prevent the pop-up figure
plt.show()

# ###########
# SAVE RESULT
# ###########
OPenv.save_sol()