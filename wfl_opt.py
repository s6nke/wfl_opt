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
Dmin  = 2            # minimum distance between turbines
k       = 0.05       # wake decay coefficient
V       = 15         # m/s
D       = 1          # rotor diameter
v_wind  = 15
w_dirs  = np.array([[1,1]])#, [0,1], [1,1], [-1,0], [-1,-1], [0,-1]])

# -----------
# Set the optimization environment
OPenv = layout_optimization(xAxis, yAxis, Dmin, k, V, D, v_wind, w_dirs)


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
w_vars = {(i): OP.continuous_var(name="w_{0}".format(i)) for i in OPenv.setV}
print("Decision variables set!")

# -----------
# CONSTRAINTS
# -----------
# min max constraint of number of turbines
OP.add_constraint(OP.sum(x_vars[i] for i in OPenv.setV) <= Nmax)
OP.add_constraint(OP.sum(x_vars[i] for i in OPenv.setV) >= Nmin)
print("Number of turbines constraints set!")

# geometric constraint
# requires the set_E[i] to be defined as the set of nodes that are in the distance
# of Dmin around the chosen node [i]
# this adds a constraint for from every node [i] to every node[j] for all i,j el of grid
for cnt in range(0,OPenv.n):
    for el in OPenv.setE[cnt]:
        OP.add_constraint(x_vars[cnt] + x_vars[el] <= 1)
print("Minimal distance constraint set!")


bigM = []
for i in OPenv.setV:
    Mi = 0
    for j in OPenv.not_setE[i]:
        xj, yj = OPenv.grid[j]
        Mi = Mi + OPenv.infer_matrix[i,xj,yj]
    bigM.append(Mi)

for i in OPenv.setV:
    for j in OPenv.setV:
        if i != j:
            # look up coordinates
            xj, yj = OPenv.grid[j]
            # add constraint
            if OPenv.infer_matrix[i,xj,yj] < bigM[i]:
                OP.add_constraint(OPenv.infer_matrix[i,xj,yj] <= w_vars[i] + bigM[i]*(1-x_vars[i]))
print("Big-m constraint set!")

# wi positive constraint
for i in OPenv.setV:
    OP.add_constraint(w_vars[i] >= 0)

# ---------
# OBJECTIVE
# ---------

# objective funtion
# P() is the (linear) power function of the turbine depending on the wind speed
obj = OP.sum(OPenv.P(15)*x_vars[i] - w_vars[i] for i in OPenv.setV)
obj += OP.sum( OPenv.dist_matrix[ OPenv.grid[i][0], OPenv.grid[i][1] ]*1/10*x_vars[i] for i in OPenv.setV)
obj += OP.sum( OPenv.geo_matrix[ OPenv.grid[i][0], OPenv.grid[i][1] ]*x_vars[i] for i in OPenv.setV)

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

print("---------------------------------")
print("Start Optimizer: \n\n")
OP.solve()

print("Number of turbines built: ", sum(x_vars[i].solution_value for i in OPenv.setV))
OPenv.sol = [OP.solution.get_value(x_vars[element]) for element in OPenv.setV]


sol_indices = []
cnt = 0
for element in OPenv.setV:
    if OP.solution.get_value(x_vars[element]) == 1:
        sol_indices.append(cnt)
    cnt += 1
OPenv.sol_indices = sol_indices


# ###########
# SAVE RESULT
# ###########
OPenv.sol_interference()
OPenv.save_sol()