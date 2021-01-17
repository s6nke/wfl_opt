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
import os
from  wfl_opt_data import *         # data file; containing all relevant data; place in ./
import docplex.mp.model as cpx      # docplex requiring cplex installation
import matplotlib.pyplot as plt     # pyplot required for plotting
import numpy as np                  # numpy package

# ---------
# VARIABLES AND LISTS
# ---------
xAxis = 20           # size of x axis in m
yAxis = 20           # size of y axis in m
Nmin  = 5            # minimal number of turbines
Nmax  = 20           # maximal number of turbines
Dmin  = 3            # minimum distance between turbines
k       = 0.05       # wake decay coefficient
D       = 2          # rotor diameter

# -----------
# Set the optimization environment
OPenv = layout_optimization(xAxis, yAxis, 1, Dmin, k, D)



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
OP.add_constraint(OP.sum(x_vars[i] for i in OPenv.setV) <= Nmax, ctname="minimum number constraint")
OP.add_constraint(OP.sum(x_vars[i] for i in OPenv.setV) >= Nmin, ctname="maximum number constraint")
print("Number of turbines constraints set!")

# geometric constraint
# requires the set_E[i] to be defined as the set of nodes that are in the distance
# of Dmin around the chosen node [i]
# this adds a constraint for from every node [i] to every node[j] for all i,j el of grid
for i in OPenv.setV:
    for j in OPenv.setE[i]:
        OP.add_constraint(x_vars[i] + x_vars[j] <= 1, ctname="distance constraint")
print("Minimal distance constraint set!")

bigM = []
for i in OPenv.setV:
    Mi = 0
    for j in OPenv.not_setE[i]:
        xj, yj = OPenv.grid[j]
        Mi += OPenv.infer_matrix[i,yj,xj]
    bigM.append(Mi)

for i in OPenv.setV:
    inf_ij = 0
    for j in OPenv.setV:
        if j != i:
            xj, yj = OPenv.grid[j]
            inf_ij += OPenv.infer_matrix[i,yj,xj]*x_vars[j]
    OP.add_constraint(inf_ij <= (w_vars[i] + bigM[i]*(1-x_vars[i])), ctname="big M constraint")
print("Big-m constraint set!")

# wi positive constraint
for i in OPenv.setV:
    OP.add_constraint(w_vars[i] >= 0)

# ---------
# OBJECTIVE
# ---------

# objective funtion
# P() is the (linear) power function of the turbine depending on the wind speed
ppMWh = 300

obj = OP.sum((OPenv.Pi*x_vars[i] - w_vars[i])#*ppMWh
            #- OPenv.dist_matrix[ OPenv.grid[i][0], OPenv.grid[i][1] ]*x_vars[i]
            #+ OPenv.geo_matrix[  OPenv.grid[i][0], OPenv.grid[i][1] ]*x_vars[i]
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

print("---------------------------------")
print("Start Optimizer: \n\n")
OP.solve()

print("Number of turbines built: ", sum(x_vars[i].solution_value for i in OPenv.setV))
OPenv.sol = [OP.solution.get_value(x_vars[element]) for element in OPenv.setV]

OPenv.sol_obj = OP.solution.get_objective_value()
print(OPenv.sol_obj)

sol_indices = []
cnt = 0
for element in OPenv.setV:
    if OP.solution.get_value(x_vars[element]) == 1:
        sol_indices.append(cnt)
    cnt += 1
OPenv.sol_indices = sol_indices

for i in OPenv.setV:
    if OP.solution.get_value(x_vars[i]) == 1:
        print(OP.solution.get_value(w_vars[i]))

# ###########
# SAVE RESULT
# ###########
OPenv.sol_interference()
OPenv.save_sol()