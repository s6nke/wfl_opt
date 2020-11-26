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

#---------
# IMPORT FROM WFL_OPT_DATA FILE:
# n
# grid
# set_V
# Nmin
# Nmax
# Dmin
# --------

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
x_vars = {(i): OP.binary_var(name="x_{0}".format(i)) for i in set_V}

# -----------
# CONSTRAINTS
# -----------
# min max constraint of number of turbines
OP.add_constraint(OP.sum(x_vars[i] for i in set_V) <= Nmax)
OP.add_constraint(OP.sum(x_vars[i] for i in set_V) >= Nmin)

# geometric constraint
# requires the set_E[i] to be defined as the set of nodes that are in the distance
# of Dmin around the chosen node [i]
# this adds a constraint for from every node [i] to every node[j] for all i,j el of grid
for cnt in range(0,n):
    for el in set_E[cnt]:
        OP.add_constraint(x_vars[cnt] + x_vars[el] <= 1)

# ---------
# OBJECTIVE
# ---------
# simplified version of the objective
# adaptable to different scenarios
# P() is the power output function of the wind turbine
obj = OP.sum(P(i)*x_vars[i] for i in set_V)

# -----
# SENSE
# -----
# the sense can be minimize or maximize;
# we want our power to be maximized
OP.maximize(obj)

# ----
# SOLVE
# ----
# solve OP, output is forced in the initialization, outherwise use print()
OP.solve()

# --------
# PLOTTING
# --------
fig, ax = plt.subplots()            # create new figure
grid_points = plot_grid(grid)       # plot_grid() in data file; plots grid nodes
for i in range(0,n):
    # plot a red star on top of the grid node if the decsision variable is set 1
    # also adds a Dmin-diameter circle around this point to show forbidden zones
    if OP.solution.get_value(x_vars[i]) == 1.0:
        plt.scatter(grid[i][0], grid[i][1], s=100, c="red", marker="*")
        ax.add_artist(plt.Circle((grid[i][0], grid[i][1]), Dmin, alpha=0.1))

# show the plot, can be removed to prevent the pop-up figure
plt.show()
