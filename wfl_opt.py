from  wfl_opt_data import *
import docplex.mp.model as cpx
import matplotlib.pyplot as plt

#---------
# IMPORT FROM DATA FILE:
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
# INIT
# ----
OP = cpx.Model(name="Wind Farm Layout", log_output=True)

# -----
# DVARS
# -----
x_vars = {(i): OP.binary_var(name="x_{0}".format(i)) for i in set_V}

# -----------
# CONSTRAINTS
# -----------
# min max constraint of number of turbines
OP.add_constraint(OP.sum(x_vars[i] for i in set_V) <= Nmax)
OP.add_constraint(OP.sum(x_vars[i] for i in set_V) >= Nmin)

# iteravily add constraints of Dmin
for cnt in range(0,n):
    for el in set_E[cnt]:
        OP.add_constraint(x_vars[cnt] + x_vars[el] <= 1)

# ---------
# OBJECTIVE
# ---------
obj = OP.sum(P(i)*x_vars[i] for i in set_V)

# -----
# SENSE
# -----
OP.maximize(obj)


##########
# MAIN
##########
# solve OP
OP.solve()

# --------
# PLOTTING
# --------
fig, ax = plt.subplots()
grid_points = plot_grid(grid)
for i in range(0,n):
    if OP.solution.get_value(x_vars[i]) == 1.0:
        plt.scatter(grid[i][0], grid[i][1], s=100, c="red", marker="*")
        ax.add_artist(plt.Circle((grid[i][0], grid[i][1]), Dmin, alpha=0.1))
plt.show()
