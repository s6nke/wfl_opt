# import cplex as cpx
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

# INIT
OP = cpx.Model(name="Wind Farm Layout", log_output=True)

# define decision variables
x_vars = {(i): OP.binary_var(name="x_{0}".format(i)) for i in set_V}

# add constraints
OP.add_constraint(OP.sum(x_vars[i] for i in set_V) <= Nmax)
OP.add_constraint(OP.sum(x_vars[i] for i in set_V) >= Nmin)
constraint3 = {i: OP.add_constraint(OP.sum(x_vars[i] + x_vars[j] for j in set_E[i]) <= 1)}

print(set_E[0])
# create objective
obj = OP.sum(P(i)*x_vars[i] for i in set_V)

# set sense of the obj func
OP.maximize(obj)

# solve OP
print(OP.solve())

# plotting the solution
fig, ax = plt.subplots()
grid_points = plot_grid(grid)

for i in range(0,n):
    if OP.solution.get_value(x_vars[i]) == 1.0:
        plt.scatter(grid[i][0], grid[i][1], s=100, c="red", marker="*")
        ax.add_artist(plt.Circle((grid[i][0], grid[i][1]), Dmin, fill=False))

plt.show()
