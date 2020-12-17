# ###########################################
# OPTIMIZING THE CABLE ROUTING OF A WIND FARM
# ###########################################
#
# TODO:
# - adapt 6.6 to cnt j
# - check if substation node is not equal to turbine node
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
xAxis = 5           # size of x axis in m
yAxis = 5           # size of y axis in m
V0 = [5,yAxis*xAxis-2]
Ph = 1



# -----------
# Set the optimization environment
# loads solution file automatically
CRenv = cable_routing_optimization(xAxis, yAxis,V0)

CRenv.calc_setA()
CRenv.calc_cost()

print(CRenv.setVT0)
print(CRenv.setA)
##############
# CPLEX SOLVER
##############
# ----
# INITIALISE MODEL
# ----
# docplex instance with name and force output after solving
OP = cpx.Model(name="Cable Routing", log_output=True)

# ----
# VARIABLES
# ----
x_vars = OP.binary_var_matrix(keys1=CRenv.setVT0, keys2=CRenv.setVT0, name="x_%i_%i")
f_vars = OP.continuous_var_matrix(keys1=CRenv.setVT0, keys2=CRenv.setVT0, name="f_%f_%f")
# y_vars = {(i)} if more than one cable type


# ----
# CONSTRAINTS
# ----
print("Adding constraints!")
# (6.3) flow conservation constraint
print("(6.4)")
for h in CRenv.setVT:
    temp_c3 = []
    for i in CRenv.setVT0:
        if i!=h:
            temp_c3.append(f_vars[h,i] - f_vars[i,h])
    OP.add_constraint(OP.sum(temp_c3) == Ph,"c3_{0}".format(h))

print("(6.5)")
# (6.5) only one cable can exit a turbine
for h in CRenv.setVT:
    temp_c4 = []
    for j in CRenv.setVT0:
        if h!=j:
            temp_c4.append(x_vars[h,j])
    OP.add_constraint(OP.sum(temp_c4) == 1, "c5_h{0}".format(h))

print("(6.6)")
# (6.6) no cable can exit a substation
for h in CRenv.setV0:
    temp_c6 = []
    for j in CRenv.setVT0:
        if j!=h:
            temp_c6.append(x_vars[h,j])
    OP.add_constraint(OP.sum(temp_c6) == 0, "c6_h{0}".format(h))


print("(6.13)")
# (6.13) 
for arc in CRenv.setA:
    OP.add_constraint(f_vars[arc[0],arc[1]] >=0, "c13_{0}_{1}".format(arc[0],arc[1]))


# ----
# OBJECTIVE
# ----
print("Create objective")
obj = []
for arc in CRenv.setA:
    obj.append(CRenv.cost(arc)*x_vars[arc[0],arc[1]])
OP.minimize(OP.sum(obj))
OP.solve()

sol_arcs = []
for i in CRenv.setVT0:
    for j in CRenv.setVT0:
        if OP.solution.get_value(x_vars[i,j]) == 1:
            sol_arcs.append([i,j])
CRenv.sol_arcs = sol_arcs


##########
# PLOTTING
##########
fig,ax = plt.subplots()
CRenv.plot_turbines(CRenv.layout_sol, col="red")
CRenv.plot_substations()
CRenv.plot_arcs()
plt.show()