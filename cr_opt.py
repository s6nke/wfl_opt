# ###########################################
# OPTIMIZING THE CABLE ROUTING OF A WIND FARM
# ###########################################
#
# TODO:
# - check if substation node is not equal to turbine node
# - add different cable types
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
V0 = [5,yAxis*xAxis-2]  # nodes of substations
Ph = 1               # power coeff of turbine
k_cap = [2,10]           # capacity of the cable in one arc -> array when multiple cable types

# -----------
# Set the optimization environment
# loads solution file automatically if existant
CRenv = cable_routing_optimization(xAxis, yAxis,V0)

# calculate neccessary set of arcs
CRenv.calc_setA()


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
# decision variable if cable on the arc is built
x_vars = OP.binary_var_matrix(keys1=CRenv.setVT0, keys2=CRenv.setVT0, name="x_%i_%i")
# variable defining the power present on the arc
f_vars = OP.continuous_var_matrix(keys1=CRenv.setVT0, keys2=CRenv.setVT0, name="f_%i_%i")


# ----
# CONSTRAINTS
# ----
print("Adding constraints!")
print("(6.3)")
for h in CRenv.setVT:
    """ 
    The power outflow minus the power inflow has to coincide with
    the actual production of the node h
    """
    temp_c3 = []
    for i in CRenv.setVT0:
        if i!=h:
            temp_c3.append(f_vars[h,i] - f_vars[i,h])
    OP.add_constraint(OP.sum(temp_c3) == Ph,"c3_{0}".format(h))

print("(6.4)")
for arc in CRenv.setA:
    """
    upper bounds the power f that can be on the arc
    -> pushes the arcs to the sub stations; mandatory
    """
    OP.add_constraint(k_cap*x_vars[arc[0],arc[1]] >= f_vars[arc[0],arc[1]])

print("(6.5)")
for h in CRenv.setVT:
    """
    Ensures that only one cable can exit the turbine node;
    it is not possible to split the power exit
    """
    temp_c4 = []
    for j in CRenv.setVT0:
        if h!=j:
            temp_c4.append(x_vars[h,j])
    OP.add_constraint(OP.sum(temp_c4) == 1, "c5_h{0}".format(h))

# (6.6) no cable can exit a substation
print("(6.6)")
for h in CRenv.setV0:
    """
    Ensures that no arc can exit a substation; only entering is possible
    """
    temp_c6 = []
    for j in CRenv.setVT0:
        if j!=h:
            temp_c6.append(x_vars[h,j])
    OP.add_constraint(OP.sum(temp_c6) == 0, "c6_h{0}".format(h))

print("(6.13)")
for arc in CRenv.setA:
    """
    Enforces positivity of the power -> no negativ power
    """
    OP.add_constraint(f_vars[arc[0],arc[1]] >=0, "c13_{0}_{1}".format(arc[0],arc[1]))

# ----
# OBJECTIVE
# ----
print("Create objective")
obj = []
for arc in CRenv.setA:
    """
    simply cost times if-build;
    adapt to multiple cable types!
    """
    obj.append(CRenv.cost(arc)*x_vars[arc[0],arc[1]])
OP.minimize(OP.sum(obj))
OP.solve()

# save solution to CRenv and file
# findinf indices...
sol_arcs = []
for i in CRenv.setVT0:
    for j in CRenv.setVT0:
        if OP.solution.get_value(x_vars[i,j]) == 1:
            sol_arcs.append([i,j])
# save to env
CRenv.sol_arcs = sol_arcs
# save to file
CRenv.save_sol()

##########
# PLOTTING
##########
fig,ax = plt.subplots()
CRenv.plot_turbines(CRenv.layout_sol, col="red")
CRenv.plot_substations()
CRenv.plot_arcs()
plt.show()