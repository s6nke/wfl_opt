# ###########################################
# OPTIMIZING THE CABLE ROUTING OF A WIND FARM
# ###########################################
#
# TODO:
# - everything
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
Nmin  = 0            # minimal number of turbines
Nmax  = 200           # maximal number of turbines
Dmin  = 2            # minimum distance between turbines
k       = 0.05       # wake decay coefficient
V       = 15         # m/s
D       = 1          # rotor diameter
v_wind  = 15
w_dirs  = np.array([[1,1]])#, [0,1], [1,1], [-1,0], [-1,-1], [0,-1]])

V0      = [0,yAxis*xAxis-1]
Ph = 1




# -----------
# Set the optimization environment
# loads solution file automatically
CRenv = cable_routing_optimization(xAxis, yAxis,V0)

CRenv.calc_setA()
CRenv.calc_cost()

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
x_vars = OP.binary_var_matrix(keys1=CRenv.setV, keys2=CRenv.setV, name="x_%i_%i")
f_vars = OP.continuous_var_matrix(keys1=CRenv.setV, keys2=CRenv.setV, name="f_%f_%f")
# y_vars = {(i)} if more than one cable type


# ----
# CONSTRAINTS
# ----
print("Adding constraints!")
# (6.4) flow conservation constraint
print("(6.4)")
for h in np.unique(np.append(CRenv.setVT,CRenv.setV0)):
    for i in CRenv.setV:
        if i!=h:
            arcHI = CRenv.setA[h][i]
            arcIH = CRenv.setA[i][h]
            OP.add_constraint(f_vars[arcHI[0],arcHI[1]] - f_vars[arcIH[0],arcIH[1]] == Ph,"c4_{0}_{1}".format(arcIH[0],arcIH[1]))

print("(6.5)")
# (6.5) only one cable can exit a turbine
for h in CRenv.setVT:
    temp_c4 = []
    for node in CRenv.setA[h]:
        temp_c4.append(x_vars[node[0], node[1]])
    OP.add_constraint(OP.sum(temp_c4) == 1, "c5_h{0}".format(h))

print("(6.6)")
# (6.6) no cable can exit a substation
for h in CRenv.setV0:
    temp_c6 = []
    for node in CRenv.setA[h]:
        temp_c6.append(x_vars[node[0], node[1]])
    OP.add_constraint(OP.sum(temp_c6) == 0, "c6_h{0}".format(h))


print("(6.13)")
# (6.13) 
for i in CRenv.setV:
    for arc in CRenv.setA:
        for el in arc:
            OP.add_constraint(f_vars[el[0],el[1]] >=0, "c13_{0}_{1}".format(el[0],el[1]))


print("Create objective")
# opjective
obj = []
for i in CRenv.setV:
    for arcs in CRenv.setA:
            for arc in arcs:
                obj.append(CRenv.costs(arc) * x_vars[arc[0],arc[1]])

OP.minimize(OP.sum(obj))
OP.solve()

print(OP.solution)

##########
# PLOTTING
##########
fig,ax = plt.subplots()
CRenv.plot_turbines(CRenv.layout_sol, col="red")
plt.show()