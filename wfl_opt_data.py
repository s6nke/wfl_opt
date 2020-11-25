#-------------------------------------
# Data file for optimizing a offshore wind farm layout
#-------------------------------------

# ------------
# DEPENDENCIES
# ------------
import random
import numpy as np
import matplotlib.pyplot as plt

# ---------
# VARIABLES AND LISTS
# ---------
n     = 100          # number of nodes
xAxis = 10           # size of x axis in m
yAxis = 10           # size of y axis in m
grid  = []          # empty grid knodes
set_V = range(0,n)  # enumerate decsion variables
set_E = []          # store not allower area for each node
Nmin  = 1           # minimal number of turbines
Nmax  = 50          # maximal number of turbines
Dmin  = 1           # minimum distance between turbines

# --------------------------
# POWER FUNCTION OF TURBINES
# --------------------------
# i.e. approximation of the nonlinear power delivery function of each turbine
# depending on the wind speed
def P(node):
    return 200

# -------------------
# ACCESSORY FUNCTIONS
# -------------------
# calculate kartesian difference to ensure minimal distance
def dist(node1, node2):
    dx = node1[0] - node2[0]
    dy = node1[1] - node2[1]
    return np.sqrt(dx**2 + dy**2)

# --------
# PLOTTING
# --------
def plot_grid(grid):
    for i in range(0,n):
        plt.scatter(grid[i][0], grid[i][1], s=10, c="black")
        plt.text(grid[i][0], grid[i][1], "{0}".format(i))
    return None


# ---------------
# NODE GENERATION
# ---------------
# generate list of (x,y)-coordinates
for ii in range(0,xAxis):
    for jj in range(0,yAxis):
        grid.append([ii,jj])

# generate the set E_i
# elements of each subset are indexes of the set_V that are in
# a distance of Dmin around the index i
for l in range(0,n):
    set_E.append([])
    for k in range(0,n):
        if dist(grid[l], grid[k]) <= Dmin:
            set_E[l].append(k)
    set_E[l].remove(l)
