#-------------------------------------
# Data file for optimizing a offshore wind farm layout
#-------------------------------------
#
# TODO:
# - adapt power function to represent nonlinear behaviour depending on wind speed
# - interference calculation or integration from different file (precalculation)
#
# ------------
# DEPENDENCIES
# ------------
import random
import numpy as np
import matplotlib.pyplot as plt
#
# ---------
# VARIABLES AND LISTS
# ---------
n     = 100          # number of nodes
xAxis = 10           # size of x axis in m
yAxis = 10           # size of y axis in m
grid  = []           # empty grid knodes
set_V = range(0,n)   # enumerate decsion variables
set_E = []           # store not allowed area for each node
Nmin  = 1            # minimal number of turbines
Nmax  = 50           # maximal number of turbines
Dmin  = 2            # minimum distance between turbines

# --------------------------
# POWER FUNCTION OF TURBINES
# --------------------------
# i.e. approximation of the nonlinear power delivery function of each turbine
# should be depending on the wind speed
def P(node):
    return 200

# -------------------
# ACCESSORY FUNCTIONS
# -------------------
# calculate kartesian difference of two grid nodes
def dist(node1, node2):
    # node1: first node of the grid; array-like
    # node2: second node of the grid; array-like
    dx = node1[0] - node2[0]
    dy = node1[1] - node2[1]
    return np.sqrt(dx**2 + dy**2)

# plot grid nodes
def plot_grid(grid):
    # grid: array of nodes; array of array-like
    for i in range(0,n):
        # scatter the node points
        plt.scatter(grid[i][0], grid[i][1], s=10, c="black")
        # add number to the node
        plt.text(grid[i][0], grid[i][1], "{0}".format(i))
    return None

# ---------------
# NODE GENERATION
# ---------------
# generate list of [x,y]-coordinates
for ii in range(0,xAxis):
    for jj in range(0,yAxis):
        # GRID
        grid.append([ii,jj])

# generate the set E_i
# elements of each subset are indexes of the set_V that are in
# a distance of Dmin around turbine i
for l in range(0,n):
    # create empty sub-list
    set_E.append([])
    for k in range(0,n):
        # check if the node is in Dmin to our node[l]
        if dist(grid[l], grid[k]) <= Dmin:
            # if too close add to set_E sublist
            set_E[l].append(k)
    # remove the node itself from the set to prevent false addition if set true
    set_E[l].remove(l)
