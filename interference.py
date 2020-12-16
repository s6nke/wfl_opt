# ------------
# DEPENDENCIES
# ------------
import numpy as np
import matplotlib.pyplot as plt
from wfl_opt_data import *
# ---------
# VARIABLES
# ---------
k       = 0.05  # wake decay coefficient
V       = 15    # m/s
D       = 1     # rotor diameter
v_wind  = 15
xAxis = 10           # size of x axis in m
yAxis = 10           # size of y axis in m
Dmin  = 2            # minimum distance between turbines

w_dir  = np.transpose(np.array([1,1]))

# Init the setting of the opt
OP_env = wfl_environmet(xAxis, yAxis, Dmin)
interference_matrix = OP_env.calc_interference(D, w_dir, v_wind, V, k)

