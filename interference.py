# ------------
# DEPENDENCIES
# ------------
import numpy as np
import matplotlib.pyplot as plt

# ---------
# VARIABLES
# ---------
grid    = []    # empty grid list -> np.array
xAxis   = 10    # x in m
yAxis   = 10    # y in m
k       = 0.05  # wake decay coefficient
V       = 15    # m/s
D       = 1     # rotor diameter

# ---------
# FUNCTIONS
# ---------
def dist(node1, node2):
    return np.linalg.norm(node1-node2)

def jensen_mod(upwind_stream, rotor_diam, distance):
    k = 0.05 # wake decay coefficient
    return upwind_stream*(1 - np.sqrt(1 - thrust_coeff(upwind_stream)))*(rotor_diam/(rotor_diam+2*k*distance))**2

def thrust_coeff(wind_speed):
    mod_coeff = 6 # chosen arbitrary to get quite good values
    if wind_speed <= 3:
        return 0
    else:
        return mod_coeff/(wind_speed)

def plot_grid(grid):
    for i in range(0,grid.shape[0]-1):
        plt.scatter(grid[i][0], grid[i][1], s=1, c="black")
        #plt.text(grid[i][0], grid[i][1], "{0}".format(i))
    return None

# ----
# MAIN
# ----
# node generation
for ii in range(0,xAxis):
    for jj in range(0,yAxis):
        grid.append([ii,jj])
grid = np.array(grid)

interferenz_matrix = np.zeros(shape=(xAxis*xAxis, xAxis, yAxis))
w_dir  = np.transpose(np.array([1,1]))

for i in range(0, grid.shape[0]):
    infer  = []
    for j in range(0, grid.shape[0]):
        alpha   = np.matmul(np.transpose(grid[j] - grid[i]), w_dir) / np.linalg.norm(w_dir)**2
        P_prime = grid[i] + alpha*w_dir
        delta_1 = dist(grid[i], P_prime)
        delta_2 = dist(grid[j], P_prime)
        if delta_2 <= (D+2*k*delta_1)/2:
            infer.append(jensen_mod(V, D, delta_1))
        else:
            infer.append(0)
    interferenz_matrix[i] = np.array(infer).reshape((xAxis, yAxis))
