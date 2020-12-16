from  wfl_opt_data import *   
import numpy as np    

# init the wfl object
OPenv = wfl_environmet(20,20)
OPenv.load_sol()

# create the figure
fig, ax = plt.subplots()             
OPenv.plot_grid(numbers=False)       
OPenv.plot_turbines(OPenv.sol)
OPenv.plot_interference(OPenv.sol)

plt.show()