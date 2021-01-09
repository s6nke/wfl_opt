import numpy as np
import matplotlib.pyplot as plt
from wfl_opt_data import *
import os

axx = 20
axy = 20
direc = os.path.dirname(__file__)
fig_dir = os.path.join(direc, "figures")

Env = wf_environment(axx,axy)
Env.load_layout_sol()
Env.load_cable_sol()
Env.load_dist_geo()
# - layout_sol
# - layout_sol_index
# - inf_sol

"""
Plot interference of the solution
"""
fig0, ax0 = plt.subplots()
heat = plt.imshow(Env.inf_sol, cmap='jet', interpolation='bilinear')
plt.colorbar(heat)
plt.close(fig0)

"""
Plot solution layout
"""
fig1, ax1 = plt.subplots()
Env.plot_turbines(ax1)
plt.close(fig1)

"""
Plot cables
"""
fig2, ax2 = plt.subplots()
Env.plot_arcs()
plt.close(fig2)

"""
Plot distance gradient and geometry of depth
"""
fig3, ax3 = plt.subplots(2)
x = np.arange(0,axx)
y = np.arange(0,axy)
xm, ym = np.meshgrid(x,y)
surf = ax3[0].contourf(xm,ym,Env.geo_matrix, cmap="hot")
fig3.colorbar(surf, ax=ax3[0])
surf = ax3[1].contourf(xm,ym,Env.geo_matrix, cmap="hot")
fig3.colorbar(surf, ax=ax3[1])
Env.plot_turbines(ax3[1])


"""
Plot some example interference
"""
#LO = layout_optimization(20,20,w_dirs=np.array([[1/np.sqrt(2), 1/np.sqrt(2)]]))
#fig4, ax4 = plt.subplots()
#heat1 = plt.imshow(LO.infer_matrix[0], cmap='jet', interpolation='bilinear')
#plt.colorbar(heat1)
#plt.savefig(os.path.join(fig_dir, "ex_inter.png"))
#plt.close(fig4)


# show plots
plt.show()