import numpy as np
import matplotlib.pyplot as plt
from wfl_opt_data import *
import os

# Paramters
axx = 20
axy = 20
x = np.arange(0,axx)
y = np.arange(0,axy)
xm, ym = np.meshgrid(x,y)
direc = os.path.dirname(__file__)
fig_dir = os.path.join(direc, "figures")

# Load values from solution files
Env = wf_environment(axx,axy)
Env.load_layout_sol(1)
print(Env.layout_sol_obj)
Env.load_cable_sol()
Env.load_dist_geo()
Env.load_initial_sol()


"""
Plot interference of the solution
"""
fig0, ax0 = plt.subplots()
heat = plt.imshow(Env.inf_sol, cmap='jet', interpolation='bilinear')
cbar = plt.colorbar(heat)
cbar.ax.set_title("Leistungsverlust in kW")
Env.plot_turbines(ax0)
plt.xlabel("x")
plt.ylabel("y")
plt.gca().invert_yaxis()
#plt.savefig(os.path.join(fig_dir, "inf_sol_20_1.png"))
plt.close(fig0)

"""
Plot solution layout
"""
fig1, ax1 = plt.subplots()
Env.plot_turbines(ax1)
plt.xlabel("x")
plt.ylabel("y")
#plt.savefig(os.path.join(fig_dir, "turbines.png"   ))
plt.close(fig1)

"""
Plot cables
"""
fig2, ax2 = plt.subplots()
Env.plot_arcs()
Env.plot_turbines(ax2)
Env.plot_substations(ax2)
plt.xlabel("x")
plt.ylabel("y")
#plt.savefig(os.path.join(fig_dir, "cables2.png"))
#plt.close(fig2)

"""
Plot distance geometry of depth in 3D
"""
fig3a = plt.figure()
ax3a = fig3a.add_subplot(projection="3d")
surf = ax3a.plot_surface(xm,ym,Env.geo_matrix, cmap="hot")
ax3a.set_zlabel("Tiefe [m]")
plt.xlabel("x")
plt.ylabel("y")
#plt.savefig(os.path.join(fig_dir, "geo3D.png"))
plt.close(fig3a)

"""
Plot distance geometry of depth in 2D
"""
fig3a2, ax3a2 = plt.subplots()
surf = ax3a2.imshow(Env.geo_matrix, cmap="hot", interpolation="bilinear")
fig3a2.colorbar(surf, ax=ax3a2)
plt.xlabel("x")
plt.ylabel("y")
#plt.savefig(os.path.join(fig_dir, "geo2D.png"))
plt.close(fig3a2)

"""
Plot geometry of depth with solution
"""
fig3b, ax3b = plt.subplots()
surf = ax3b.imshow(Env.geo_matrix, cmap="hot", interpolation="bilinear")
fig3b.colorbar(surf, ax=ax3b)
Env.plot_turbines(ax3b)
#
#with open(os.path.join(direc, "solutions/wfl_sol_20_20_1.npy"), "rb") as bla:
#    layout = np.load(bla)
#    layout_ind  = np.load(bla)
#for node_index in layout_ind:  
#    ax3b.scatter(Env.grid[node_index][0], Env.grid[node_index][1], s=100, alpha=0.2, facecolor="blue", edgecolor="black")
#
plt.xlabel("x")
plt.ylabel("y")
#plt.savefig(os.path.join(fig_dir, "geo_sol.png"))
plt.close(fig3b)


"""
Plot some interference
"""
LO = layout_optimization(axx,axy,WindInd=1)
fig4, ax4 = plt.subplots()
heat1 = plt.imshow(LO.infer_matrix[int(axx*axy/2+axx/2)], cmap='jet', interpolation='bicubic')
plt.colorbar(heat1)
plt.xlabel("x")
plt.ylabel("y")
plt.gca().invert_yaxis()
#plt.savefig(os.path.join(fig_dir, "ex_inter12.png"))
plt.close(fig4)

"""
Plot discretized grid
"""
fig5, ax5 = plt.subplots()
Env.plot_grid(ax5, numbers=False)
plt.xlabel("x")
plt.ylabel("y")
#plt.savefig(os.path.join(fig_dir, "grid.png"))
plt.close(fig5)

"""
Plot distance gradient from shore
"""
fig6, ax6 = plt.subplots()
heat = plt.imshow(Env.dist_matrix, cmap='hot', interpolation='bicubic')
plt.colorbar(heat)
plt.xlabel("x")
plt.ylabel("y")
#plt.savefig(os.path.join(fig_dir, "dist_gradient.png"))
plt.close(fig6)

"""
Plot initial solution
- falsche x achse
"""
fig7, ax7 = plt.subplots()
Env.plot_turbines(ax7, which="initial")
plt.xlabel("x")
plt.ylabel("y")
ax7.set_xticks([0,2.5,5,7.5,10,12.5,15,17.5])
#plt.savefig(os.path.join(fig_dir, "initial.png"))
plt.close(fig7)


plt.show()


# ###################
# CALCULATIONS
# ###################

Env.eval_obj(1,3,6,12,1)
