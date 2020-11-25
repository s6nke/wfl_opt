# ------------
# DEPENDENCIES
# ------------
import matplotlib.pyplot as plt
from matplotlib.path import Path
from matplotlib.patches import PathPatch
import numpy as np

# ---------
# VARIABLES
# ---------
V = 15      # m/s
D = 50      # rotor diameter

# for plotting
xAxis = 1000
yAxis = 1000
axx = np.linspace(0, xAxis, 101)
axy = np.linspace(0, yAxis, 101)
X,Y = np.meshgrid(axx, axy)
wind_speeds = np.linspace(0,20,101)
power = []
ct = []

# ---------
# FUNCTIONS
# ---------

# Jensens Model
def jensen_mod(upwind_stream, rotor_diam, distance):
    k = 0.05 # wake decay coefficient
    return upwind_stream - upwind_stream*(1 - np.sqrt(1 - thrust_coeff(upwind_stream)))*(rotor_diam/(rotor_diam+2*k*distance))**2

# function of thrust coefficient
def thrust_coeff(wind_speed):
    mod_coeff = 6 # chosen arbitrary to get quite good values
    if wind_speed <= 3:
        return 0
    else:
        return mod_coeff/(wind_speed)

# linear approx. of the power function
def P(v_wind):
    if v_wind <= 3:
        return 0
    elif 3 < v_wind <= 16:
        return 2.3/13*(v_wind-3)
    elif v_wind > 16:
        return 2.3

# --------
# PLOTTING
# --------
# calculate discrete values
for v in wind_speeds:
    power.append(P(v))
    ct.append(thrust_coeff(v))

def top(x_dist, yAxis, rotor_diam):
    return 0.05*x_dist + yAxis/2 + rotor_diam/2
def bottom(x_dist, yAxis, rotor_diam):
    return -0.05*x_dist + yAxis/2 - rotor_diam/2

ttop = top(axx, yAxis, D)
bbottom = bottom(axx, yAxis, D)


jensens = jensen_mod(15, D, X)

path = Path([[0, 500+D/2], [1000, max(ttop)], [1000, min(bbottom)], [0, 500-D/2]])
patch = PathPatch(path, facecolor='none', edgecolor='none')

# create figure
fig, ax = plt.subplots(2)
ax[0].plot(wind_speeds, power)
ax[0].plot(wind_speeds, ct)
ax[1].add_patch(patch)
im = plt.imshow(jensens, interpolation='bilinear', cmap=plt.cm.gray,
                origin='lower', extent=[0, 1000, 250, 750],
                clip_path=patch, clip_on=True)
im.set_clip_path(patch)
fig.colorbar(im, ax=ax[1])
plt.show()
