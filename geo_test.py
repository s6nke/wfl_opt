import numpy as np
import matplotlib.pyplot as plt

axx = 20
axy = 20
dist_matrix = np.zeros(shape=(axx, axy))

for i in range(axx):
    dist_matrix[:, axx-1-i] = i

x = np.arange(0,axx)
y = np.arange(0,axy)
xm, ym = np.meshgrid(x,y)
geo_matrix = np.sin(xm+ym)-1

#zz = z + geo_matrix/10


fig = plt.figure()
ax = fig.gca(projection='3d')
surf = ax.plot_surface(xm,ym,geo_matrix, cmap="hot")
fig.colorbar(surf, ax=ax)
plt.show()