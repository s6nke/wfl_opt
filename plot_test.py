import matplotlib.pyplot as plt
import numpy as np


def plott(ax):
    for i in range(len(y)):
        ax.scatter(t[i],y[i])
    return ax

t = np.linspace(0,10,100)
y = np.sin(t)

fig, ax = plt.subplots(2,1)
plott(ax[0])
ax[1].scatter(t,y)
plt.show()