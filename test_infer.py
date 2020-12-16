import numpy as np
import matplotlib.pyplot as plt
from wfl_opt_data import *

nodeoi = 23
plt.scatter(grid[nodeoi][0], grid[nodeoi][1])
plt.imshow(infer_matrix[nodeoi], cmap='hot', interpolation='bicubic')
plt.show()