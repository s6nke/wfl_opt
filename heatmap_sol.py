import numpy as np
import matplotlib.pyplot as plt


with open("interference_matrix.npy","rb") as interference_file:
    with open("opt_sol.npy", "rb") as sol_file:
        # load sets
        ifm = np.load(interference_file)
        sol = np.load(sol_file)

        ifm_sol = np.zeros(shape=(10,10))
        for i in range(len(sol)):
            if sol[i] == 1.0:
                ifm_sol = ifm[i]

        plt.imshow(ifm_sol, cmap='hot', interpolation='gaussian')
        #plt.show()

plt.imshow()