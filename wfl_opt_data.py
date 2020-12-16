#-------------------------------------
# Data file for optimizing a offshore wind farm layout
#-------------------------------------
#
# TODO:
# - adapt power function to represent nonlinear behaviour depending on wind speed
# - verify optimality
# - include cable-routing problem
# - make some comments
# 
# ----------------
# LIST OF METHODS:
# ----------------
# init_grid()
# init_setE()
# init_interference_matrix()
#   - calc_interference()
#   - jensens_mod()
#   - thrust_coeff()
#   - save_interference()
# warmstart()
#   - find_initial_sol()
# dist()
# plot_grid()
# plot_interference()
# plot_turbines()
# P()
# save_sol()
# load_sol()
#
# ------------
# DEPENDENCIES
# ------------
import random
import numpy as np
import matplotlib.pyplot as plt
import docplex.mp.model as cpx  

# ##############################
# BEGIN OF WFL ENVIRONMENT CLASS

class wfl_environmet:
    def __init__(self, axx, axy, Dmin=5, k=0.05, V=15, D=1, v_wind=15, w_dirs=[1,1]):
        # init variables
        self.axx = axx
        self.axy = axy
        self.n   = axx*axy
        self.Dmin = Dmin
        self.k = k
        self.V = V
        self.D = D
        self.v_wind = v_wind
        self.w_dirs  = w_dirs

        # init calculations
        self.init_grid()
        self.init_setE()
        self.init_interference_matrix()
        self.setV = range(0, self.n)
    
    def init_grid(self):
        # construct grid nodes for geometry of the ares
        # assigns self the grid numpy array
        print("Generate grid notes ....")
        grid = []
        for ii in range(0,self.axx):
            for jj in range(0,self.axy):
                grid.append([ii,jj])
        self.grid = np.array(grid)
        print("Done!")

    def init_setE(self):
        # generate the set E_i
        # elements of each subset are indexes of the set_V that are in
        # a distance of Dmin around turbine i
        print("Generate set E_i ....")
        set_E = []
        for l in range(0,self.n):
            # create empty sub-list
            set_E.append([])
            for k in range(0,self.n):
                # check if the node is in Dmin to our node[l]
                if self.dist(self.grid[l], self.grid[k]) <= self.Dmin:
                    # if too close add to set_E sublist
                    set_E[l].append(k)
            # remove the node itself from the set to prevent false addition if set true
            set_E[l].remove(l)
        self.setE = np.array(set_E)
        print("Done!")

    def init_interference_matrix(self):
        print("Open interference file ...")
        try:
            with open("interference_matrix_" + str(self.axx) + "_" + str(self.axy) + "_" + str(self.n) + ".npy", "rb") as infer_file:
                self.infer_matrix = np.load(infer_file)
            print("Done!")
        except:
            print("no file found, calulated now automatically...")
            self.calc_interference()
            print("Done!")

    def warmstart(self, OP, x_vars, Nmin, Nmax):
        self.find_initial_sol(Nmin, Nmax)
        warmstart = OP.new_solution()
        for el in self.initial_sol:
            warmstart.add_var_value(x_vars[el], el)
        OP.add_mip_start(warmstart)

    # calculate kartesian difference of two grid nodes
    def dist(self, node1, node2):
        return np.linalg.norm(node1-node2)

    # plot grid nodes
    def plot_grid(self, numbers=True):
        # grid: array of nodes; array of array-like
        for i in range(0,self.n):
            # scatter the node points
            plt.scatter(self.grid[i][0], self.grid[i][1], s=10, c="black")
            # add number to t   he node
            if numbers:
                plt.text(self.grid[i][0], self.grid[i][1], "{0}".format(i))
    
    def plot_turbines(self, sr, col="white"):
        for i in range(0, self.n):
        # plot a red star on top of the grid node if the decsision variable is set 1
        # also adds a Dmin-diameter circle around this point to show forbidden zones
            if sr[i] == 1.0:
                plt.scatter(self.grid[i][0], self.grid[i][1], s=100, c=col, marker="*")
                #ax.add_artist(plt.Circle((grid[i][0], grid[i][1]), Dmin, alpha=0.1))

    # plot heatmap of interference
    def plot_interference(self, sol):
        ifm_sol = np.zeros(shape=(self.axx,self.axy))
        for i in range(len(sol)):
            if sol[i] == 1.0:
                ifm_sol = ifm_sol + self.infer_matrix[i]

        heat = plt.imshow(ifm_sol, cmap='jet', interpolation='bilinear')
        plt.colorbar(heat)

    def P(self, v_wind):
        if v_wind <= 3:
            return 0
        elif 3 < v_wind <= 16:
            return 2.3/13*(v_wind-3)
        elif v_wind > 16:
            return 2.3

    def jensen_mod(self, upwind_stream, rotor_diam, distance):
        k = 0.05 # offshore wake decay coefficient
        return upwind_stream*(1 - np.sqrt(1 - self.thrust_coeff(upwind_stream)))*(rotor_diam/(rotor_diam+2*k*distance))**2

    def thrust_coeff(self, wind_speed):
        mod_coeff = 6 # chosen arbitrary to get quite good values
        if wind_speed <= 3:
            return 0
        else:
            return mod_coeff/(wind_speed)
    
    def calc_interference(self):
        # calculate interence matrix based on jensens wake model
        # entries are the loss of power at each node relative to node i
        print("Creating interference matrix ...")
        interferenz_matrix = np.zeros(shape=(self.n, self.axx, self.axy))
        for wind_direction in self.w_dirs:
            for i in range(0, self.n):
                infer  = []
                for j in range(0, self.n):
                    alpha   = np.matmul(np.transpose(self.grid[j] - self.grid[i]), wind_direction) / np.linalg.norm(wind_direction)**2
                    if alpha < 0:
                        infer.append(0)
                    else:
                        P_prime = self.grid[i] + alpha*wind_direction
                        delta_1 = self.dist(self.grid[i], P_prime)
                        delta_2 = self.dist(self.grid[j], P_prime)
                        if delta_2 <= (self.D+2*self.k*delta_1)/2:
                            infer.append( self.P(self.v_wind) - self.P(self.v_wind - self.jensen_mod(self.V, self.D, delta_1)) )
                        else:
                            infer.append(0)
                interferenz_matrix[i] += (1/self.w_dirs.shape[0])*np.array(infer).reshape((self.axx, self.axy)).T
            
        print("Done!")
        self.save_interference(interferenz_matrix)
        self.init_interference_matrix()
        return interferenz_matrix
    
    def save_interference(self, ifmxx):
        filename = "interference_matrix_" + str(self.axx) + "_" + str(self.axy) + "_" + str(self.n) + ".npy" 
        print("Save interference matrix to " + filename)
        with open(filename, "wb") as infer_file:
            np.save(infer_file, ifmxx)
        print("Done!")
    
    def find_initial_sol(self, Nmin, Nmax):
        OP_initial = cpx.Model(name="Wind Farm Layout", log_output=True)
        x_vars = {(i): OP_initial.binary_var(name="x_{0}".format(i)) for i in self.setV}
        OP_initial.add_constraint(OP_initial.sum(x_vars[i] for i in self.setV) <= Nmax)
        OP_initial.add_constraint(OP_initial.sum(x_vars[i] for i in self.setV) >= Nmin)
        for cnt in range(0,self.n):
            for el in self.setE[cnt]:
                OP_initial.add_constraint(x_vars[cnt] + x_vars[el] <= 1)
        obj =   OP_initial.sum(self.P(15)*x_vars[i] for i in self.setV)
        OP_initial.maximize(obj)
        OP_initial.solve()
        self.initial_sol = [OP_initial.solution.get_value(x_vars[element]) for element in self.setV]

    def save_sol(self):
        filename = "opt_sol_" + str(self.axx) + "_" + str(self.axy) + ".npy"
        with open(filename, "wb") as sol_file:
            np.save(sol_file, self.sol)

    def load_sol(self):
        filename = "opt_sol_" + str(self.axx) + "_" + str(self.axy) + ".npy"
        try:
            with open(filename, "rb") as ifm_file:
                self.sol = np.load(ifm_file, allow_pickle=True)
        except:
            print("No mathing file found!")


# END OF WFL_ENVIRONMENT CLASS
# ############################
