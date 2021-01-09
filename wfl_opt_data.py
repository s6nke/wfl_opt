#-------------------------------------
# Data file for optimizing a offshore wind farm layout
#-------------------------------------
#
# TODO:
# - adapt power function to represent nonlinear behaviour depending on wind speed
# - verify optimality
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
import os
import numpy as np
import matplotlib.pyplot as plt
import docplex.mp.model as cpx  

# ################################
# BEGIN of optimization base class
class wf_environment:
    def __init__(self,axx,axy):
        """
        Initialize Environment object with
        axx: int, dimension in x direction
        axy: int, dimension in y direction
        """
        # init variables
        self.axx = axx
        self.axy = axy
        self.n   = axx*axy

        # init calculations
        self.setV = range(0, self.n)
        self.init_grid()

        # relative path init
        self.dir = os.path.dirname(os.path.abspath(__file__))
        self.dir_sol = os.path.join(self.dir, "solutions")
    # END __INIT___

    def init_grid(self):
        """
        Construct grid nodes
        """
        print("Generate grid notes ....")
        grid = []
        for ii in range(0,self.axx):
            for jj in range(0,self.axy):
                grid.append([ii,jj])
        self.grid = np.asarray(grid)
        print("Done!")
    # INIT_GRID

    def dist(self, node1, node2):
        """
        Kartesian distance of two nodes
        """
        return np.linalg.norm(node1-node2)
    # END DIST
   
    def load_layout_sol(self):
        """
        Initialize optimal layout from file
        """
        print("Load layout file ...")
        filename = "wfl_sol_" + str(self.axx) + "_" + str(self.axy) + ".npy"
        try:
            with open(os.path.join(self.dir_sol,filename), "rb") as layout_file:
                # the hole grid with 0 and 1
                self.layout_sol = np.load(layout_file)
                # only the indices of the turbine nodes
                self.layout_sol_indices = np.load(layout_file)
                # interference matrix of the solution
                self.inf_sol = np.load(layout_file)
            print("Done!")
        except:
            print("No mathing file found!")
    # END LOAD_LAYOUT_SOL
    
    def load_cable_sol(self):
        """
        Initialize optimal cable routing from file
        """
        print("Load cable routing file!")
        filename = "cable_sol_20_20.csv"
        try:
            with open(os.path.join(self.dir_sol,filename), "rb") as cr_file:
                self.cr_arcs0 = np.load(cr_file)
                self.cr_arcs1 = np.load(cr_file)
                self.cr_arcs2 = np.load(cr_file)
            print("Done!")
        except:
            print("No matching file found!")
    # END LOAD_CABLE_SOL

    def plot_turbines(self, ax, col="black"):
        """
        Plot helper function for turbine nodes
        """
        for node_index in self.layout_sol_indices:
            # also adds a Dmin-diameter circle around this point to show forbidden zones
            ax.scatter(self.grid[node_index][0], self.grid[node_index][1], s=20, c=col, marker="o")
            #ax.add_artist(plt.Circle((grid[i][0], grid[i][1]), Dmin, alpha=0.1))
        return ax
    # END PLOT_TURBINES

    def plot_arcs(self):
        """
        plot helper function for arcs of the solution
        """
        for vec in self.cr_arcs0:
            x0,y0 = self.grid[vec[0]]
            x1,y1 = self.grid[vec[1]]
            dx = x1-x0
            dy = y1-y0
            plt.arrow(x0,y0,dx,dy, color="yellow")
        for vec in self.cr_arcs1:
            x0,y0 = self.grid[vec[0]]
            x1,y1 = self.grid[vec[1]]
            dx = x1-x0
            dy = y1-y0
            plt.arrow(x0,y0,dx,dy, color="orange")
        for vec in self.cr_arcs2:
            x0,y0 = self.grid[vec[0]]
            x1,y1 = self.grid[vec[1]]
            dx = x1-x0
            dy = y1-y0
            plt.arrow(x0,y0,dx,dy, color="red")
    # END PLOT_ARCS
    


# END of optimization base class
# ##############################


# ##############################
# BEGIN OF WFL ENVIRONMENT CLASS

class layout_optimization(wf_environment):
    def __init__(self, axx, axy, Dmin=5, k=0.05, V=15, D=1, v_wind=15, w_dirs=[1,1]):
        """
        Initialize layout object
        """
        super().__init__(axx,axy)
        self.Dmin = Dmin
        self.k = k
        self.V = V
        self.D = D
        self.v_wind = v_wind
        self.w_dirs  = w_dirs

        # init calculations
        self.init_interference_matrix()
        self.init_dist_matrix()
        self.init_geo_matrix()
        self.init_setE()
    # END __INIT__

    def init_setE(self):
        """
        Generate the set E of all nodes that are in the minimum
        distance to the node i
        """
        print("Generate set E_i ....")
        set_E = []
        not_set_E = []
        for l in range(0,self.n):
            # create empty sub-list
            set_E.append([])
            not_set_E.append([])
            for k in range(0,self.n):
                # check if the node is in Dmin to our node[l]
                if self.dist(self.grid[l], self.grid[k]) <= self.Dmin:
                    # if too close add to set_E sublist
                    set_E[l].append(k)
                else:
                    not_set_E[l].append(k)
            # remove the node itself from the set to prevent false addition if set true
            try:
                set_E[l].remove(l)
            except:
                not_set_E[l].remove(l)
        self.setE = np.array(set_E, dtype=object)
        self.not_setE = np.array(not_set_E, dtype=object)
        print("Done!")
    # END INIT_SETE

    def init_dist_matrix(self):
        """
        calculate linear distance matrix from shore (east)
        """
        print("Calculate distance matrix ...")
        dist_matrix = np.zeros(shape=(self.axx, self.axy))
        for i in range(self.axx):
            dist_matrix[:, self.axx-1-i] = i
        self.dist_matrix = dist_matrix
        print("Done!")
    # END INIT_DIST_MATRIX

    def init_geo_matrix(self):
        """
        calculate some arbirtrary ground depth
        """
        print("Calculate geo matrix ...")
        x = np.arange(0,self.axx)
        y = np.arange(0,self.axy)
        xm, ym = np.meshgrid(x,y)
        self.geo_matrix = np.sin(xm+ym)-1
        print("Done!")
    # END INIT_GEO_MATRIX

    def init_interference_matrix(self):
        print("Open interference file ...")
        try:
            filename = "interference_matrix_" + str(self.axx) + "_" + str(self.axy) + "_" + str(self.n) + ".npy"
            with open(os.path.join(self.dir_sol,filename), "rb") as infer_file:
                self.infer_matrix = np.load(infer_file)
            print("Done!")
        except:
            print("no file found, calulated now automatically...")
            self.calc_interference()
            print("Done!")
    # END INIT_WARMSTART

    def warmstart(self, OP, x_vars, Nmin, Nmax):
        print("Initialise warmstart ...")
        self.find_initial_sol(Nmin, Nmax)
        warmstart = OP.new_solution()
        for el in self.initial_sol:
            warmstart.add_var_value(x_vars[el], el)
        OP.add_mip_start(warmstart)
        print("Done!")
    # END WARMSTART

    def plot_grid(self, axs, numbers=True):
        # grid: array of nodes; array of array-like
        for i in range(0,self.n):
            # scatter the node points
            axs.scatter(self.grid[i][0], self.grid[i][1], s=10, c="black")
            # add number to the node
            if numbers:
                axs.text(self.grid[i][0], self.grid[i][1], "{0}".format(i))
        return axs
    # END PLOT_GRID
    
    def sol_interference(self):
        """
        Helper function for interference of solution -> plotting
        """
        ifm_sol = np.zeros(shape=(self.axx,self.axy))
        for i in range(len(self.sol)):
            if self.sol[i] == 1.0:
                ifm_sol = ifm_sol + self.infer_matrix[i]
        self.sol_inf = ifm_sol
    # END SOL_INTERFERENCE

    def P(self, v_wind):
        """
        Linear power function of turbine
        """
        if v_wind <= 3:
            return 0
        elif 3 < v_wind <= 16:
            return 2.3/13*(v_wind-3)
        elif v_wind > 16:
            return 2.3
    # END P

    def jensen_mod(self, upwind_stream, rotor_diam, distance):
        """
        Jensens wake model function
        """
        k = 0.05 # offshore wake decay coefficient
        return upwind_stream*(1 - np.sqrt(1 - self.thrust_coeff(upwind_stream)))*(rotor_diam/(rotor_diam+2*k*distance))**2
    # END JENSENS_MOD

    def thrust_coeff(self, wind_speed):
        """
        Approximation of the thrust coefficicent dependent of the wind speed
        """
        mod_coeff = 6 # chosen arbitrary to get quite good values
        if wind_speed <= 3:
            return 0
        else:
            return mod_coeff/(wind_speed)
    # END THURST_COEFF
    
    def calc_interference(self):
        """
        Calculate interence matrix based on jensens wake model
        entries are the loss of power at each node relative to node i
        """
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
    # CALC_INTERFERENCE
    
    def save_interference(self, ifmxx):
        """
        Save the interference matrix to a file
        """
        filename = "interference_matrix_" + str(self.axx) + "_" + str(self.axy) + "_" + str(self.n) + ".npy" 
        print("Save interference matrix to " + filename)
        with open(os.path.join(self.dir_sol,filename), "wb") as infer_file:
            np.save(infer_file, ifmxx)
        print("Done!")
    # END SAVE_INTERFERENCE
    
    def find_initial_sol(self, Nmin, Nmax):
        """
        Search for an unconstrained solution to warmstart the solver
        """
        OP_initial = cpx.Model(name="Wind Farm Layout", log_output=True)
        x_vars = {(i): OP_initial.binary_var(name="x_{0}".format(i)) for i in self.setV}
        OP_initial.add_constraint(OP_initial.sum(x_vars[i] for i in self.setV) <= Nmax)
        OP_initial.add_constraint(OP_initial.sum(x_vars[i] for i in self.setV) >= Nmin)
        for cnt in range(0,self.n):
            for el in self.setE[cnt]:
                OP_initial.add_constraint(x_vars[cnt] + x_vars[el] <= 1)
        obj = OP_initial.sum(self.P(15)*x_vars[i] for i in self.setV)
        OP_initial.maximize(obj)
        OP_initial.solve()
        self.initial_sol = [OP_initial.solution.get_value(x_vars[element]) for element in self.setV]
    # END FIND_INITIAL_SOL

    def save_sol(self):
        """
        Save solution of optimizaton method to file
        """
        filename = "../solutions/wfl_sol_" + str(self.axx) + "_" + str(self.axy) + ".npy"
        with open(os.path.join(self.dir_sol, filename), "wb") as sol_file:
            np.save(sol_file, self.sol)
            np.save(sol_file, self.sol_indices)
            np.save(sol_file, self.sol_inf)
    # END SAVE_SOL



# END OF WFL_ENVIRONMENT CLASS
# ############################

# ###################################
# START OF CABLE ROUTING OPTIMIZATION

class cable_routing_optimization(wf_environment):
    def __init__(self,axx,axy,V0):
        """ 
        Init of cr object
        axx: int, size of x axis
        axy: int, size of y axis
        V0 : list of substation nodes
        """
        # call init from environment
        super().__init__(axx,axy)
        
        # load sol file
        self.load_layout_sol()
        
        # setVT is the set of node indices from sol turbines
        self.setVT = self.layout_sol_indices

        # init substation indices 
        self.setV0 = V0

        # VT0 is turbines and substations; has to be unique
        self.setVT0 = np.unique(np.append(self.setVT,self.setV0))
    # END __INIT__

    def calc_setA(self):
        """
        define a set of all possible arcs between turbine nodes
        """
        setA = []
        for i in self.setVT0:
            for j in self.setVT0:
                if i!=j:
                    setA.append([i,j])
        self.setA = setA
    # END CALC_SETA

    def cost(self,arc):
        """
        deine a cost function based on the distance of the nodes
        -> i.e. length of the arc
        """
        node0 = self.grid[arc[0]]
        node1 = self.grid[arc[1]]
        return self.dist(node0, node1)
    # END COST
    
    def plot_substations(self):
        """
        plot helper functon for substations
        """
        for ss in self.setV0:
            plt.scatter(self.grid[ss][0],self.grid[ss][1], c='blue')
    # END PLOT_SUBSTATIONS

    def save_sol(self):
        """
        save the solution to a numpy file
        """
        print("Save solution to file ...")
        filename = "cable_sol_" + str(self.axx) + "_" + str(self.axy) + ".csv"
        #pd.DataFrame(self.sol_arcs).to_csv(os.path.join(self.dir_sol, filename))
        with open(os.path.join(self.dir_sol,filename), "wb") as cr_sol:
            np.save(cr_sol, self.sol_arcs0)
            np.save(cr_sol, self.sol_arcs1)
            np.save(cr_sol, self.sol_arcs2)
        print("Done!")
    # END SAVE_SOL

# END OF CABLE ROUTING OPTIMIZATION
# ###################################

