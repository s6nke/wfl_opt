#-------------------------------------
# MODULE file for optimizing a offshore wind farm layout
#-------------------------------------
#
# TODO:
# - verify optimality
# - make some comments
# - set A new calc, when j > i
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
        self.dir_sol  = os.path.join(self.dir, "solutions")
        self.dir_wind = os.path.join(self.dir, "wind_files")
    # END __INIT___

    def init_grid(self):
        """
        Construct grid nodes
        """
        print("Generate grid notes ....")
        grid = []
        for i in range(0,self.axx):
            for j in range(0,self.axy):
                grid.append([i,j])
        self.grid = np.asarray(grid)
        print("Done!")
    # INIT_GRID

    def dist(self, node1, node2):
        """
        Kartesian distance of two nodes
        """
        return np.linalg.norm(node1-node2)
    # END DIST

    def eval_obj(self, lay1, lay2, lay3, lay4, windcase):
        """
        Evaluation function to compare the objective value of different layouts to different
        wind scenarios
        """
        layout_1, P1, obj1 = self.load_layout_sol(lay1, ret=True)
        layout_2, P2, obj2 = self.load_layout_sol(lay2, ret=True)
        layout_3, P3, obj3 = self.load_layout_sol(lay3, ret=True)
        layout_4, P4, obj4 = self.load_layout_sol(lay4, ret=True)
        inf_data = self.load_infer_matrix(windcase, ret=True)

        inf1 = np.zeros((self.axx,self.axy))
        for node_i in layout_1:
            xi, yi = self.grid[node_i]
            inf1 += inf_data[node_i]
            inf1[yi,xi] = inf1[yi,xi] - inf_data[node_i, yi,xi] 

        inf2 = np.zeros((self.axx,self.axy))
        for node_i in layout_2:
            xi, yi = self.grid[node_i]
            inf2 += inf_data[node_i]
            inf2[yi,xi] = inf2[yi,xi] - inf_data[node_i, yi,xi]

        inf3 = np.zeros((self.axx,self.axy))
        for node_i in layout_3:
            xi, yi = self.grid[node_i]
            inf3 += inf_data[node_i]
            inf3[yi,xi] = inf3[yi,xi] - inf_data[node_i, yi,xi]

        inf4 = np.zeros((self.axx,self.axy))
        for node_i in layout_4:
            xi, yi = self.grid[node_i]
            inf4 += inf_data[node_i]
            inf4[yi,xi] = inf4[yi,xi] - inf_data[node_i, yi,xi]   

        w1 = []
        for node_i in layout_1:
            xj,yj = self.grid[node_i]
            w1.append(inf1[yj,xj])
        w2 = []
        for node_i in layout_2:
            xj,yj = self.grid[node_i]
            w2.append(inf2[yj,xj])
        w3 = []
        for node_i in layout_3:
            xj,yj = self.grid[node_i]
            w3.append(inf3[yj,xj])
        w4 = []
        for node_i in layout_4:
            xj,yj = self.grid[node_i]
            w4.append(inf4[yj,xj])

        print(np.round(np.sum(P1 - w1),4), np.round(np.sum(P1 - w1)/obj1, 4), np.round(np.sum(P1 - w1)/obj4, 4))
        print(np.round(np.sum(P2 - w2),4), np.round(np.sum(P2 - w2)/obj2, 4), np.round(np.sum(P2 - w2)/obj4, 4))
        print(np.round(np.sum(P3 - w3),4), np.round(np.sum(P3 - w3)/obj3, 4), np.round(np.sum(P3 - w3)/obj4, 4))
        print(np.round(np.sum(P4 - w4),4), np.round(np.sum(P4 - w4)/obj4, 4), np.round(np.sum(P4 - w4)/obj4, 4))    

        temp = 0
        for node_i in layout_1:
            xi,yi = self.grid[node_i]
            temp += self.geo_matrix[yi,xi]*20
        print(temp, np.round(np.sum(P1 - w1),4), np.round(np.sum(P1 - w1),4)+temp)
    # END EVAL_OBJ
   
    def load_layout_sol(self, WindInd, ret=False):
        """
        Initialize optimal layout from file
        """
        print("Load layout file ...")
        filename = "wfl_sol_dual_" + str(self.axx) + "_" + str(self.axy) + "_" + str(WindInd) + ".npy"
        try:
            with open(os.path.join(self.dir_sol,filename), "rb") as layout_file:
                # the hole grid with 0 and 1
                self.layout_sol = np.load(layout_file)
                # only the indices of the turbine nodes
                self.layout_sol_indices = np.load(layout_file)
                # interference matrix of the solution
                self.inf_sol = np.load(layout_file)
                # objective value
                self.layout_sol_obj = np.load(layout_file)
                # 
                self.layout_Pi = np.load(layout_file)
            print("Done!")
        except:
            print("No mathing file found!")
        if ret:
            return [self.layout_sol_indices, self.layout_Pi, self.layout_sol_obj]
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
                self.setV0 = np.load(cr_file)
            print("Done!")
        except:
            print("No matching file found!")
    # END LOAD_CABLE_SOL

    def load_infer_matrix(self,WindInd,ret=False):
        print("Open interference file ...")
        try:
            filename = "interference_matrix_" + str(self.axx) + "_" + str(self.axy) + "_" + str(self.n) + "_" + str(WindInd) + ".npy"
            with open(os.path.join(self.dir_sol,filename), "rb") as infer_file:
                self.infer_matrix = np.load(infer_file)
            print("Done!")
        except:
            print("Failed!")
        if ret:
            return self.infer_matrix
    # END INIT_INTERFERENCE_MATRIX

    def load_dist_geo(self):
        print("Load distance and geo file")
        filename = "dist_geo_20_20.npy"
        try:
            with open(os.path.join(self.dir_sol,filename), "rb") as dg_file:
                self.dist_matrix = np.load(dg_file)
                self.geo_matrix  = np.load(dg_file)
        except:
            print("Unable to read file!")
    # END LOAD_DIST_GEO

    def load_initial_sol(self):
        print("Loading initial solution")
        filename = "wfl_initial_sol_20_20.npy"
        try:
            with open(os.path.join(self.dir_sol,filename), "rb") as ini_file:
                self.initial_sol = np.load(ini_file)
        except:
            print("something failed")

    def plot_turbines(self, ax, which="sol", col="black"):
        """
        Plot helper function for turbine nodes
        """
        if which == "sol":
            for node_index in self.layout_sol_indices:  
                ax.scatter(self.grid[node_index][0], self.grid[node_index][1], s=100, facecolor="white", edgecolor=col)
        if which == "initial":
            for node_index in self.initial_sol:
                # also adds a Dmin-diameter circle around this point to show forbidden zones
                ax.scatter(self.grid[node_index][0], self.grid[node_index][1], s=20, facecolor="white", edgecolor=col)
                ax.add_artist(plt.Circle((self.grid[node_index][0], self.grid[node_index][1]), 2, alpha=0.1))
        return ax
    # END PLOT_TURBINES

    def plot_grid(self, axs, numbers=True):
        # grid: array of nodes; array of array-like
        for i in range(0,self.n):
            # scatter the node points
            axs.scatter(self.grid[i][0], self.grid[i][1], marker="o", s=10, facecolor="white", edgecolor="black")
            # add number to the node
            if numbers:
                axs.text(self.grid[i][0], self.grid[i][1], "{0}".format(i))
        return axs
    # END PLOT_GRID

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

    def plot_substations(self,ax):
        """
        plot helper functon for substations
        """
        for ss in self.setV0:
            ax = plt.scatter(self.grid[ss][0],self.grid[ss][1], c='blue')
        return ax
    # END PLOT_SUBSTATIONS
#
#
#
# END of optimization base class
# ##############################


# ##############################
# BEGIN OF WFL ENVIRONMENT CLASS

class layout_optimization(wf_environment):
    def __init__(self, axx, axy, WindInd, Dmin=5, k=0.05, D=1):
        """
        Initialize layout object
        """
        super().__init__(axx,axy)
        self.Dmin = Dmin
        self.k = k
        self.D = D
        self.init_wind_data(WindInd)
        self.WindInd = WindInd

        # init calculations
        self.Pi()
        self.init_interference_matrix()
        self.init_dist_matrix()
        self.init_geo_matrix()
        self.save_dist_geo()
        self.init_setE()
    # END __INIT__

    def Pi(self):
        vi = 0
        for i in range(self.v_wind.shape[0]):
            vi += self.v_wind[i]*(1/self.v_wind.shape[0])
        self.Pi = self.P(vi)
    # END PI

    def init_wind_data(self,ind):
        if ind == 1:
            self.w_dirs = np.asarray([[-0.64278761, -0.766044443]], dtype=object)
            self.v_wind = np.asarray([8])
        if ind == 3:
            self.w_dirs = np.asarray([[0.5, 0.8660254], [-1,0], [-0.5, -0.866025404]], dtype=object)
            self.v_wind = np.asarray([8, 8, 8])
        if ind == 6:
            self.w_dirs = np.asarray([[0.866025404, 0.5], [0,1], [-0.866025404, 0.5], [-0.866025404, -0.5], [0,-1], [0.866025404, -0.5]], dtype=object)
            self.v_wind = np.asarray([8]*6)
        if ind == 12:
            self.w_dirs = np.asarray([[0.965925826, 0.258819045], [0.707106781,	0.707106781], [-0.258819045,	0.965925826], [-0.707106781,	0.707106781],
                            [-0.965925826,	0.258819045],  [-0.965925826,	-0.258819045], [-0.707106781,	-0.707106781], [-0.258819045,	-0.965925826], 
                            [0.258819045,	-0.965925826], [0.707106781,	-0.707106781], [0.965925826,	-0.258819045]], dtype=object) 
            self.v_wind = np.asarray([8]*12)
    # END INIT_WIND_DATA

    def init_setE(self):
        """
        Generate the set E of all nodes that are in the minimum
        distance to the node i
        """
        print("Generate set E_i ....")
        set_E = []
        not_set_E = []
        for i in self.setV:
            # create empty sub-list
            set_E.append([])
            not_set_E.append([])
            for j in self.setV:
                if j > i:
                    # check if the node is in Dmin to our node[i]
                    if self.dist(self.grid[i], self.grid[j]) < self.Dmin:
                        # if too close add to set_E sublist
                        set_E[i].append(j)
                    else:
                        not_set_E[i].append(j)
                else:
                    not_set_E[i].append(j)
        self.setE = np.array(set_E, dtype=object)
        self.not_setE = np.array(not_set_E, dtype=object)
        print("Done!")
    # END INIT_SETE

    def init_dist_matrix(self):
        """ 
        calculate linear distance matrix from shore (east)
        """
        print("Calculate distance matrix ...")
        ppm_offshore = 5/self.axx
        dist_matrix = np.zeros(shape=(self.axx, self.axy))
        for i in range(self.axx):
            dist_matrix[:, self.axx-1-i] = i*ppm_offshore
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
        self.geo_matrix = (np.sin((xm)/5)-1)*10
        print("Done!")
    # END INIT_GEO_MATRIX

    def init_interference_matrix(self):
        print("Open interference file ...")
        try:
            filename = "interference_matrix_" + str(self.axx) + "_" + str(self.axy) + "_" + str(self.n) + "_" + str(self.WindInd) + ".npy"
            with open(os.path.join(self.dir_sol,filename), "rb") as infer_file:
                self.infer_matrix = np.load(infer_file)
            print("Done!")
        except:
            print("no file found, calulated now automatically...")
            self.calc_interference()
            print("Done!")
    # END INIT_INTERFERENCE_MATRIX

    def save_dist_geo(self):
        print("Save distance and geo matrices!")
        try:
            filename = "dist_geo_" + str(self.axx) + "_" + str(self.axy) + ".npy"
            with open(os.path.join(self.dir_sol,filename), "wb") as dg_file:
                np.save(dg_file, self.dist_matrix)
                np.save(dg_file, self.geo_matrix)
        except:
            print("unable to solve file")
    # END SAVE_DIST_GEO

    def warmstart(self, OP, x_vars, Nmin, Nmax):
        print("Initialise warmstart ...")
        self.find_initial_sol(Nmin, Nmax)
        warmstart = OP.new_solution()
        for el in self.initial_sol:
            warmstart.add_var_value(x_vars[el], el)
        OP.add_mip_start(warmstart)
        print("Done!")
    # END WARMSTART
    
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
        if v_wind <= 1:
            return 0
        elif 1 < v_wind <= 16:
            return 2300/15*(v_wind-1)
        elif v_wind > 16:
            return 2300
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
        mod_coeff =  4 # chosen arbitrary to get quite good values
        if wind_speed <= 3:
            return 0
        else:
            #return mod_coeff/(wind_speed)
            return 0.8
    # END THRuST_COEFF
    
    def calc_interference(self):
        """
        Calculate interence matrix based on jensens wake model
        entries are the loss of power at each node relative to node i
        """
        print("Creating interference matrix ...")
        interferenz_matrix = np.zeros(shape=(self.n, self.axx, self.axy))
        pos = 0
        for wind_direction in self.w_dirs:
            for i in self.setV:
                infer  = []
                for j in self.setV:
                    alpha   = np.matmul(np.transpose(self.grid[j] - self.grid[i]), wind_direction) / np.linalg.norm(wind_direction)**2
                    if alpha < 0:
                        infer.append(0)
                    else:
                        P_prime = self.grid[i] + alpha*wind_direction
                        delta_1 = self.dist(self.grid[i], P_prime)
                        delta_2 = self.dist(self.grid[j], P_prime)
                        if delta_2 <= (self.D+2*self.k*delta_1)/2:
                            infer.append( self.P(self.v_wind[pos]) - self.P(self.v_wind[pos] - self.jensen_mod(self.v_wind[pos], self.D, delta_1)) )
                        else:
                            infer.append(0)
                interferenz_matrix[i] += (1/self.w_dirs.shape[0])*(1/self.w_dirs.shape[0])*np.array(infer).reshape((self.axx, self.axy)).T
            pos += 1

        filename = "interference_matrix_" + str(self.axx) + "_" + str(self.axy) + "_" + str(self.n) + "_" + str(self.WindInd) +".npy" 
        print("Save interference matrix to " + filename)
        with open(os.path.join(self.dir_sol,filename), "wb") as infer_file:
            np.save(infer_file, interferenz_matrix)
  
        self.init_interference_matrix()
        print("Done!") 
    # END CALC_INTERFERENCE
    
    def find_initial_sol(self, Nmin, Nmax):
        """
        Search for an unconstrained solution to warmstart the solver
        """
        print("Doing init with warmstart")
        OP_initial = cpx.Model(name="Wind Farm Layout Warmstart", log_output=True)
        x_vars = {(i): OP_initial.binary_var(name="x_{0}".format(i)) for i in self.setV}
        OP_initial.add_constraint(OP_initial.sum(x_vars[i] for i in self.setV) <= Nmax)
        OP_initial.add_constraint(OP_initial.sum(x_vars[i] for i in self.setV) >= Nmin)
        for cnt in range(0,self.n):
            for el in self.setE[cnt]:
                OP_initial.add_constraint(x_vars[cnt] + x_vars[el] <= 1)
        obj = OP_initial.sum(self.Pi*x_vars[i] for i in self.setV)
        OP_initial.maximize(obj)
        OP_initial.solve()
        self.initial_sol = [OP_initial.solution.get_value(x_vars[element]) for element in self.setV]
        
        initial_sol_indices = []
        cnt = 0
        for element in self.setV:
            if OP_initial.solution.get_value(x_vars[element]) == 1:
                initial_sol_indices.append(cnt)
            cnt += 1

        # save this
        filename = "wfl_initial_sol_" + str(self.axx) + "_" + str(self.axy) + ".npy"
        with open(os.path.join(self.dir_sol, filename), "wb") as sol_file:
            np.save(sol_file, initial_sol_indices)
    # END FIND_INITIAL_SOL

    def save_sol(self):
        """
        Save solution of optimizaton method to file
        """
        filename = "../solutions/wfl_sol_dual_" + str(self.axx) + "_" + str(self.axy) + "_" + str(self.WindInd) +".npy"
        with open(os.path.join(self.dir_sol, filename), "wb") as sol_file:
            np.save(sol_file, self.sol)
            np.save(sol_file, self.sol_indices)
            np.save(sol_file, self.sol_inf)
            np.save(sol_file, self.sol_obj)
            np.save(sol_file, self.Pi)
        print("Done!")
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
        self.load_layout_sol(1)
        
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
            np.save(cr_sol, self.setV0)
        print("Done!")
    # END SAVE_SOL

# END OF CABLE ROUTING OPTIMIZATION
# ###################################

