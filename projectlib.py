import numpy as np

class Settings:
    """
    A class representing the settings for a simple opinion model simulation.
    """
    def __init__(self, duration=500, alpha_pro=0.25, alpha_neg=0.25, lambda_dist=1, beta_update=0.01, beta_spread=0.01, gamma_extr=0.005, 
                 N_l=None, N_n=500, 
                 G_l=[[(0,0),(0,1),(1,0),(1,1)],[(0,0.5),(0.5,0),(1,0.5),(0.5,1)]], G_n=[4,4]):

        if (not isinstance(duration, int)):
            raise ValueError(
                "Duration is not valid"
            )
        self.duration = duration
        self.alpha_pro = alpha_pro
        self.alpha_neg = alpha_neg
        self.lambda_dist = lambda_dist
        self.beta_update = beta_update
        self.beta_spread = beta_spread
        self.gamma_extr = gamma_extr

        #
        # Initialise individuals and their locations in N_l
        """
        if (N_l is None):
            N_l = []
            # Randomly assign activity positions in each activity period
            for j in range(1,N_n+1):
                Ind = Individual(alpha_neg, alpha_pro)
                #(x, y) = np.random.uniform(0,1,2)
                N_l.append(tuple([Ind.x, Ind.y]))
        """
        self.N_n = N_n
        self.N_l = N_l

        #
        # Initialise the activity locations if G_l has been specified as None
        """
        if (G_l is None):
            G_l = []
            # Randomly assign activity positions in each activity period
            for j in G_n:
                aList = []
                for k in range(1,j+1):
                    (x, y) = np.random.uniform(0,1,2)
                    aList.append(tuple([x,y]))
                G_l.append(aList)
        """
        self.G_n = G_n
        self.G_l = G_l
            

class Individual:

    def __init__(self, alpha_neg, alpha_pos):

        #
        # Each individual has
        # an id - initially None this will later be set to the graph node id
        # a type - unbiased / negative / positive
        # an opinion - dependent on their type
        # a location - dependent on their type

        self.id = None

        self.type = np.random.choice(["neg","pos","unbiased"], 1, True, [alpha_neg, alpha_pos, 1-alpha_neg-alpha_pos])[0]
        
        if (self.type == "unbiased"):
            self.opinion = np.random.uniform()
            (self.x, self.y) = np.random.uniform(0, 1, 2)
        
        elif (self.type == "neg"):
            self.opinion = np.random.uniform(0, 0.25)
            (self.x, self.y) = np.random.uniform(0, 0.25, 2)
        
        else:  # (self.type == "pos")
            self.opinion = np.random.uniform(0.75, 1)
            (self.x, self.y) = np.random.uniform(0.75, 1, 2)

class Activity:

    def __init__(self, period, x=None, y=None):

        #
        # Each activity has
        # an id - initially None this will later be set to the graph node id
        # a period - this is synonymous with group in the multipartite graph
        # a location 

        self.id = None

        self.period = period
        
        # If location coordinates are not given randomly assign them
        # otherwise use the coordinates
        if (x is None or y is None):
            (self.x, self.y) = np.random.uniform(0, 1, 2)
        else:
            self.x = x
            self.y = y
        

"""
IGNORE
class Individuals:

    def __init__(self, N_l, N_n=500):

        if (not (N_l is None) and len(N_l) != N_n):
            raise ValueError(
                "List of locations is not valid"
            )
        self.N_n = N_n
        if (N_l is None):
            # generate the list of individuals locations
            self.N_l = N_l
"""