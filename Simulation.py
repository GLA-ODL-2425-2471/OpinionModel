import numpy as np
import pandas as pd
from base import Multipartite
from projectlib import Settings
from projectlib import Individual
from projectlib import Activity
#import projectlib as pjl

class Simulation:
    """
    A class representing a simple opinion model simulation.
    """
    def __init__(self, settings):
        #
        # Parameter settings
        self.duration = settings.duration
        self.alpha_pro = settings.alpha_pro
        self.alpha_neg = settings.alpha_neg
        self.lambda_dist = settings.lambda_dist
        self.beta_update = settings.beta_update
        self.beta_spread = settings.beta_spread
        self.gamma_extr = settings.gamma_extr

        #
        # Individual settings
        self.N_n = settings.N_n
        if (settings.N_l is None):
            self.N_l = []
            self.N_l_supplied = False
        else:
            self.N_l = settings.N_l
            self.N_l_supplied = True

        #
        # Activity settings
        self.G_n = settings.G_n
        if (settings.G_l is None):
            self.G_l = []
            self.G_l_supplied = False
        else:
            self.G_l = settings.G_l
            self.G_l_supplied = True

        #
        # A flag to determine if simulation has been run. Initially False
        self._has_simulation_run = False

        # The graph will hold the relationships between Individuals and Activities
        # but data structures are also necessary to hold information regarding 
        # Individuals and Activities
        # Data structure to hold information about Individuals (node id, type, opinion, location)
        self.individuals = {}
        # Data structure to hold information about Activities (node_id, location)
        self.activities = {}

        #
        # Create the multipartite graph
        self.graph = Multipartite(1+len(self.G_n))
        print("graph created")
        self.graph = self.get_multipartite_graph()



    def __str__(self):
        """ Print string representation of the object """
        # return "A simulation with {} groups, {} nodes, and {} edges.".format(len(self.groups.keys()), len(self.graphNodes.keys()), len(self.graphEdges))
        return "A simulation with {} Individuals. {} Activities. Has run? {}".format(self.N_n, len(self.G_n), self._has_simulation_run) 
    
    """
    Helper function to compute Eucliden distance between locations
    """
    def _distance(self, indX, indY, activityX, activityY):
        #
        # May need some error checking????
        #
        return (np.sqrt( (indX-activityX)**2 + (indY-activityY)**2))
    
    def _compute_activity_probabilities(self, x, y, group):

        nodes = self.graph.get_nodes(group)

        probabilities = []
        distances = {}
        denom = 0
        for node in nodes:
            dist = self._distance(x, y, self.activities.loc[node].x, self.activities.loc[node].y)
            distances = { **distances, node : dist}
            denom += np.exp(-self.lambda_dist * dist)

        for node in nodes:
            prob = np.exp(-self.lambda_dist * distances[node]) / denom
            probabilities.append((node, prob))

        return probabilities

    def get_opinion(self, time):

        #
        # Missing the time element for now
        return pd.DataFrame(self.individuals, columns=["opinion"])
        
    def get_multipartite_graph(self):
        

        # Create individuals, add them to the graph and supporting data structures
        id_list=[]
        type_list=[]
        opinion_list=[]
        x_list=[]
        y_list=[]
        for j in range(self.N_n):
            # 'Make' an Individual
            person = Individual(self.alpha_neg, self.alpha_pro)
            
            # Add the Individual to the graph in group 0 as per spec
            person.id = self.graph.add_node(0)
            
            # Add the Individual's location to the N_l list
            id_list.append(person.id)
            type_list.append(person.type)
            opinion_list.append(person.opinion)
            if (self.N_l_supplied):
                # Add the Individual to the dictionary of individuals identified by their node_id
                details = {"type" : person.type, "opinion" : person.opinion, "x" : self.N_l[j][0], "y" : self.N_l[j][1]}
                x_list.append(self.N_l[j][0])
                y_list.append(self.N_l[j][1])
            else:
                self.N_l.append(tuple([person.x, person.y]))
                # Add the Individual to the dictionary of individuals identified by their node_id
                details = {"type" : person.type, "opinion" : person.opinion, "x" : person.x, "y" : person.y}
                x_list.append(person.x)
                y_list.append(person.y)

            self.individuals = {**self.individuals, person.id : details}
            #self.individuals = {**self.individuals, "node_id" : person.id, "type" : person.type, "opinion" : person.opinion, "x" : person.x, "y" : person.y}
        
        individual_dict = {"id" : id_list, "type" : type_list, "opinion" : opinion_list, "x" : x_list, "y" : y_list}
        print(individual_dict) 
        individual_df = pd.DataFrame.from_dict(individual_dict)
        print(individual_df)

        # Create the Activities, add them to the graph and supporting data structures
        #
        # Initialise the activity locations if G_l has been specified as None

        for i in range(len(self.G_n)):

            aList = []
            # Randomly assign activity positions in each activity period
            for j in range(self.G_n[i]):
                    
                # 'Make' an Activity.
                # Location was not provided - it wil be generated randomly
                activity = Activity(i+1)

                # Add the Activity to the graph with group = period
                activity.id = self.graph.add_node(i+1)

                if (self.G_l_supplied):
                    # Add the Activity to the dictionary of activities identified by their node_id
                    details = {"period" : activity.period, "x" : self.G_l[i][j][0], "y" : self.G_l[i][j][1]}
                else:
                    # Add the Activity's location to the list of activity locations for the period
                    aList.append(tuple([activity.x, activity.y]))
                    # Add the Activity to the dictionary of activities identified by their node_id
                    details = {"period" : activity.period, "x" : activity.x, "y" : activity.y}

                self.activities = {**self.activities, activity.id : details}
                #self.activities = {**self.activities, "node_id" : activity.id, "period" : activity.period, "x" : activity.x, "y" : activity.y}

            if (not self.G_l_supplied):
                self.G_l.append(aList)
        
        # Convert the dictionaries to panda data frames and transpose them
        self.individuals = (pd.DataFrame(self.individuals)).T
        self.activities = (pd.DataFrame(self.activities)).T


        """
        print("----------------------")
        print("Individuals")
        print(self.individuals)
        print("----------------------")

        print("----------------------")
        print("Activities")
        print(self.activities)
        print("----------------------")
        """

        for person_id in list(self.individuals.index.values):
        #for person_id in self.individuals.keys():
            #person = self.individuals.loc[person_id]
            #person = self.individuals.get(person_id)
            for group in range(1,len(self.G_n)+1):
                probability_list = self._compute_activity_probabilities(person.x, person.y, group)
                #probability_list = self._compute_activity_probabilities(person["x"], person["y"], group)
                #print("Group {} Probabilities {}".format(group, probability_list))
                choice = np.random.uniform()
                run_sum = 0
                cum_prob = [(i, (run_sum := run_sum + j)) for (i,j) in probability_list]
                #cum_prob = ((run_sum := run_sum + j) for (i,j) in probability_list)
                activity_id = [ i for (i, prob) in cum_prob if prob > choice]
                #print("Cumulative probability {}. Random choice {}. Selected Activity {}.".format(cum_prob, choice, activity_id[0]))
                #print("Individual {}. Selected Activity {}.".format(person_id, activity_id[0]))
                #
                # Now create edges between individuals and activities based on returned probabilities
                #if self.graph.node_exists(int(person_id)):
                #    print("Individual node exists {}. Group {}".format(person_id, self.graph.get_group(int(person_id))))
                #if self.graph.node_exists(activity_id[0]):
                #    print("Activity node exists {}. Group {}".format(activity_id[0], self.graph.get_group(activity_id[0])))
                #print("Activity node esists {}.").format(self.graph.node_exists(activity_id))
                self.graph.add_multipartite_edge(int(person_id), activity_id[0])

        #for person in individual_df.iterrows():
        #    if self.graph.node_exists(person["id"]):
        #        print("Individual node exists {}. Group {}".format(person["id"], self.graph.get_group(person["id"])))

        [ print(self.graph.node_exists(person)) for person in individual_df['id'] ]
        #print(s)


        #print(self.activities.loc[7,["x"]])
        #print(self.activities.loc[7]["x"])
        #
        # Compute activity probability for each individual
        # use that to create edges between individuals and activities

        return(self.graph)

    def _perform_opinion_update(self, opinions):
        pass

    def _perform_opinion_activity(self, activity):
        pass

    def run(self):
        #
        # run the simulation

        #
        # Flag the simulation has run
        self._has_simulation_run = True

    def plot_network(self):
        pass

    def chart(self):
        pass

    def most_polarised(self):
        pass

    def most_polarised_day(self):
        pass

    def activity_summary(self, t):
        pass

    def individual_summary(self):
        pass

    """ needs to be a class or static method"""
    def ensemble_statistics(self):
        pass
    
    def get_elbow_plot(self):
        pass

    def kmeans_clustering(self):
        pass

    def fit_regression_model(self):
        pass


"""
settings = Settings(100, G_l=None, N_n=5)
s = Simulation(settings)
print(s)
#print(s.get_opinion(0))

settings = Settings(100, N_n=5)
s = Simulation(settings)
print(s)
"""

settings = Settings(100, N_l=[(0,0),(0.2,0.2),(0.4,0.4),(0.6,0.6),(0.8,0.8)], N_n=5)
s = Simulation(settings)
print(s)
