class Multipartite:
    """
    A class representing a simple multipartite graph 
    with no weights and directions.
    """

    def __init__(self, number_of_groups):
        """Creates an multipartite graph

        Graph is created with a fixed number of groups/collections 
        of independent sets.

        If the parameter is set to 2 this is equivalent to a bipartite graph,
        if 3 it is a tripartite graph.

        Args:
            number_of_groups: Number of groups/independent sets

        """
        # Test that the number of groups is reasonable
        # Chose to exclude graphs with 1 group as then no edges are possible
        # but could allow this.
        if not isinstance(number_of_groups, int) or number_of_groups < 1:
            raise ValueError("Not a valid number of groups")
        self.number_of_groups = number_of_groups
        # graph datastructure constructed using a dictionary for quick lookups
        # maps node_identifier -> set of connections
        # dictionary is used for quick/easy deletions,
        # set is used for quick lookups and filtering
        self.graph = {}
        self.node2group = {}  # Dictionary to map a node to its groups
        self.next_node_id = 0  # next node id for node creation
        self._num_edges = 0  # a counter of the number of edges

    def get_number_of_groups(self):
        """Returns the number of groups in the graph."""
        return self.number_of_groups

    def _check_group(self, group):
        """Internal function for checking group is valid

        Function confirms that the group is of valid type (i.e. an int)
        and is within valid values for the group.

        Arg:
            group: Group to check

        """
        if (not isinstance(group, int)) or (group < 0) or group >=self.number_of_groups:
            raise ValueError(
                "Group is not valid"
            )  # Generic exception could be specialised

    def _check_node(self, node_identifier):
        """Internal function for checking a node identifier is valid

        Function confirms that the node identifier is of valid type (i.e. an int)
        and is a valid node identifiers (i.e. the node has been assigned).

        Arg:
            node_identifier: Group

        """
        # can just use the node exists method, and raise an exception if required
        if not self.node_exists(node_identifier):
            raise ValueError("Node not known")

    def add_node(self, group):
        """Adds a new node to the graph in the given group

        Args:
             group: The group which the new node is added

        Returns:
             Node identifer correcting to the new node. Note this ID can vary between
             different implementations of this class.

         Raises:
             ValueError: If the provided group is not valid
        """
        self._check_group(group)  # use check group to confirm the group is valid
        new_node = self.next_node_id  # assign the next_node_id
        self.node2group[new_node] = group  # assign the group in the relevant dict
        self.graph[new_node] = set()  # add set to store edges
        self.next_node_id += 1  # increase node id for next node
        return new_node

    def node_exists(self, node_identifier):
        """Tests if a node exists in the graph

        Args:
             node_identifier: The node identifier to be checked

        Returns:
             Boolean (True/False) indicating the presence of the node.
        """
        # if identifier is not an int (node_identifer format) it cannot be in the graph
        # This gives desired behaviour if we attempt to check for non-hashable data
        if not isinstance(node_identifier, int):
            return False
        return node_identifier in self.node2group

    def get_nodes(self, group=None):
        """Returns a list of nodes in the full/group of the graph.

        If the optional argument group is set, the method will only
        return nodes in this group.

        Note, changing the output of this method does not change
        the nodes that are present in the graph.

        Args:
            group: Group to limit the returned nodes (Optional)

        Returns:
            A list of nodes that are in the group

        Raises:
             ValueError: If the provided group is not valid
        """
        if group is None:  # is group is node return all notes
            return list(self.node2group.keys())
        else:
            self._check_group(group)  # confirm group is valid
            # filter nodes in this group
            return [x for x, y in self.node2group.items() if y == group]

    def num_nodes(self, group=None):
        """Returns number of nodes in the full/group of the graph.

        If the optional argument group is set, the method will only
        count nodes in this group.

        Args:
            group: Group to limit the returned nodes (Optional)

        Returns:
            Number of nodes in the graph/requested group as an integer.

        Raises:
             ValueError: If the provided group is not valid
        """
        if group is None:  # if group is none can use the size of node2group
            return len(self.node2group)
        else:
            self._check_group(group)  # confirm group is valid
            # count number of nodes with this group
            # could take len(self.get_nodes(group)) for code simplicy over efficiency
            return sum(1 for y in self.node2group.values() if y == group)

    def get_group(self, node_identifier):
        """Returns number of nodes in the full/group of the graph.

        If the optional argument group is set, the method will only
        count nodes in this group.

        Args:
            group: Group to limit the returned nodes (Optional)

        Returns:
            Number of nodes in the graph/requested group as an integer.

        Raises:
             ValueError: If the provided group is not valid
        """
        self._check_node(node_identifier)
        return self.node2group[node_identifier]

    def _add_edge_checks(self, node_identifier1, node_identifier2):
        """Checks that an edge can be added

        Method takes two node identifiers and checks if:
        - Node identifiers are valid
        - Nodes are from different groups
        - Edge does not already exist

        Raises and exception (ValueError) any of the checks fails.

        Implemented as a separate method to help with inheritance for directed/other graph constructions

        Args:
             node_identifier1: node identifier corresponds to the first node in 
                               the potential edge
             node_identifier2: node identifier corresponds to the second node in 
                               the potential edge

        Raises:
             ValueError: If the edge cannot be added to the graph
        """
        self._check_node(node_identifier1)
        self._check_node(node_identifier2)
        g1 = self.get_group(node_identifier1)
        g2 = self.get_group(node_identifier2)
        if g1 == g2:
            raise ValueError("not multipartite")
        # use has edge to tell if the edge already exists
        if self.has_edge(node_identifier1, node_identifier2):
            raise ValueError("edge already exists")

    def add_multipartite_edge(self, node_identifier1, node_identifier2):
        """Add a multipartite edge to the graph

        Method takes two node identifiers, checks if an edge between them
        can be added to the graph. If it cannot an exception is raised.
        If the edge can be added, it is added to graph.

        Args:
             node_identifier1: node identifier corresponds to the first node in the potential edge
             node_identifier2: node identifier corresponds to the second node in the potential edge

        Raises:
             ValueError: If the edge cannot be added to the graph
        """
        # check if the edge can be added
        self._add_edge_checks(node_identifier1, node_identifier2)
        # add other node to the set of connections for each node
        self.graph[node_identifier1].add(node_identifier2)
        self.graph[node_identifier2].add(node_identifier1)
        self._num_edges += 1  # increase the edge counter

    def remove_multipartite_edge(self, node_identifier1, node_identifier2):
        """Removes a multipartite edge to the graph

        Method takes two node identifiers, checks if an edge between them
        can be removed from the graph. If it can, this operation is performed, if it cannot,
        an exception is raised.

        Args:
             node_identifier1: node identifier corresponds to the first node in the potential edge
             node_identifier2: node identifier corresponds to the second node in the potential edge

        Raises:
             ValueError: If the edge cannot be added to the graph
        """
        self._check_node(node_identifier1)  # check that each node exists
        self._check_node(node_identifier2)
        if not self.has_edge(
            node_identifier1, node_identifier2
        ):  # confirm the edge exists
            raise ValueError("There is no edge between the provided nodes")
        # remove other node to the set of connections for each node
        self.graph[node_identifier1].remove(node_identifier2)
        self.graph[node_identifier2].remove(node_identifier1)
        self._num_edges -= 1

    def get_edges(self):
        """Returns a list of all edges in the graph

        Returns:
            A list of all edges in the graph. Edge edge is represented by a tuple of (node_id1,node_id2).

        """
        # create list using a comprehension as each edge will exist twice
        # (in the set for each its constituent nodes) we filter to only include cases with x<y
        return [(x, y) for x in self.graph for y in self.graph[x] if x < y]

    def get_edges_node(self, node_identifier):
        """Returns a list of edges in the graph which involve a specified node

        Args:
             node_identifier: Identifier for the node that must be involved in all returned edges

        Returns:
            A list of all edges involving the given node. Edge edge is represented by a tuple of (node_id1,node_id2).

        Raises:
            ValueError: If the node_identifier is not valid

        """
        self._check_node(node_identifier)  # check that the node is valid
        # Construct edge list from the graph datastructure
        # only return edges (x,y) where x<y to give consistent behaviour to get_edges
        return [
            (node_identifier, y) if node_identifier < y else (y, node_identifier)
            for y in self.graph[node_identifier]
        ]

    def get_edges_nodes(self, collection_of_node_identifiers, exclusive=False):
        """Returns a list of edges in the graph which involve a collection of nodes

        Args:
             collection_of_node_identifiers: A collection (iterable) of node identifier which must be involved in all returned edges
             exclusive: A boolean flag, if set to True it restricts returned edges to only those with both nodes in the supplied collection

        Returns:
            A list of all edges involving the given nodes. Edge edge is represented by a tuple of (node_id1,node_id2).

        Raises:
            ValueError: If the node_identifier is not valid, or the exclusive flag is not True/False 

        """
        if not ((exclusive is True) or (exclusive is False)): 
            raise ValueError("Exclusive is not True or False")
        # make a list in case the provided collection is a generator/iterator
        collection_of_node_identifiers = list(collection_of_node_identifiers)
        # check that each node is valid
        for node in collection_of_node_identifiers:
            self._check_node(node)
        if exclusive:
            # form a set to use efficient set methods
            s1 = set(collection_of_node_identifiers)  
            return [(x, y) for x in s1 for y in s1.intersection(self.graph[x]) if x < y]
        else:
            # Construct edge list from the graph datastructure
            # only return edges (x,y) where x<y to give consistent behaviour to get_edges
            res = [
                (x, y) if x < y else (y, x)
                for x in collection_of_node_identifiers
                for y in self.graph[x]
            ]
            return sorted(list(set(res)))

    def get_edges_group(self, group1=None, group2=None):
        """Returns a list of edges in the graph that involve upto 2 groups

        This method returns all edges that involve the groups given.

        Depending on the number of groups defined the method will return a different set of edges. This can be specified as follows:

        - Group1 and Group2 are set as None: Returns all edges in the graph
        - Only Group1/Group2 is set to a valid group, other is None: Returns all edges involving nodes in this group to any other node in the group
        - Group1 and Group2 is set to valid groups: Returns all edges that involve both groups.

        Args:
            group1: Group to limit the returned edges (Optional)
            group2: Group to limit the returned nodes (Optional)

        Returns:
            A list of all edges involving the given groups (see docstring for details). Edge edge is represented by a tuple of (node_id1,node_id2).

        Raises:
            ValueError: If the groups are not valid

        """
        if group1 is None and group2 is None:
            return self.get_edges()  # if both None can use get_edges method
        elif group1 is None or group2 is None:
            g1 = group1 if group2 is None else group2  # extract group that is not None
            self._check_group(g1)  # check this group
            tempEdges = [self.get_edges_node(x) for x in self.get_nodes(group=g1)]
            # filter for edges where x<y to give unique edges
            return [(x, y) for s1 in tempEdges for x, y in s1 if x < y]
        else:
            # test if nodes in edge below to the correct groups
            self._check_group(group1)
            self._check_group(group2)
            edges = [
                (x, y) if x < y else (y, x)
                for x in self.get_nodes(group=group1)
                for y in self.graph[x]
                if self.node2group[y] == group2
            ]
            return edges

    def has_edge(self, node_identifier1, node_identifier2):
        """Checks if an edge exists in the graph

        Method takes two node identifiers, checks if an edge exists between them.

        Args:
             node_identifier1: node identifier corresponds to the one node in the potential edge
             node_identifier2: node identifier corresponds to the other node in the potential edge

        Returns:
            A boolean indicating if the edge exists.

        Raises:
             ValueError: If the node identifiers are not valid
        """
        self._check_node(node_identifier1)  # check that the nodes are valid
        self._check_node(node_identifier2)
        # check if the node exists
        return node_identifier2 in self.graph[node_identifier1]

    def remove_node(self, node_identifier):
        """Remove a node from the graph.

        Method will remove all information about a node from the graph,
        including its group and all edges associated with it.

        Args:
             node_identifier1: node identifier to be removed

        Raises:
             ValueError: If the node identifier is not valid
        """
        self._check_node(node_identifier)  # check ID is valid
        del self.node2group[node_identifier]  # remove node from group mapping
        for x in self.graph[node_identifier]:
            # remove node from the edge lists of all connected nodes
            self.graph[x].remove(node_identifier)
        # reduce the number of edges in the graph by the corresponding amount
        self._num_edges -= len(self.graph[node_identifier])
        # finally remove the set containing the connections for the node
        del self.graph[node_identifier]

    def num_edges(self):
        """Number of edges in the graph.

        Returns:
            The current number of edges in the graph
        """
        return self._num_edges

    def _number_possible_edges(self):
        """Number of possible edges in the graph.

        Method is mostly used for inheritance with other graph constructions.

        Returns:
            The current number of possible edges in the graph
        """
        count = [
            0,
        ] * self.number_of_groups  # list of number of nodes in each group

        # count number of nodes in each group
        for x in self.node2group.values():
            count[x] += 1
        t1 = sum(count)  # total number of nodes
        # compute the number of possible edges as
        # sum of num_nodes_in_group_i * num_nodes_not_in_group_i
        # divide by 2 due to the undirected nature
        return sum(x * (t1 - x) for x in count) / 2

    def density(self):
        """Density of the graph.

        Density is defined as the percentage of the possible
        edges in the graph that exist.
        As this is a multipartite structure, this method
        only counts edges that are allowed by the structure.

        Returns:
            The density of the graph as a float
        """
        return self.num_edges() / self._number_possible_edges()
