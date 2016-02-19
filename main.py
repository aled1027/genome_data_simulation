#import collections
# GenomeGraphSimulation (master) conda create -n compbio python=3 anaconda
import copy
import matplotlib.pyplot as plt
import numpy as np
import networkx as nx

class Fragment:
    def __init__(self, start=0, end=0):
        self.start = start
        self.end = end

    def intersect(self, other):
        """ Returns true if self and other intersect
        :param other: a fragment object
        """
        if self.start <= other.start <= self.end or \
                self.start <= other.end <= self.end:
            return True
        else:
            return False

    def __len__(self):
        return self.end - self.start

    def __str__(self):
        return "Fragment(start={}, end={})".format(self.start, self.end)

    def __repr__(self):
        return str(self)

class GraphSimulator:
    def __init__(self):
        self.longer_seq = 10
        self.deletion = Fragment(5, 7)
        self.coverage = 1
        self.poisson_lambda = 2
        self.bases_per_fragment = 2
        self.longer_fragments = []
        self.shorter_fragments = []
        self.adj_mat = {}

    def run(self):
        """ main argument"""
        self.generate_fragments_for_longer_sequence()
        self.create_shorter_fragments()

    def give_me_a_networkx_graph(self):
        self.generate_adjacency_matrix()
        return nx.Graph(self.adj_mat)


    def generate_adjacency_matrix(self):
        """ Using self.longer_fragments and self.shorter_fragments,
        generate a graph. """
        # the worst way this could possibly be written.
        # whatever.
        self.adj_mat = {}
        for u in self.longer_fragments:
            for v in self.longer_fragments:
                if u.intersect(v):
                    if u in self.adj_mat:
                        self.adj_mat[u].append(v)
                    else:
                        self.adj_mat[u] = [v]

    def create_shorter_fragments(self):
        """creates shorter_fragment based on longer_fragments.
        assumes longer_fragments is already populated"""
        self.shorter_fragments = copy.deepcopy(self.longer_fragments)
        no_intersect = lambda x : False if self.deletion.intersect(x) else True
        return filter(no_intersect, self.shorter_fragments)

    def generate_fragments_for_longer_sequence(self):
        for _ in range(self.coverage):
            prev = 0
            while True:
                start = prev + np.random.poisson(self.poisson_lambda)
                end = start + self.bases_per_fragment
                if end <= self.longer_seq:
                    self.longer_fragments.append(Fragment(start, start + self.bases_per_fragment))
                    prev = start
                else:
                    break

    def __str__(self):
        return "longer: " + str(self.longer_fragments) + "\nshorter: " + str(self.shorter_fragments)

    def __repr__(self):
        return str(self)

gs = GraphSimulator()
gs.run()
gr = gs.give_me_a_networkx_graph()
nx.draw(gr)
plt.savefig("fig.png")

