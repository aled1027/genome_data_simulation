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
        :param other: Fragment
        """
        if self.start <= other.start <= self.end or \
                self.start <= other.end <= self.end:
            return True
        else:
            return False

    def is_contained_in(self, other):
        """ Returns True if self is totally inside of other
        :param other: Fragment"""
        if other.start <= self.start and self.end <= other.end:
            return True
        else:
            return False

    def __len__(self):
        return self.end - self.start

    def __str__(self):
        return str(self.start)
#        return "Fragment(start={}, end={})".format(self.start, self.end)

    def __repr__(self):
        return str(self)

class StructuralSVGraphSimulator:
    def __init__(self):
        self.longer_seq = 600
        self.deletion = Fragment(200, 400)
        self.coverage = 30
        self.poisson_lambda = 10
        self.bases_per_fragment = 50
        self.longer_fragments = []
        self.shorter_fragments = []
        self.adj_mat = None
        self.nx_graph = None

    @property
    def nx_density(self):
        if self.nx_graph is None:
            return -1
        return nx.density(self.nx_graph)

    def generate_fragments(self):
        """ main argument"""
        self.generate_fragments_for_longer_sequence()
        self.create_shorter_fragments()

    def print_edgelist(self):
        if self.nx_graph is None:
            self.make_networkx_graph()
        for line in nx.generate_edgelist(self.nx_graph, data=False):
            print(line)

    def make_networkx_graph(self):
        if self.adj_mat is None:
            self.generate_adjacency_matrix()
        self.nx_graph = nx.Graph(self.adj_mat)

    def save_nx_fig(self, filename):
        if self.nx_graph is None:
            self.make_networkx_graph()
        pos = nx.spring_layout(self.nx_graph)
        nx.draw(self.nx_graph, pos)

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
        # TODO: this argument needs to be more robust, where the fragments
        # are adjusted based on the deletion, not just removed

        """creates shorter_fragment based on longer_fragments.
        assumes longer_fragments is already populated"""
        # again, writing this in the most naiive, unpythonic way
        self.shorter_fragments = copy.deepcopy(self.longer_fragments)
        not_contained = lambda x : not x.is_contained_in(self.deletion)
        self.shorter_fragments = list(filter(not_contained, self.shorter_fragments))

        for frag in self.shorter_fragments:
            if frag.intersect(self.deletion):
                if frag.start < self.deletion.start:
                    frag.end = self.deletion.start
                else:
                    frag.start = self.deletion.end

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
        return "longer: " + str(len(self.longer_fragments)) + "\nshorter: " + str(len(self.shorter_fragments))

    def __repr__(self):
        return str(self)

for _ in range(10):
    G = StructuralSVGraphSimulator()
    G.generate_fragments()
    G.make_networkx_graph()
    print(G)
    print(G.nx_density)
    print('')
#G.print_edgelist()
