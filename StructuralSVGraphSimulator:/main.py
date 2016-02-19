#import collections
# GenomeGraphSimulation (master) conda create -n compbio python=3 anaconda
import copy
import matplotlib.pyplot as plt
import numpy as np
import networkx as nx
import pprint
import os

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
        return False

    def is_before(self, other):
        """ returns true if self is finishes before other"""
        if self.end < other.start:
            return True
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
        simple = True
        if simple:
            self.longer_seq = 12
            self.deletion = Fragment(6, 6)
            self.coverage = 20
            self.poisson_lambda = 1
            self.bases_per_fragment = 3
            self.longer_fragments = []
            self.shorter_fragments = []
            self.sources = []
            self.targets = []
            self.adj_mat = None
            self.nx_graph = None
        else:
            self.longer_seq = 600
            self.deletion = Fragment(200, 400)
            self.coverage = 30
            self.poisson_lambda = 10
            self.bases_per_fragment = 50
            self.longer_fragments = []
            self.shorter_fragments = []
            self.sources = []
            self.targets = []
            self.adj_mat = None
            self.nx_graph = None

        self.generate_fragments()
        self.make_networkx_graph()

    def generate_fragments(self):
        """ main argument"""
        self.generate_fragments_for_longer_sequence()
        self.create_shorter_fragments()
        self.find_sources_targets()

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

    def create_shorter_fragments(self):
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

    def find_sources_targets(self):
        """populates self.sources and self.targets with the
        fragments are before and after the deletion"""
        li = copy.deepcopy(self.longer_fragments)
        not_intersecting_deletion = lambda x : not x.intersect(self.deletion)
        self.sources = []
        self.targets = []
        for frag in filter(not_intersecting_deletion, li):
            if frag.is_before(self.deletion):
                self.sources.append(frag)
            else:
                self.targets.append(frag)

    def make_networkx_graph(self):
        # TODO: there's a bug here where it won't add two nodes
        # for two fragments that start in the same place
        if self.adj_mat is None: self.generate_adjacency_matrix()
        self.nx_graph = nx.Graph(self.adj_mat)

    def nx_density(self):
        if self.nx_graph is None:
            return -1
        return nx.density(self.nx_graph)

    def print_edgelist(self):
        if self.nx_graph is None:
            self.make_networkx_graph()
        for line in nx.generate_edgelist(self.nx_graph, data=False):
            print(line)

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
        for i,u in enumerate(self.longer_fragments):
            u_name = str(i) + "_" + str(u)
            for j,v in enumerate(self.longer_fragments):
                v_name = str(j) + "_" + str(v)
                if u.intersect(v):
                    if u in self.adj_mat:
                        self.adj_mat[u_name].append(v_name)
                    else:
                        self.adj_mat[u_name] = [v_name]

    def find_k_shortest_paths(self, k):
        """ finds k shortest paths using PathLinker.
        This uses a roundabout, verbose method by writing to files,
        and calling PathLinker via the os.
        Assumes all edges have weight 1"""

        if self.nx_graph is None:
            self.make_networkx_graph()

        fn_net = "tmp_net.txt" # a file specifying the network structre
        fn_src_tgt = "tmp_src_tgt.txt" # a file specifying the sources and targets

        with open(fn_net, 'w') as f:
            f.write("#Node1\tNode2\n")
            for line in nx.generate_edgelist(self.nx_graph, data=False):
                nodes = line.strip().split()
                f.write(str(nodes[0]) + "\t" + str(nodes[1]) + "\t1")
                f.write("\n")

        with open(fn_src_tgt, 'w') as f:
            f.write("#Node Node type\n")
            for src in self.sources:
                f.write(str(src) + "\tsource")
                f.write("\n")
            for tgt in self.targets:
                f.write(str(tgt) + "\ttarget")
                f.write("\n")

        # run PathLinker
        s = "python ../PathLinker/PathLinker.py " + fn_net + " " + fn_src_tgt
        os.system(s)

    def print_fragments(self):
        print("longer fragments:")
        pprint.pprint(self.longer_fragments)
        print("shorter fragments:")
        pprint.pprint(self.shorter_fragments)
        print("source fragments:")
        pprint.pprint(self.sources)
        print("target fragments:")
        pprint.pprint(self.targets)

    def __str__(self):
        return "longer: " + str(len(self.longer_fragments)) + \
                "\nshorter: " + str(len(self.shorter_fragments)) + \
                "\nsources: " + str(len(self.sources)) + \
                "\ntargets: " + str(len(self.targets))

    def __repr__(self):
        return str(self)


G = StructuralSVGraphSimulator()
print(G)
G.find_k_shortest_paths(100)

