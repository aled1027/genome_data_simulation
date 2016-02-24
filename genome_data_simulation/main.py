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

    def __hash__(self):
        return int(self.start) * 10000 + self.end

    def __eq__(self, other):
        return self.start == other.start and self.end == other.end

    def __len__(self):
        return self.end - self.start

    def __str__(self):
        return str(self.start) + ":" + str(self.end)
#        return "Fragment(start={}, end={})".format(self.start, self.end)

    def __repr__(self):
        return str(self)

class StructuralSVGraphSimulator:
    def __init__(self):
        simple = True
        if simple:
            print("using simple parameters")
            self.longer_seq = 12
            self.deletion = Fragment(6, 7)
            self.coverage = 3
            self.poisson_lambda = 1
            self.bases_per_fragment = 3
            self.longer_fragments = []
            self.shorter_fragments = []
            self.sources = []
            self.targets = []
            self.adj_mat = None
            self.nx_graph = None
        else:
            print("using real parameters")
            self.longer_seq = 600
            self.deletion = Fragment(200, 400)
            self.coverage = 10
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
        self.find_sources_targets()

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
        for i, frag in enumerate(li):
            if frag.intersect(self.deletion):
                continue
            frag_name = str(i) + "_" + str(frag)
            if frag.is_before(self.deletion):
                self.sources.append(frag_name)
            else:
                self.targets.append(frag_name)

        if not self.sources:
            print("No source nodes")
        if not self.targets:
            print("No target nodes")


    def make_networkx_graph(self):
        """ based on self.longer_fragments and self.shorter_fragments,
        make a networkx graph"""
        self.nx_graph = nx.Graph()
        all_nodes = list(set(self.longer_fragments).union(set(self.shorter_fragments)))
        for i,u in enumerate(all_nodes):
            u_name = str(i) + "_" + str(u)

            f = None
            if u in self.longer_fragments and u in self.shorter_fragments:
                f = "both"
            elif u in self.shorter_fragments:
                f = "shorter"
            elif u in self.longer_fragments:
                f = "longer"
            else:
                raise RuntimeError("messed up")
            self.nx_graph.add_node(u_name, frag=f)

            for j,v in enumerate(all_nodes):
                if u.intersect(v):
                    v_name = str(j) + "_" + str(v)
                    self.nx_graph.add_edge(u_name, v_name, frag=f)

    def nx_density(self):
        assert(self.nx_graph is not None)
        return nx.density(self.nx_graph)

    def draw_nx(self, filename):
        print("drawing")
        assert(self.nx_graph is not None)

        # settings
        node_size = 450 #default is 300
        font_size = 4.5

        pos = nx.spring_layout(self.nx_graph)
        #pos = nx.spectral_layout(self.nx_graph)

        blue_nodes = [n for n,d in self.nx_graph.nodes_iter(data=True) if d['frag'] == 'longer']
        red_nodes = [n for n,d in self.nx_graph.nodes_iter(data=True) if d['frag'] == 'shorter']
        purple_nodes = [n for n,d in self.nx_graph.nodes_iter(data=True) if d['frag'] == 'both']
        nx.draw_networkx_nodes(self.nx_graph, pos, nodelist=blue_nodes, \
                node_size=node_size, node_color='b')
        nx.draw_networkx_nodes(self.nx_graph, pos, nodelist=red_nodes, \
                node_size=node_size, node_color='r')
        nx.draw_networkx_nodes(self.nx_graph, pos, nodelist=purple_nodes, \
                node_size=node_size, node_color='#aa8cc5')

        nx.draw_networkx_edges(self.nx_graph, pos)
        nx.draw_networkx_labels(self.nx_graph, pos, font_size=font_size)

        # For adding legend -- covering up some nodes
        # import matplotlib.patches as mpatches
        #blue_patch = mpatches.Patch(color='blue', label='only longer')
        #red_patch = mpatches.Patch(color='red', label='only shorter')
        #purple_patch = mpatches.Patch(color='purple', label='both')
        #plt.legend(handles=[red_patch, blue_patch, purple_patch])

        plt.title("red=shorter, blue=longer, purple=both")
        plt.savefig(filename)

    def find_k_shortest_paths(self, k):
        """ finds k shortest paths using PathLinker.
        This uses a roundabout, verbose method by writing to files,
        and calling PathLinker via the os.
        Assumes all edges have weight 1"""

        assert(self.nx_graph is not None and "should be made by now")

        fn_net = "tmp_net.txt" # a file specifying the network structre
        fn_src_tgt = "tmp_src_tgt.txt" # a file specifying the sources and targets


        with open(fn_net, 'w') as f:
            f.write("#Node1\tNode2\n")
            for line in nx.generate_edgelist(self.nx_graph, data=False):
                nodes = line.strip().split()
                f.write(str(nodes[0]) + "\t" + str(nodes[1]) + "\t1")
                f.write("\n")

        # generate sources and targets by grabing all nodes which have attribute
        # "both"

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

    def print_edgelist(self):
        if self.nx_graph is None:
            self.make_networkx_graph()
        for line in nx.generate_edgelist(self.nx_graph, data=False):
            print(line)

    def print_fragments(self):
        print("longer fragments:")
        pprint.pprint(self.longer_fragments)
        print("shorter fragments:")
        pprint.pprint(self.shorter_fragments)
        print("source fragments:")
        pprint.pprint(self.sources)
        print("target fragments:")
        pprint.pprint(self.targets)
        print("deletion: ", self.deletion)

    def print_adj_mat(self):
        assert(self.adj_mat is not None)
        pprint.pprint(self.adj_mat)

    def __str__(self):
        return "longer: " + str(len(self.longer_fragments)) + \
                "\nshorter: " + str(len(self.shorter_fragments)) + \
                "\nsources: " + str(len(self.sources)) + \
                "\ntargets: " + str(len(self.targets))

    def __repr__(self):
        return str(self)


G = StructuralSVGraphSimulator()
print(G)
print(G.nx_density())
G.find_k_shortest_paths(100)
G.draw_nx("graph.pdf")
#G.find_k_shortest_paths(100)

