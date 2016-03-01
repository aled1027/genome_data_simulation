#import collections
# GenomeGraphSimulation (master) conda create -n compbio python=3 anaconda
import copy
import matplotlib.pyplot as plt
import numpy as np
import networkx as nx
import pprint
import os

class SimpleFragment:
    """ a simple simple fragment with a start and stop"""
    def __init__(self, start=0, end=0, key=0, string=None):
        if string is None:
            self.start = start
            self.end = end
            self.key = key
        else:
            s = string.strip().split("_")
            self.start = int(s[0])
            self.end = int(s[1])
            self.key = int(s[2])

    def intersect(self, other):
        """ Returns true if self and other intersect
        :param other: SimpleFragment
        """
        if self.start <= other.start <= self.end or \
                self.start <= other.end <= self.end:
            return True
        else:
            return False

    def is_contained_in(self, other):
        """ Returns True if self is totally inside of other
        :param other: SimpleFragment"""
        if other.start <= self.start and self.end <= other.end:
            return True
        return False

    def is_before(self, other):
        """ returns true if self is finishes before other"""
        if self.end <= other.start:
            return True
        return False

    def __hash__(self):
        return int(self.start) * 10000 + self.end

    def __eq__(self, other):
        return self.start == other.start and self.end == other.end

    def __len__(self):
        return self.end - self.start

    def __str__(self):
        return str(self.start) + "_" + str(self.end) + "_" + str(self.key)

    def __repr__(self):
        return str(self)

class Fragment:
    def __init__(self, simple_fragments=None):
        self.simple_fragments = simple_fragments

    def intersect(self, other):
        """
        :param other: Fragment
        """
        return any(s.intersect(o) for s in self.simple_fragments for o in other.simple_fragments)

    def __str__(self):
        s = ""
        for sf in self.simple_fragments:
            s += str(sf) + "\n"
        return s[0:-1]

    def __repr__(self):
        return str(self)

    def is_before(self, other):
        """
        returns true if all self.simple_fragments are before other.simple_fragments else false
        :param other: Fragment
        """
        return all(s.is_before(o) for s in self.simple_fragments for o in other.simple_fragments)

    def is_contained_in(self, other):
        s = sorted(self.simple_fragments, key=lambda f : f.start)
        o = sorted(other.simple_fragments, key=lambda f : f.start)
        return o[0].start <= s[0].start and s[-1].end <= o[-1].end

    def __hash__(self):
        return sum(hash(s)*(100**i) for i,s in enumerate(self.simple_fragments))

    def __len__(self):
        return sum(len(s) for s in self.simple_fragments)

    def __eq__(self, other):
        sf = sorted(self.simple_fragments, key=lambda f : f.start)
        of = sorted(other.simple_fragments, key=lambda f : f.start)
        for s,o in zip(sf, of):
            if s != o:
                return False
        return False

    def starts_before(self, other):
        sf = sorted(self.simple_fragments, key=lambda f : f.start)
        of = sorted(other.simple_fragments, key=lambda f : f.start)
        return sf[0].start < of[0].start


class StructuralSVGraphSimulator:
    def __init__(self):
        simple = False
        if simple:
            print("using simple parameters")
            self.longer_seq = 12
            self.deletion = Fragment([SimpleFragment(6, 7)])
            self.coverage = 3
            self.poisson_lambda = 1
            self.bases_per_fragment = 3
        else:
            print("using real parameters")
            self.longer_seq = 600
            self.deletion = Fragment([SimpleFragment(200, 400)])
            self.coverage = 10
            self.poisson_lambda = 10
            self.bases_per_fragment = 50

        self.longer_fragments = []
        self.shorter_fragments = []
        self.nx_graph = None

        self.generate_fragments()
        self.make_networkx_graph()

    def generate_fragments(self):
        """ main argument"""
        self.generate_fragments_for_longer_sequence()
        self.create_shorter_fragments()

    def generate_fragments_for_longer_sequence(self):
        count = 0
        for _ in range(self.coverage):
            prev = 0
            while True:
                start = prev + np.random.poisson(self.poisson_lambda)
                end = start + self.bases_per_fragment
                if end <= self.longer_seq:
                    self.longer_fragments.append(Fragment([SimpleFragment(start, start + self.bases_per_fragment, count)]))
                    prev = start
                    count += 1
                else:
                    break

    def create_shorter_fragments(self):
        """creates shorter_fragment based on longer_fragments.
        assumes longer_fragments is already populated"""
        # again, writing this in the most naive, unpythonic way
        count = len(self.longer_fragments)
        for _ in range(self.coverage):
            prev = 0
            while True:
                start = prev + np.random.poisson(self.poisson_lambda)
                if self.deletion.simple_fragments[0].start <= start <= self.deletion.simple_fragments[0].end:
                    prev = self.deletion.simple_fragments[0].end + 1
                    continue

                end = start + self.bases_per_fragment

                if end <= self.longer_seq:
                    self.shorter_fragments.append(Fragment([SimpleFragment(start, start + self.bases_per_fragment, count)]))
                    prev = start
                    count += 1
                else:
                    break



    def make_networkx_graph(self):
        """ based on self.longer_fragments and self.shorter_fragments,
        make a networkx graph"""

        """
        TODO work on this.
        our nodes are really Fragments (not SimpleFragments)
        labels should be: before, after, long, short (as opposed to current labels)
        """

        self.nx_graph = nx.Graph()

        long_simple_frags = []
        for frag in self.longer_fragments:
            long_simple_frags.extend(frag.simple_fragments)
        long_simple_frags = set(long_simple_frags)

        short_simple_frags = []
        for frag in self.shorter_fragments:
            short_simple_frags.extend(frag.simple_fragments)
        short_simple_frags = set(short_simple_frags)

        both_simple_frags = set(long_simple_frags.intersection(short_simple_frags))
        long_simple_frags = set(long_simple_frags.difference(both_simple_frags))
        short_simple_frags = set(short_simple_frags.difference(both_simple_frags))

        for u in list(both_simple_frags):
            self.nx_graph.add_node(str(u), frag="both")

        for u in list(long_simple_frags):
            self.nx_graph.add_node(str(u), frag="longer")

        for u in list(short_simple_frags):
            self.nx_graph.add_node(str(u), frag="shorter")

        all_nodes = list(set(both_simple_frags.union(long_simple_frags).union(short_simple_frags)))

        # need to make this based on fragments, not simple fragments
        for u in :
            for v in all_nodes:
                if u.intersect(v):
                    self.nx_graph.add_edge(str(u), str(v))

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
        nx.draw_networkx_nodes(self.nx_graph, pos, nodelist=purple_nodes, \
                node_size=node_size, node_color='#aa8cc5')
        nx.draw_networkx_nodes(self.nx_graph, pos, nodelist=blue_nodes, \
                node_size=node_size, node_color='b')
        nx.draw_networkx_nodes(self.nx_graph, pos, nodelist=red_nodes, \
                node_size=node_size, node_color='r')

        nx.draw_networkx_edges(self.nx_graph, pos)
        nx.draw_networkx_labels(self.nx_graph, pos, font_size=font_size)

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
        both_fragments = [Fragment(string=s) for s,d in self.nx_graph.nodes_iter(data=True) if d['frag'] == 'both']
        sources = list(filter(lambda x: x.is_before(self.deletion), both_fragments))
        targets = list(filter(lambda x: self.deletion.is_before(x), both_fragments))

        with open(fn_src_tgt, 'w') as f:
            f.write("#Node Node type\n")
            for src in sources:
                f.write(str(src) + "\tsource")
                f.write("\n")
            for tgt in targets:
                f.write(str(tgt) + "\ttarget")
                f.write("\n")

        # run PathLinker
        s = "python ../PathLinker/PathLinker.py " + fn_net + " " + fn_src_tgt
        os.system(s)

        # read from file
        fn = "out_k_{}-ranked-edges.txt".format(k)
        with open(fn, 'r') as f:
            for line in f:
                if line[0] == '#':
                    continue
                l = line.split()
                tail = Fragment(l[0])
                head = Fragment(l[1])
                index = l[2]
                print(tail, head, index)

    def draw_lines(self):
        pass

    def print_edgelist(self):
        assert(self.nx_graph is not None)
        for line in nx.generate_edgelist(self.nx_graph, data=False):
            print(line)

    def print_fragments(self):
        print("longer fragments:")
        pprint.pprint(self.longer_fragments)
        print("shorter fragments:")
        pprint.pprint(self.shorter_fragments)
        print("deletion: ", self.deletion)

    def __str__(self):
        return "longer: " + str(len(self.longer_fragments)) + \
                "\nshorter: " + str(len(self.shorter_fragments))

    def __repr__(self):
        return str(self)


G = StructuralSVGraphSimulator()
print(G)
print(G.nx_density())
G.draw_nx("graph.pdf")
#G.print_fragments()
##G.find_k_shortest_paths(100)

