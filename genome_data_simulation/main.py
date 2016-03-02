#import collections
# GenomeGraphSimulation (master) conda create -n compbio python=3 anaconda
import copy
import matplotlib.pyplot as plt
import numpy as np
import networkx as nx
import pprint
import os
import random

class SimpleFragment:
    """ a simple simple fragment with a start and stop"""
    def __init__(self, start=0, end=0, string=None):
        if string is None:
            self.start = start
            self.end = end
            self.key = random.randint(0,10000)
        else:
            s = string.strip().split("_")
            self.start = int(s[0])
            self.end = int(s[1])
            self.key = int(s[2])

    def add(self, other):
        self.end += len(other)

    def intersect(self, other):
        """ Returns true if self and other intersect
        :param other: SimpleFragment
        """
        if self.start <= other.start <= self.end or \
                self.start <= other.end <= self.end or \
                other.start <= self.start <= other.end or \
                other.start <= self.end <= other.end:
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

    def frag_thats_before(self, other):
        """ returns portion of frag thats before other"""
        assert(len(self.simple_fragments) == len(other.simple_fragments) == 1\
                and "only works for simple case")
        if self.simple_fragments[0].start >= other.simple_fragments[0].start:
            return SimpleFragment(0,0)
        else:
            return SimpleFragment(self.simple_fragments[0].start , other.simple_fragments[0].start)

    def frag_thats_intersecting(self, other):
        """ returns portion of frag thats before other"""
        assert(len(self.simple_fragments) == len(other.simple_fragments) == 1\
                and "only works for simple case")
        e = None
        s = None
        if self.simple_fragments[0].intersect(other.simple_fragments[0]):
            if self.simple_fragments[0].start <= other.simple_fragments[0].start:
                s = other.simple_fragments[0].start
            else:
                s = self.simple_fragments[0].start

            if self.simple_fragments[0].end <= other.simple_fragments[0].end:
                e = self.simple_fragments[0].end
            else:
                e = other.simple_fragments[0].end
        else:
            s = 0
            e = 0
        return SimpleFragment(s,e)

    def frag_thats_after(self, other):
        """ returns portion of frag thats after other"""
        assert(len(self.simple_fragments) == len(other.simple_fragments) == 1\
                and "only works for simple case")
        if other.simple_fragments[0].end >= self.simple_fragments[0].end:
            return SimpleFragment(other.simple_fragments[0].end, other.simple_fragments[0].end)
        else:
            return SimpleFragment(other.simple_fragments[0].end, self.simple_fragments[0].end)



    def get_overlap(self, other):
        assert(len(self.simple_fragments) == len(other.simple_fragments) == 1\
                and "only works for simple case")
        if self.intersect(other) == False:
            return SimpleFragment(0,0)
        elif self.is_contained_in(other):
            return SimpleFragment(self.simple_fragment[0].start, self.simple_fragment[0].end)
        elif other.is_contained_in(self):
            return SimpleFragment(other.simple_fragment[0].start, other.simple_fragment[0].end)
        elif self.start < other.start:
            return SimpleFragment(0, self.end - other.start)
        else:
            return SimpleFragment(0, other.end - self.start)





class StructuralSVGraphSimulator:
    def __init__(self):
        simple = True
        if simple:
            print("using simple parameters")
            self.longer_seq = 50
            self.deletion = Fragment([SimpleFragment(20, 30)])
            self.coverage = 1 # coverage per shorter/longer
            self.poisson_lambda = 4
            self.bases_per_fragment = 6
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
                    self.longer_fragments.append(Fragment([SimpleFragment(start, start + self.bases_per_fragment)]))
                    prev = start
                    count += 1
                else:
                    break

    def create_shorter_fragments(self):
        """creates shorter_fragment based on longer_fragments.
        assumes longer_fragments is already populated"""
        # again, writing this in the most naive, unpythonic way
        for _ in range(self.coverage):
            print('break')
            prev = 0
            while True:
                start = prev + np.random.poisson(self.poisson_lambda)
                end = start + self.bases_per_fragment
                frag = Fragment([SimpleFragment(start, end)])

                if end >= self.longer_seq: # and end of sequence
                    print('here')
                    break

                if self.deletion.intersect(frag):
                    # find how much overlap
                    before = frag.frag_thats_before(self.deletion)
                    intersecting = frag.frag_thats_intersecting(self.deletion)
                    after = frag.frag_thats_after(self.deletion)
                    after.add(intersecting)

                    if len(before) == 0 and len(after) == 0:
                        print('case 1')
                        prev = self.deletion.simple_fragments[0].end
                        pass
                    elif len(before) == 0 and len(after) != 0:
                        print('case 2')
                        self.shorter_fragments.append(Fragment([after]))
                        prev = after.end

                    elif len(before) != 0 and len(after) == 0:
                        print('case 3')
                        self.shorter_fragments.append(Fragment([before]))
                        prev = self.deletion.simple_fragments[0].end

                    elif len(before) != 0 and len(after) != 0:
                        print('case 4')
                        print(frag)
                        self.shorter_fragments.append(Fragment([before, after]))
                        prev = after.end
                    else:
                        print('bug')
                    print(self.shorter_fragments[-1].simple_fragments)

                else:
                    self.shorter_fragments.append(frag)
                    prev = start
                print(prev)



    def make_networkx_graph(self):
        """ based on self.longer_fragments and self.shorter_fragments,
        make a networkx graph"""

        """
        our nodes are really Fragments (not SimpleFragments)
        labels should be: before, after, long, short (as opposed to current labels)
        """

        self.nx_graph = nx.Graph()


        for frag in self.longer_fragments:
            if frag.is_before(self.deletion):
                self.nx_graph.add_node(str(frag), frag="before")
            elif self.deletion.is_before(frag):
                self.nx_graph.add_node(str(frag), frag="after")
            else:
                self.nx_graph.add_node(str(frag), frag="long")

        for frag in self.shorter_fragments:
            if frag.is_before(self.deletion):
                self.nx_graph.add_node(str(frag), frag="before")
            elif self.deletion.is_before(frag):
                self.nx_graph.add_node(str(frag), frag="after")
            else:
                self.nx_graph.add_node(str(frag), frag="short")

        # need to make this based on fragments, not simple fragments


        longer = copy.deepcopy(self.longer_fragments)
        shorter = copy.deepcopy(self.shorter_fragments)
        all_nodes = list(set(longer).union(set(shorter)))

        for u in all_nodes:
            for v in all_nodes:
                if u != v and u.intersect(v):
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

        blue_nodes = [n for n,d in self.nx_graph.nodes_iter(data=True) if d['frag'] == 'long']
        red_nodes = [n for n,d in self.nx_graph.nodes_iter(data=True) if d['frag'] == 'short']
        purple_nodes = [n for n,d in self.nx_graph.nodes_iter(data=True) if d['frag'] == 'before']
        green_nodes = [n for n,d in self.nx_graph.nodes_iter(data=True) if d['frag'] == 'after']
        nx.draw_networkx_nodes(self.nx_graph, pos, nodelist=purple_nodes, \
                node_size=node_size, node_color='#aa8cc5')
        nx.draw_networkx_nodes(self.nx_graph, pos, nodelist=blue_nodes, \
                node_size=node_size, node_color='b')
        nx.draw_networkx_nodes(self.nx_graph, pos, nodelist=red_nodes, \
                node_size=node_size, node_color='r')
        nx.draw_networkx_nodes(self.nx_graph, pos, nodelist=green_nodes, \
                node_size=node_size, node_color='g')

        nx.draw_networkx_edges(self.nx_graph, pos)
        #nx.draw_networkx_labels(self.nx_graph, pos, font_size=font_size)

        plt.title("red=shorter, blue=longer, purple=before, green=after")
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

    def draw_degree_histogram(self):
        assert(self.nx_graph is not None and "should be made by now")
        degree_sequence=sorted(nx.degree(self.nx_graph).values(),reverse=True) # degree sequence
        dmax=max(degree_sequence)
        plt.clf()
        plt.loglog(degree_sequence,'b-',marker='o')
        plt.title("Degree rank plot")
        plt.ylabel("degree")
        plt.xlabel("rank")
        plt.savefig("degree_histogram.pdf")

    def draw_lines(self):
        """ draws the cool (or lame, depending on your perspective) plot of all reads"""
        plt.clf()

        # draw longer fragments
        longer = copy.deepcopy(self.longer_fragments)
        shorter = copy.deepcopy(self.shorter_fragments)
        tot_frags = len(list(set(longer).union(set(shorter))))
        for i, frag in enumerate(self.longer_fragments):
            for f in frag.simple_fragments:
                xs = [f.start, f.end]
                ys = [(i / tot_frags) * 400] * 2
                plt.plot(xs, ys, 'b-', lw=1)

        # draw shorter fragments
        for i, frag in enumerate(self.shorter_fragments):
            j = len(self.longer_fragments) + i
            print(frag.simple_fragments)
            for f in frag.simple_fragments:
                xs = [f.start, f.end]
                ys = [(j / tot_frags) * 400] * 2
                plt.plot(xs, ys, 'r-', lw=1)

        # draw deletion
        xs = [self.deletion.simple_fragments[0].start, self.deletion.simple_fragments[0].start]
        ys = [0, 400]
        plt.plot(xs, ys, 'k-', lw=1)
        xs = [self.deletion.simple_fragments[0].end, self.deletion.simple_fragments[0].end]
        ys = [0, 400]
        plt.plot(xs, ys, 'k-', lw=1)
        plt.title("blue = longer, red = shorter")

        #plt.plot(xs, ys, 'r-', lw=1)

        plt.savefig("line.pdf")

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
#G.draw_nx("graph.pdf")
#G.draw_degree_histogram()
G.draw_lines()
#G.print_fragments()
##G.find_k_shortest_paths(100)

