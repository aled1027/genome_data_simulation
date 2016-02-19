#import collections
# GenomeGraphSimulation (master) conda create -n compbio python=3 anaconda
import copy
import numpy as np

class Fragment:
    def __init__(self, start=0, end=0):
        self.start = start
        self.end = end

    def intersect(self, other)

        """ Other is a fragment object.
        returns true if the self and other fragment object intersect each other"""

    def __len__(self):
        return self.end - self.start

    def __str__(self):
        return "Fragment(start={}, end={})".format(self.start, self.end)

    def __repr__(self):
        return str(self)

class GraphSimulator:

    # also using a poisson distribution

    def __init__(self):
        self.longer_seq = 300
        self.deletion = Fragment(50,150)
        self.coverage = 1 # TODO change to 30
        self.poisson_lambda = 50
        self.bases_per_fragment = 50
        self.longer_fragments = []
        self.shorter_fragments = []


    def run(self):
        """ main argument"""
        self.generate_fragments_for_longer_sequence()
        self.shorter_fragments = copy.deepcopy(self.longer_fragments)
        self.shorter_fragments = self.filter_around_deletion(self.shorter_fragments)

    def generate_graph(self):
        """ Using self.longer_fragments and self.shorter_fragments,
        generate a graph. """




    def filter_around_deletion(self, frags):
        def _avoids_deletion(x):
            """ takes a Fragment object as input.
            Returns True if Fragment avoids the deletion, else False"""
            if self.deletion.start <= x.start <= self.deletion.end or \
                    self.deletion.start <= x.end <= self.deletion.end:
                return False
            else:
                return True
        return filter(_avoids_deletion, frags)

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

G = GraphSimulator()
G.run()
print(G)

