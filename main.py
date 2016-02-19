#import collections
# GenomeGraphSimulation (master) ➜ conda create -n compbio python=3 anaconda
import random

class Fragment:
    def __init__(self, start=0, end=0):
        self.start = start
        self.end = end

    def __len__(self):
        return self.end - self.start

    def __str__(self):
        return "Fragment(start={}, end={})".format(self.start, self.end)

class GraphSimulator:
    len_seq_1 = 1000
    len_seq_2 = 1200
    deletion = Fragment(400,600)
    coverage = 1 # TODO change to 30
    bases_per_fragment = 50
    fragments = []
    # also using a poisson distribution

    def __init__(self):
        pass

    def run(self):
        prev = 0
        #for _ in range(self.coverage):
           #start = prev + random.randint(2,10) # TODO sample from poisson
           #self.fragments.append(Fragment(start, start + self.bases_per_fragment))
        self.fragments.append(Fragment(1,2))

    def __str__(self):
        return str(self.fragments)
        pass

G = GraphSimulator()
G.run()
print(G)


