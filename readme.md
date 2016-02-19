# Genome Graph Simulation

Simulates some data so that we can test out some algorithmms - vague, I know. 


## TODO
- fix putting fragments in networkx. 
    - fragments that start at the same point will be considered the same node, even when they shouldn't be. 
    - One soln is to name the node : `str(list.index(node)) + str(node)`, so that each node is guaranteed to have a unique name.
