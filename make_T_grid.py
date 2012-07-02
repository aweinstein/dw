import networkx as nx
import numpy as np
from scipy.io import savemat

from diffusion import get_T

def make_T(N, add_selfloops=False):
    '''Compute the diffusion operator of an N-by-N grid graph.

    The difussion operator is defined as

    T = D^0.5 * P * D^-0.5, where D is the degree matrix, P=D^-1W is the random
    walk matrix, and W is the adjancency matrix.
    '''
    G = nx.grid_2d_graph(N,N)
    if add_selfloops:
        for n in G:
            G.add_edge(n, n)
            
    T, _ = get_T(G)
    
    return T

for N in [4, 8, 16]:
    T = make_T(N)
    fn = 'T_grid_%d' % N
    savemat(fn, {'T':T}, oned_as='column')
    print 'Saved as %s' % fn

    T = make_T(N, True)
    fn = 'T_grid_sl_%d' % N
    savemat(fn, {'T':T}, oned_as='column')
    print 'Saved as %s' % fn

