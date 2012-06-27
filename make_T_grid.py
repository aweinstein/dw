import networkx as nx
import numpy as np
from scipy.io import savemat

def make_T(N):
    '''Compute the diffusion operator of an N-by-N grid graph.

    The difussion operator is defined as

    T = D^0.5 * P * D^-0.5, where D is the degree matrix, P=D^-1W is the random
    walk matrix, and W is the adjancency matrix.
    '''
    G = nx.grid_2d_graph(N,N)
    W = nx.to_numpy_matrix(G, nodelist=sorted(G.nodes()))
    L = nx.laplacian(G, nodelist=sorted(G.nodes()))
    D = W.sum(1)
    P = W / D
    D = np.array(D).flatten()
    Dsq = np.sqrt(D)
    Disq = (1 / Dsq)
    T = np.dot(np.diag(Dsq), P)
    T = np.dot(T, np.diag(Disq))
    T = (T + T.T) / 2 # Iron out numerical wrinkles
    
    return T, L

for N in [4, 8, 16]:
    T, L = make_T(N)
    fn = 'T_grid_%d' % N
    savemat(fn, {'T':T, 'L':L}, oned_as='column')
    print 'Saved as %s' % fn
