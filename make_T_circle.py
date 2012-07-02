'''
Generate a circulant matrix with first column [0.25 0.5 0 ... 0.25].
'''

import numpy as np
import networkx as nx
import matplotlib.pyplot as plt
from scipy.linalg import circulant
from scipy.io import savemat

def get_matrices(G):
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

    return T, L, W

def circle(N):
    G = nx.cycle_graph(N)
    T, L, W = get_matrices(G)

    fn = 'T_circle_%d.mat' % N
    savemat(fn, {'T':T, 'L':L, 'W':W}, oned_as='column')
    print 'T saved as %s' % fn

    # Add self loops to all nodes
    for node in G.nodes():
        G.add_edge(node, node)

    T, L, W = get_matrices(G)

    fn = 'T_circle_selfloop_%d.mat' % N
    savemat(fn, {'T':T, 'L':L, 'W':W}, oned_as='column')
    print 'T saved as %s' % fn


N = 256
circle(N)
## r = np.zeros(N)
## r[[0, 1, -1]] = [0.5, 0.25, 0.25]
## Tp = circulant(r).T
## savemat('Tp.mat', {'Tp':Tp}, oned_as='column')
## print 'Tp saved as Tp.mat'

## plot = False
## if plot:
##     plt.matshow(Tp)
##     plt.show()
