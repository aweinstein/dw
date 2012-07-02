'''
Generate a circulant matrix with first column [0.25 0.5 0 ... 0.25].
'''

import numpy as np
import networkx as nx
import matplotlib.pyplot as plt
from scipy.linalg import circulant
from scipy.io import savemat

from diffusion import get_T

def get_matrices(G):
    '''

    Parameters
    ----------
    G : NetworkX graph

    Returns
    -------
    T : Numpy array

    L : NumPy array
        Combinatorial Laplacian of G
    W :

    '''
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

    return T, L, np.array(W)

def circle(N):
    G = nx.cycle_graph(N)
    T, L = get_T(G)

    fn = 'T_circle_%d.mat' % N
    savemat(fn, {'T':T, 'L':L}, oned_as='column')
    print 'T saved as %s' % fn

    # Add self loops to all nodes
    for node in G.nodes():
        G.add_edge(node, node)

    T_loop, L_loop = get_T(G)

    fn = 'T_circle_selfloop_%d.mat' % N
    savemat(fn, {'T':T_loop, 'L':L_loop}, oned_as='column')
    print 'T saved as %s' % fn

    return T, L, T_loop, L_loop

if __name__ == '__main__':
    N = 256
    T, L, T_loop, L_loop = circle(N)
