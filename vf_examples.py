'''Create graphs for environments described in 

Maggioni, M., and Mahadevan, S. (2006). A Multiscale Framework for Markov
Decision Processes using Diffusion Wavelets ( No. 2006-36) (pp. 1-44).
'''
import networkx as nx
import numpy as np
import matplotlib.pyplot as plt
from scipy.io import savemat

from nxutils import plot_grid
from diffusion import get_T

def ex_1():
    '''Example 1 of the paper.

    A 2-by-20 grid with reflecting boundary condition and an obstacle in the
    middle.
    '''
    N = 20
    G = nx.grid_2d_graph(N,N)

    n_middle = (N/2, N/2)

    # Add reflecting boundary condition
    for i,j in G:
        if i == 0 or i == N-1 or j == 0 or j == N-1:
            G.add_edge((i,j), (i,j))
    
    # Remove edges of the node at the middle
    # (keep one edge, to simplify the bookiping
    for u,v in sorted(G.edges(n_middle))[:-1]:
        G.add_edge(v,v) # Add self loop
        G.remove_edge(u, v)

    T, _ = get_T(G)
    savemat('vf_ex1', {'T':T}, oned_as='column')
    return G

if __name__ == '__main__':
    G = ex_1()
    #plot_grid(G)
    #plt.show()
    
