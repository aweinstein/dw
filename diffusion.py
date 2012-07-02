import networkx as nx
import numpy as np

def get_T(G):
    ''' Return diffusion operator of a graph.

    The diffusion operator is defined as T = I - L, where L is the normalized
    Laplacian.

    Parameters
    ----------
    G : NetworkX graph
    
    Returns
    -------
    T : NumPy array
    Ln : NumPy array
      Normalized Laplacian of G.

    Notes
    -----
    Computing the normalized laplacian by hand. It seems there are some
    inconsistencies using nx.normalized_laplacian when G has selfloops.
    '''
    A = nx.adj_matrix(G, nodelist=sorted(G.nodes()))
    D = np.array(np.sum(A,1)).flatten()

    Disqrt = np.array(1 / np.sqrt(D))
    Disqrt = np.diag(Disqrt)
    L = np.diag(D) - A
    Ln = np.dot(np.dot(Disqrt, L), Disqrt)
    T =  np.eye(len(G)) - Ln
    T = (T + T.T) / 2 # Iron out numerical wrinkles
    
    return T, L

## def get_T_(G):
##     ''' Return diffusion operator of a graph.

##     The diffusion operator is defined as T = I - L, where L is the normalized
##     Laplacian.

##     Parameters
##     ----------
##     G : NetworkX graph
    
##     Returns
##     -------
##     T : NumPy array
##     L : NumPy array
##        Normalized Laplacian of G.
##     '''
##     L = nx.normalized_laplacian(G, nodelist=sorted(G.nodes()))
##     T =  np.eye(len(G)) - L
##     T = (T + T.T) / 2 # Iron out numerical wrinkles
    
##     return T, L
