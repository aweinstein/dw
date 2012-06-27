'''
Generate a circulant matrix with first column [0.25 0.5 0 ... 0.25].
'''

import numpy as np
import matplotlib.pyplot as plt
from scipy.linalg import circulant
from scipy.io import savemat

N = 512
r = np.zeros(N)
r[[0, 1, -1]] = [0.5, 0.25, 0.25]
Tp = circulant(r).T
savemat('Tp.mat', {'Tp':Tp}, oned_as='column')
print 'Tp saved as Tp.mat'

plot = False
if plot:
    plt.matshow(Tp)
    plt.show()
