"""
Created on #DATA
@author: Yida Yin
"""

import numpy as np


def calc_var_diff_bet_pair(design_matrix):
    xx = np.dot(design_matrix.T, design_matrix)
    if np.linalg.matrix_rank(xx) < min(xx.shape):
        print "Matrix Not Invertible"
        return None
    xxi = np.linalg.inv(xx)
    # note that xxi is a square matrix
    n = xxi.shape[0]
    return n

a = np.array([[1,0,0],[0,5,0], [0,0,1]])
print calc_var_diff_bet_pair(a)
