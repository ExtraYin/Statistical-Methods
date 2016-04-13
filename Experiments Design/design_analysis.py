"""
Created on 2016/4/12
@author: Yida Yin
"""

import numpy as np
import itertools


def calc_var_diff_bet_pair(design_matrix):
    """
    return the variance/sigma^2 of difference of every two treatments given a known design matrix
    """
    # calculate X'X
    xx = np.dot(design_matrix.T, design_matrix)
    if np.linalg.matrix_rank(xx) < min(xx.shape):
        print "Matrix Not Invertible"
        return None
    xxi = np.linalg.inv(xx)
    # note that xxi is a square matrix
    n = xxi.shape[0]
    var = {}
    for pair in itertools.combinations(range(n), 2):  # TODO: the last treatment is not in the \beta vector
        # generate pairwise contrast
        contrast = [0 for i in range(n)]
        contrast[pair[0]], contrast[pair[1]] = 1, -1
        contrast = np.array(contrast)
        # calculate A'(X'X)^(-1)A
        v = reduce(np.dot, [contrast.T, xxi, contrast])
        v = round(v, 3)
        if v in var.keys():
            var[v].append(pair)
        else:
            var[v] = [pair]
    return var


def generate_design_matrix(design):
    """
    Find the design matrix for the given design
    The constrain condition set the last treatment=0, that is the \beta vector has the length #blocks + (#treatments-1)
    :param design:  A list of list, each sub-list stands for a block, each element in the sub-list is a treatment
    :return:  A design matrix
    """
    assert(isinstance(design[0], list))
    b = len(design)
    # count total treatment numbers and add their index to the dictionary
    treatments = set()
    for i, blk in enumerate(design):
        for trt in blk:
            treatments.add(trt)
    n = len(treatments)
    treatments = list(treatments)
    treatments.sort()
    treatments = {trt:i for i,trt in enumerate(treatments)}
    # create matrix
    design_matrix = []
    for bi, blk in enumerate(design):
        for trt in blk:
            row = [0 for i in range(b+n-1)]
            row[bi] = 1
            if treatments[trt] != n-1:
                row[b+treatments[trt]] = 1
            design_matrix.append(row)
    return design_matrix


def occurrence(design):
    """
    return the frequency of every treatment pairs which appear in the same block
    """
    pairs = []
    pair_occurance = []
    for blk in design:
        for pair in itertools.combinations(blk, 2):
            if set(pair) not in pairs:
                pairs.append(set(pair))
                pair_occurance.append(1)
            else:
                pair_occurance[pairs.index(set(pair))] += 1
    return {tuple(key): value for key, value in zip(pairs, pair_occurance)}



if __name__ == "__main__":
    a = np.array([[1,0,0,0,0,0,0,0,0,0,0,0,0,0,   0,1,0,0,0,0],
                 [1,0,0,0,0,0,0,0,0,0,0,0,0,0,   0,0,1,0,0,0],
                 [1,0,0,0,0,0,0,0,0,0,0,0,0,0,   0,0,0,0,0,1],
                 [1,0,0,0,0,0,0,0,0,0,0,0,0,0,   0,0,0,0,0,0],

                 [0,1,0,0,0,0,0,0,0,0,0,0,0,0,   0,0,0,1,0,0],
                 [0,1,0,0,0,0,0,0,0,0,0,0,0,0,   0,0,0,0,1,0],
                 [0,1,0,0,0,0,0,0,0,0,0,0,0,0,   0,0,0,0,0,1],
                 [0,1,0,0,0,0,0,0,0,0,0,0,0,0,   0,0,0,0,0,0],

                 [0,0,1,0,0,0,0,0,0,0,0,0,0,0,   0,0,0,1,0,0],
                 [0,0,1,0,0,0,0,0,0,0,0,0,0,0,   0,0,0,0,1,0],
                 [0,0,1,0,0,0,0,0,0,0,0,0,0,0,   0,0,0,0,0,1],
                 [0,0,1,0,0,0,0,0,0,0,0,0,0,0,   0,0,0,0,0,0],

                 [0,0,0,1,0,0,0,0,0,0,0,0,0,0,   1,0,0,0,0,0],
                 [0,0,0,1,0,0,0,0,0,0,0,0,0,0,   0,0,1,0,0,0],
                 [0,0,0,1,0,0,0,0,0,0,0,0,0,0,   0,0,0,1,0,0],
                 [0,0,0,1,0,0,0,0,0,0,0,0,0,0,   0,0,0,0,1,0],

                 [0,0,0,0,1,0,0,0,0,0,0,0,0,0,   1,0,0,0,0,0],
                 [0,0,0,0,1,0,0,0,0,0,0,0,0,0,   0,0,1,0,0,0],
                 [0,0,0,0,1,0,0,0,0,0,0,0,0,0,   0,0,0,1,0,0],
                 [0,0,0,0,1,0,0,0,0,0,0,0,0,0,   0,0,0,0,0,0],

                 [0,0,0,0,0,1,0,0,0,0,0,0,0,0,   1,0,0,0,0,0],
                 [0,0,0,0,0,1,0,0,0,0,0,0,0,0,   0,0,1,0,0,0],
                 [0,0,0,0,0,1,0,0,0,0,0,0,0,0,   0,0,0,0,0,1],
                 [0,0,0,0,0,1,0,0,0,0,0,0,0,0,   0,0,0,0,0,0],

                 [0,0,0,0,0,0,1,0,0,0,0,0,0,0,   1,0,0,0,0,0],
                 [0,0,0,0,0,0,1,0,0,0,0,0,0,0,   0,0,1,0,0,0],
                 [0,0,0,0,0,0,1,0,0,0,0,0,0,0,   0,0,0,0,1,0],
                 [0,0,0,0,0,0,1,0,0,0,0,0,0,0,   0,0,0,0,0,1],

                 [0,0,0,0,0,0,0,1,0,0,0,0,0,0,   0,1,0,0,0,0],
                 [0,0,0,0,0,0,0,1,0,0,0,0,0,0,   0,0,1,0,0,0],
                 [0,0,0,0,0,0,0,1,0,0,0,0,0,0,   0,0,0,0,1,0],
                 [0,0,0,0,0,0,0,1,0,0,0,0,0,0,   0,0,0,0,0,1],

                 [0,0,0,0,0,0,0,0,1,0,0,0,0,0,   0,1,0,0,0,0],
                 [0,0,0,0,0,0,0,0,1,0,0,0,0,0,   0,0,1,0,0,0],
                 [0,0,0,0,0,0,0,0,1,0,0,0,0,0,   0,0,0,1,0,0],
                 [0,0,0,0,0,0,0,0,1,0,0,0,0,0,   0,0,0,0,0,0],

                 [0,0,0,0,0,0,0,0,0,1,0,0,0,0,   1,0,0,0,0,0],
                 [0,0,0,0,0,0,0,0,0,1,0,0,0,0,   0,1,0,0,0,0],
                 [0,0,0,0,0,0,0,0,0,1,0,0,0,0,   0,0,0,0,1,0],
                 [0,0,0,0,0,0,0,0,0,1,0,0,0,0,   0,0,0,0,0,0],

                 [0,0,0,0,0,0,0,0,0,0,1,0,0,0,   1,0,0,0,0,0],
                 [0,0,0,0,0,0,0,0,0,0,1,0,0,0,   0,1,0,0,0,0],
                 [0,0,0,0,0,0,0,0,0,0,1,0,0,0,   0,0,0,0,1,0],
                 [0,0,0,0,0,0,0,0,0,0,1,0,0,0,   0,0,0,0,0,0],

                 [0,0,0,0,0,0,0,0,0,0,0,1,0,0,   1,0,0,0,0,0],
                 [0,0,0,0,0,0,0,0,0,0,0,1,0,0,   0,1,0,0,0,0],
                 [0,0,0,0,0,0,0,0,0,0,0,1,0,0,   0,0,0,1,0,0],
                 [0,0,0,0,0,0,0,0,0,0,0,1,0,0,   0,0,0,0,0,1],

                 [0,0,0,0,0,0,0,0,0,0,0,0,1,0,   1,0,0,0,0,0],
                 [0,0,0,0,0,0,0,0,0,0,0,0,1,0,   0,1,0,0,0,0],
                 [0,0,0,0,0,0,0,0,0,0,0,0,1,0,   0,0,0,1,0,0],
                 [0,0,0,0,0,0,0,0,0,0,0,0,1,0,   0,0,0,0,0,1],

                 [0,0,0,0,0,0,0,0,0,0,0,0,0,1,   0,1,0,0,0,0],
                 [0,0,0,0,0,0,0,0,0,0,0,0,0,1,   0,0,1,0,0,0],
                 [0,0,0,0,0,0,0,0,0,0,0,0,0,1,   0,0,0,1,0,0],
                 [0,0,0,0,0,0,0,0,0,0,0,0,0,1,   0,0,0,0,1,0]])

    vd = calc_var_diff_bet_pair(a)
    for key in vd:
        print key, vd[key]


    design = [['A', 'C', 'D', 'E'],
              ['A', 'C', 'D', 'E'],
              ['D', 'E', 'F', 'G'],
              ['D', 'E', 'F', 'G'],
              ['A', 'C', 'F', 'G'],
              ['A', 'C', 'F', 'G'],
              ['B', 'C', 'E', 'F'],
              ['B', 'C', 'E', 'G'],
              ['A', 'B', 'E', 'G'],
              ['B', 'C', 'D', 'G'],
              ['A', 'B', 'E', 'F'],
              ['A', 'B', 'D', 'G'],
              ['A', 'B', 'D', 'F'],
              ['B', 'C', 'D', 'F']]

    b = generate_design_matrix(design)
    b = np.matrix(b)
    print np.linalg.det(np.dot(a.T,a))
    print np.linalg.det(np.dot(b.T,b))

    print occurrence(design)

    design2 = [['A', 'C', 'D', 'E'],
              ['D', 'E', 'F', 'F']]
    print occurrence(design2)