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
    for pair in itertools.combinations(range(n), 2):  # ALREADY FIXED: the last treatment is not in the \beta vector
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
    # compare with the least treatment
    for trt_1 in range(n):
        contrast = [0 for i in range(n)]
        contrast[trt_1] = 1
        contrast = np.array(contrast)
        # calculate A'(X'X)^(-1)A
        v = reduce(np.dot, [contrast.T, xxi, contrast])
        v = round(v, 3)
        if v in var.keys():
            var[v].append(tuple([trt_1, n]))
        else:
            var[v] = [tuple([trt_1, n])]
    return var


def calc_trt_eff(design_matrix):
    """
    calculate treatment effects
    return (X'X)^{-1}X'Y     but we don't know Y
    """
    pass


def control_first(t):
    """
    usually 'O' stands for 'control group'
    """
    if t == 'O':
        return '0'
    return t


def generate_design_matrix(design):
    """
    Find the design matrix for the given design
    The constrain condition set the last treatment=0, that is the \beta vector has the length #blocks + (#treatments-1)
    :param design:  A list of list, each sub-list stands for a block, each element in the sub-list is a treatment
    :return:  A design matrix
    """
    assert (isinstance(design[0], list))
    b = len(design)
    # count total treatment numbers and add their index to the dictionary
    treatments = set()
    for i, blk in enumerate(design):
        for trt in blk:
            treatments.add(trt)
    n = len(treatments)
    treatments = list(treatments)
    treatments.sort(key=control_first)
    treatments = {trt: i for i, trt in enumerate(treatments)}
    # create matrix
    design_matrix = []
    for bi, blk in enumerate(design):
        for trt in blk:
            row = [0 for i in range(b + n - 1)]
            row[bi] = 1
            if treatments[trt] != n - 1:
                row[b + treatments[trt]] = 1
            design_matrix.append(row)
    return np.array(design_matrix)


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
    design75 = [['O', 'A', 'B', 'C'],
                ['O', 'A', 'B'],
                ['O', 'A', 'C'],
                ['O', 'B', 'C']]
    print occurrence(design75)
    print calc_var_diff_bet_pair(generate_design_matrix(design75))
    print '-' * 80

    design77a = [['B', 'C', 'E', 'F'],
                 ['A', 'D', 'E', 'F'],
                 ['A', 'C', 'D', 'E'],
                 ['A', 'C', 'D', 'F'],
                 ['B', 'C', 'D', 'F'],
                 ['B', 'C', 'D', 'E'],
                 ['A', 'B', 'E', 'F'],
                 ['A', 'B', 'D', 'E'],
                 ['A', 'B', 'C', 'F']]
    print occurrence(design77a)
    print calc_var_diff_bet_pair(generate_design_matrix(design77a))
    print '-' * 80

    design77b = [['A', 'B', 'D', 'F'],
                 ['A', 'B', 'E', 'F'],
                 ['A', 'D', 'E', 'F'],
                 ['A', 'B', 'D', 'E'],
                 ['A', 'B', 'C', 'D'],
                 ['A', 'C', 'D', 'E'],
                 ['A', 'C', 'E', 'F'],
                 ['A', 'B', 'C', 'E'],
                 ['A', 'C', 'D', 'F'],
                 ['A', 'B', 'C', 'E']]
    print occurrence(design77b)
    print calc_var_diff_bet_pair(generate_design_matrix(design77b))
    print '-' * 80

    design78 = [['A', 'B', 'C'],
                ['A', 'B', 'C'],
                ['A', 'D', 'E'],
                ['B', 'D', 'E'],
                ['C', 'D', 'E'], ]
    print occurrence(design78)
    print calc_var_diff_bet_pair(generate_design_matrix(design78))



