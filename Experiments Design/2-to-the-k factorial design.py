"""
@author: Yida Yin

Copyright (c) 16/5/3 Yida Yin. All rights reserved.
"""
import string
import itertools


def standard_order(k):
    # return the standard order of treatment A,B,......
    # which is (1), a, b, ab, c, ac, bc, abc, ......
    if k <= 0:
        return []
    order = []
    for i in range(2 ** k):
        s = ""
        for j in range(k + 1):
            if i % (2 ** (j + 1)) >= 2 ** j:
                s += string.lowercase[j]
        order.append(s)
    order[0] = '(1)'
    return order


def show_block(k, contrasts):
    # show factor-level combinations that go into each block
    block = {}
    for trt in standard_order(k):
        oe = ''  # odd-even string
        for ctt in contrasts:
            oe += str(len([letter for letter in trt if letter in ctt.lower()]) % 2)  # odd-even method
        if oe in block.keys():
            block[oe].append(trt)
        else:
            block[oe] = [trt]
    return block


def show_confounded(contrasts):
    # show al the confounded effects
    confounded = []
    for i in range(1, len(contrasts)+1):
        for combination in set(itertools.combinations(contrasts, i)):
            r = set()
            for single in combination:
                for letter in single:
                    if letter in r:
                        r -= set(letter)
                    else:
                        r.add(letter)
                confounded.append(''.join(r))
    return set(confounded)


def show_plus_minus(k, write=False):    # TODO
    print "     ",
    for i in range(1, k+1):
        for combination in list(itertools.combinations(string.uppercase[:k], i)):
            print ''.join(combination),
    print
    sd_order = standard_order(k)
    for trt in sd_order:
        print trt,

        print


if __name__ == "__main__":
    print standard_order(0) == []
    print standard_order(3) == ['(1)', 'a', 'b', 'ab', 'c', 'ac', 'bc', 'abc']
    print standard_order(4) == ['(1)', 'a', 'b', 'ab', 'c', 'ac', 'bc', 'abc', 'd', 'ad', 'bd', 'abd', 'cd', 'acd',
                                'bcd', 'abcd']
    print show_block(4, ['AD', 'BCD']) == {'11': ['ab', 'ac', 'd', 'bcd'], '10': ['a', 'abc', 'bd', 'cd'],
                                           '00': ['(1)', 'bc', 'abd', 'acd'], '01': ['b', 'c', 'ad', 'abcd']}
    print show_confounded(['ABCF', 'ABDE', 'ACDE', 'BCDH']) == {'BEFH', 'CEDF', 'ACEH', 'ACBF', 'BEDF', 'AF', 'CB',
                                                                'CEFH', 'ACBDFH', 'ACED', 'HD', 'ABED', 'ADFH', 'HCBD',
                                                                'ABEH'}
    show_plus_minus(3)
