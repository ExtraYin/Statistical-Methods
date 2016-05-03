"""
@author: Yida Yin

Copyright (c) 16/5/3 Yida Yin. All rights reserved.
"""
import string


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
        b = ''   # odd_even method
        for ctt in contrasts:
            b += str(len([letter for letter in trt if letter in ctt]) % 2)   # TODO
        if b in block.keys():
            block[b] += trt
        else:
            block[b] = [trt]
    return block


def show_confounded(k, contrasts):
    # show al the confounded effects
    pass

if __name__ == "__main__":
    print standard_order(0) == []
    print standard_order(3) == ['(1)', 'a', 'b', 'ab', 'c', 'ac', 'bc', 'abc']
    print standard_order(4) == ['(1)', 'a', 'b', 'ab', 'c', 'ac', 'bc', 'abc', 'd', 'ad', 'bd', 'abd', 'cd', 'acd',
                                'bcd', 'abcd']
    print show_block(5, ['bcd', 'ace'])
