## Provide some miscellaneous functions to main script or window statistics script.
## Kung-Ping Lin, 20221117.

import math
import numpy as np
from itertools import combinations
import random

def calc_row_idx(k, n):
    return int(math.ceil((1/2.) * (- (-8*k + 4 *n**2 -4*n - 7)**0.5 + 2*n -1) - 1))

def elem_in_i_rows(i, n):
    return i * (n - 1 - i) + (i*(i + 1))//2

def calc_col_idx(k, i, n):
    return int(n - elem_in_i_rows(i + 1, n) + k)

def condensed_to_square(k, n):
    i = calc_row_idx(k, n)
    j = calc_col_idx(k, i, n)
    return i, j

def randomComb(groupIDs, seed=42689):
# Take a list with two sublists, return a iterator of randomized new combinations, used for permuation test
# groupIDs = [[0, 1, 2, 3], [4, 5, 6, 7]] for two group IDs.
    group1Len = len(groupIDs[0])
    mergedGroupIDs = groupIDs[0] + groupIDs[1]
    random.Random(seed).shuffle(mergedGroupIDs)
    return combinations(mergedGroupIDs, group1Len)

def permuP(score, stats, combIter, maxPerm, maxExtm, maxP=5e-7, *args, **kwargs):
# Return the p-value of the given score using a permutation test
# score: the given score
# stats: the function calculating the divergence score
# combIter: an iterator running randomized combinations of groupIDs for permutation test
# maxPerm: maximum number of permutations
# maxExtm: output the p-value if this number of score that is more extreme than the given score is observed
# maxP: the p-value that will be return if maxPerm is reached
    permNum = 1
    extremeNum = 0
    decompGroupIDs = kwargs['groupIDs'][0] + kwargs['groupIDs'][1]
    for g1IDs in combIter:
        if permNum > maxPerm:
        # If the number of permutation loops exceeds the given threshold, abort and report maxP
            return maxP

        if extremeNum > maxExtm:
        # If the number of extreme score is higher than the set thrshold, return the p-value of the permutation test
            return maxExtm / permNum

        g1IDs = list(g1IDs)
        g2IDs = [x for i, x in enumerate(decompGroupIDs) if x not in g1IDs]
        newGroupIDs = [g1IDs, g2IDs]
        kwargs['groupIDs'] = newGroupIDs

        # Call the divergence score function
        permScore = stats(*args, **kwargs)

        if permScore >= score:
        # if the permutation score is as extreme as or more extreme than the observed score
            extremeNum += 1

        permNum += 1
