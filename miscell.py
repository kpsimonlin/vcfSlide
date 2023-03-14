## Provide some miscellaneous functions to main script or window statistics script.
## Kung-Ping Lin, 20230313.

import math
import numpy as np
from itertools import combinations
from random import shuffle
import time

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

def randomComb(groupIDs):
# Take a list with two sublists, return a iterator of randomized new combinations, used for permuation test
# groupIDs = [[0, 1, 2, 3], [4, 5, 6, 7]] for two group IDs.
    group1Len = len(groupIDs[0])
    mergedGroupIDs = groupIDs[0] + groupIDs[1]
    shuffle(mergedGroupIDs)
    return combinations(mergedGroupIDs, group1Len)

def permuP(scores, stats, combIter, maxPerm, maxExtm, maxP=5e-7, *args, **kwargs):
# Return the p-value of the given score using a permutation test
# scores: a list of given divergence scores
# stats: a list of functions calculating the divergence score
# combIter: an iterator running randomized combinations of groupIDs for permutation test
# maxPerm: maximum number of permutations
# maxExtm: output the p-value if this number of score that is more extreme than the given score is observed
# maxP: the p-value that will be return if maxPerm is reached
    permNum = 1
    extremeNum = [0] * len(stats)
    returnP = [-1] * len(stats)    # the p-values that is going to be returned
    decompGroupIDs = kwargs['groupIDs'][0] + kwargs['groupIDs'][1]
    for g1IDs in combIter:
        if permNum > maxPerm:
        # If the number of permutation loops exceeds the given threshold, abort and report maxP for those statistics that haven't received maxExtm number of extreme values
            for i, p in enumerate(returnP):
                if p == -1:
                    returnP[i] = maxP
            return returnP

        g1IDs = list(g1IDs)
        g2IDs = [x for i, x in enumerate(decompGroupIDs) if x not in g1IDs]
        newGroupIDs = [g1IDs, g2IDs]
        kwargs['groupIDs'] = newGroupIDs

        sTime = time.time()

        # Call the divergence score function
        for i, stat in enumerate(stats):
            if extremeNum[i] < maxExtm:
                permScore = stat(*args, **kwargs)

                if permScore >= scores[i]:
                # if the permutation score is as extreme as or more extreme than the observed score
                    extremeNum[i] += 1
                    if extremeNum[i] == maxExtm:
                    # If the number of extreme score is equal to the set threshold, return the p-value of the permutation test
                        returnP[i] = maxExtm / permNum
            print(str(i) + ' ' + str(time.time() - sTime))
            sTime = time.time()

        permNum += 1

        print(permNum)
        print(returnP)
