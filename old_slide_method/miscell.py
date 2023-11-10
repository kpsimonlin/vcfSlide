## Provide some miscellaneous functions to main script or window statistics script.
## Kung-Ping Lin, 20230313.

import math
import numpy as np
from itertools import combinations
import random
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
    random.shuffle(mergedGroupIDs)
    return combinations(mergedGroupIDs, group1Len)

def permuP(scores, stats, maxPerm, maxExtm, minExtm=0.5, *args, **kwargs):
# Return the p-value of the given score using a permutation test
# scores: a list of given divergence scores
# stats: a list of functions calculating the divergence score
# maxPerm: maximum number of permutations
# maxExtm: output the p-value if this number of score that is more extreme than the given score is observed
# minExtm: minimum number of extreme values if no extreme values were observed in permutations, minExtm / maxPerm = the given p-value of this situation
    sTime = time.time()

    extremeNum = [0] * len(stats)   # nuber of extreme values observed in permutations
    returnP = [-1] * len(stats)     # the p-values that is going to be returned
    decompGroupIDs = kwargs['groupIDs'][0] + kwargs['groupIDs'][1]
    g0Len = len(kwargs['groupIDs'][0])
    g1Len = len(kwargs['groupIDs'][1])
    for i in range(maxPerm):
        print(i)
        print(extremeNum)
        print(returnP)
        print(time.time() - sTime)
        print()
        g0IDs = random.sample(decompGroupIDs, g0Len)
        g1IDs = list(set(decompGroupIDs) - set(g0IDs))
        newGroupIDs = [g0IDs, g1IDs]
        kwargs['groupIDs'] = newGroupIDs

        # Call the divergence score function
        for s, stat in enumerate(stats):
            if extremeNum[s] < maxExtm:
                permScore = stat(*args, **kwargs)

                if permScore >= scores[s]:
                # if the permutation score is as extreme as or more extreme than the observed score
                    extremeNum[s] += 1

                    if extremeNum[s] == maxExtm:
                    # If the number of extreme score is equal to the set threshold, return the p-value of the permutation test
                        returnP[s] = maxExtm / (i + 1)

                        if returnP.count(-1) == 0:
                        # If all test reach maximum extreme values
                            return returnP

    # If the number of permutation loops exceeds the given threshold, report the given p-value for those statistics that have 0 extreme values, report p for the rest of the statistics with < maxExtm extreme values
    for i, p in enumerate(returnP):
        if p == -1:
            en = extremeNum[i]
            if en == 0:
            # if no extreme values were observed, given p as there's minExtm extreme value observed
                returnP[i] = minExtm / maxPerm
            else:
                returnP[i] = en / maxPerm
    return returnP
