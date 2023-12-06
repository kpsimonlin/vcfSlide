## Module for methods combining statistics from several non-overlapping windows
## Kung-Ping Lin, 20231115.
## These methods take a list of "results" for input, the "results" are a two-element nested list structured for a single sliding window like this:
#### [[snp_pos1, snp_pos2 ... ], [snp_number, stats1, stats2, ... ]] where stats are reported statistics like CSS

import numpy as np
import window_stats as ws
import random
import allel as al

## Main function, of combining results from several sliding windows
def mergeResults(resultList, resultHeader):
# resultHeader is a list recording the statistics specified, starting with snp_number
    mergedResults = [[], []]                # the merged results list: [[snp_poss], [snp_number, stats1, ... ]]
    statsNum = len(resultList[0][1]) - 1    # the number of statistics ready to be combined w/o snp_number
    statsHeader = resultHeader[1:]          # statistics to combine
    hiddens = {}                            # dictionary for hiddens, if applicable will be filled

    ## Merge the number of SNPs and the statistics
    mergedNum = np.sum(r[1][0] for r in resultList)     # total number of SNPs in the combined window
    if mergedNum == 0:
    # if there is no SNPs in the combined window, return directly
        mergedResults[1].append(mergedNum)
        mergedResults[1] += ['NA']*statsNum
        return mergedResults

    if statsNum > 0:
        statsList = []
        for i in range(statsNum):
            statsList.append([])    # the list for storing each stats from each sub-window, NAs would also be stored
        snpNumList = []             # the list for storing SNP number from each sub-window

    ## Main loop for combining statistics and concatenate SNP positions
    for r, results in enumerate(resultList):
        # concatenate SNP positions
        mergedResults[0] += results[0]

        if statsNum > 0:
        # if there are statsitics to combine, construct statsList and snpNumeList
            for i, s in enumerate(results[1][1:]):
                statsList[i].append(s)
                # note that NAs could present in this list, will deal with this issue when specific merging method is specified later
            snpNumList.append(results[1][0])

        if len(results[2]) > 0:
        # if there is hiddens in this window
            for k in results[2]:
                if k not in hiddens:
                    hiddens[k] = []
                    if r > 0:
                    # fill up NAs for previous windows without hiddens
                        for i in range(r):
                            hiddens[k].append('NA')
                hiddens[k].append(results[2][k])    # fill up value for current window

        elif len(hiddens) > 0:
        # if no hiddens in this window but there were hiddens in previous window(s)
            for k in hiddens:
                hiddens[k].append('NA') # Add NA to fill up to correct records

    # Constructing the second sublist of mergedResults
    mergedResults[1].append(mergedNum)
    for i, h in enumerate(statsHeader):
    # loop through each statistics, determine which merging method to use, ignore the _p statistics, these will be dealt with in permPMerge
        # if all NA, output NA
        if statsList[i].count('NA') == len(statsList[i]):
            mergedResults[1].append('NA')
            continue

        if h == 'Allele_freq_diff':
            afd = weightedAverage(statsList[i], snpNumList)
            mergedResults[1].append(afd)
        
        elif h == 'dXY':
            dxy = simpleAverage(statsList[i])
            mergedResults[1].append(dxy)
    
        elif h == 'CSS_mod':
            css = simpleSum(statsList[i])
            mergedResults[1].append(css)

        elif h == 'FST':
            varcomps = hiddens['varcomp']
            totalvars = hiddens['totalvar']
            fst = simpleSum(varcomps) / simpleSum(totalvars)
            mergedResults[1].append(fst)

        elif h == 'Variance_comp':
            varcomp = simpleSum(statsList[i])
            mergedResults[1].append(varcomp)

    return mergedResults


##### Methods for merging statistics

## Average statistics weighted by SNP number
def weightedAverage(statsList, snpNumList):
    notNA = np.array(statsList) != 'NA'     # selection array indicating the element that is not NA
    fStatsList = np.array(np.array(statsList)[notNA], dtype=float)
    fSnpNumList = np.array(snpNumList)[notNA]

    weights = fSnpNumList / np.sum(fSnpNumList)
    weightedStats = fStatsList * weights

    return np.sum(weightedStats)


## Average statistics (simple mean) over all sub-windows, don't remove windows when no SNPs were present
def simpleAverage(statsList):
    nanStatsList = [np.nan if s=='NA' else s for s in statsList]
    return np.nansum(nanStatsList) / len(statsList)

## Sum up statistics over all sub-windows, ignore the windows with no SNPs
def simpleSum(statsList):
    nanStatsList = [np.nan if s=='NA' else s for s in statsList]
    return np.nansum(nanStatsList)


#### Method for permutating multiple windows
## ABANDONED due to potential non-randomness
## Permutation test on selected statistics but for inherent method, report the permutated sets and scores for inheritance, and the p-value
def permuPMerge(resultList, resultHeader, mergedResults, maxPerm, maxExtm, minExtm=0.5, *args, **kwargs):
# Return the p-value of the given score using a permutation test
# resultList: a list of results from each subwindow: [[result1], [result2] ... ], each result follows the structure: [[snppos], [snp_number, stats ... ], [hiddens], [genoArray]]
# resultHeader: a list of strings of the headers of the result output file
# mergedResults: a list of merged result following the structure [[snppos], [snp_number, stats ... ]]
# maxPerm: maximum number of permutations
# maxExtm: output the p-value if this number of score that is more extreme than the given score is observed
# inheritSets: inherited sets of groupIDs for each permutation from previous permuPMerge, 0 if none (the first window)
# inheritScores: inherited statistics scores for each permutation from previous permuPMerge, 0 if none (the first window)
# minExtm: minimum number of extreme values if no extreme values were observed in permutations, minExtm / maxPerm = the given p-value of this situation
    # Extracting the statistics name need to be calculated
    stats = []      # a list of strings, e.g., ['AFD_p', 'CSSm_p']
    for i, h in enumerate(resultHeader):
        if '_p' in h:
            stats.append(h)         # a list of strings, e.g., ['AFD_p', 'CSSm_p']

    # Thresholds for halting the permutation
    extremeNum = [0] * len(stats)   # number of extreme values observed in permutations
    returnP = [-1] * len(stats)     # the p-values that is going to be returned

    # Concatenate groupIDs and get sample number in each group
    decompGroupIDs = kwargs['groupIDs'][0] + kwargs['groupIDs'][1]
    g0Len = len(kwargs['groupIDs'][0])
    g1Len = len(kwargs['groupIDs'][1])

    # Extract genoArrays from each subwindow
    genoArrays = [r[3] for r in resultList]
    genoArrays = [i for i in genoArrays if type(i)!=int]      # remove all of the 0s (no variants)
    ctnGenoArray = al.GenotypeArray(np.concatenate(genoArrays))

    # Empirical merged scores
    empScores = mergedResults[1][1:]

    # Permutation loop starts
    nPerm = 0
    while True:
        nPerm += 1

        g0IDs = random.sample(decompGroupIDs, g0Len)
        g1IDs = list(set(decompGroupIDs) - set(g0IDs))
        newGroupIDs = [g0IDs, g1IDs]    # the permutated ID sets

        # Test if at least one of the group is completely missing
        isMiss, newGenoArray, keptList = ws.genoArrayMiss(ctnGenoArray, newGroupIDs, threshold=1.0)
        if isMiss:
            # reset if missed
            continue

        ## Call the divergence score functions
        hiddens = [[], []]    # in case of calculating FST, this will be [[varcomps ... ], [totalvars ... ]]

        # Loop the statistics first as it's easier to combine them
        for s, stat in enumerate(stats):
            # Check if the statistics still need to be calculated
            if extremeNum[s] >= maxExtm:
                continue

            ## Calculate statistics on concatenated windows
            if stat == 'AFD_p':
                permScore = ws.meanAlleleFreqDiffAllel(newGenoArray, newGroupIDs)
            elif stat == 'FST_p':
                permScore, varcomp, totalvar = ws.wcFst(newGenoArray, newGroupIDs)
                hiddens[0].append(varcomp)
                hiddens[1].append(totalvar)
            elif stat == 'VRC_p':
                if len(hiddens[0]) != 0:
                # if varcomp already calculated when dealing with FST
                    permScore = hiddens[0][-1]
                else:
                    fst, permScore, totalvar = ws.wcFst(newGenoArray, newGroupIDs)
            elif stat == 'CSSm_p':
                permScore = ws.css(newGenoArray, newGroupIDs)
            elif stat == 'dXY_p':
                permScore = ws.dxy(newGenoArray, newGroupIDs, kwargs['wSize'])

            ## Check if the permutation score is more extreme than empirical score
            if permScore >= empScores[s]:
            # if the permutation score is as extreme as or more extreme than the observed score
                extremeNum[s] += 1

                if extremeNum[s] == maxExtm:
                # if the number of extreme score is equal to the set threshold, return the p-value of the permutation test
                    returnP[s] = maxExtm / nPerm

                    if returnP.count(-1) == 0:
                    # If all test reach maximum extreme values
                        return returnP

        ## If maximum number of permutation is reached
        if nPerm == maxPerm:
            for i, p in enumerate(returnP):
                if p == -1:
                    en = extremeNum[i]
                    if en == 0:
                    # if no extreme values were observed, given p as there's minExtm extreme value observed
                        returnP[i] = minExtm / maxPerm
                    else:
                        returnP[i] = en / maxPerm
            return returnP
