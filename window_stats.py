## Module for statistics funcitons applied in a single sliding window
## Kung-Ping Lin, 20221110.
## All functions take a list of records (PyVCF object) in the sliding window as input

import numpy as np
import allel as al
import re
import scipy
import miscell
import random
import time
from itertools import combinations
from scipy.spatial.distance import pdist, cdist

## Return the positions of SNPs in the window
def snpPositions(records):
    poss = [r.POS for r in records]
    return poss


## Return the number of SNPs in the window
def snpNumber(records):
    num = len(records)
    return num


## Return the indices of the group individuals in the VCF sample names
def getGroupID(record, group):
    # get sample names for two defined groups
    with open(group[0], 'r') as readG1:
        group1Ind = readG1.read().splitlines()     # list of sample names for group1
    with open(group[1], 'r') as readG2:
        group2Ind = readG2.read().splitlines()     # list of sample names for group2

    # get sample names for the VCF
    vcfInd = [s.sample for s in record.samples]
    
    # get the indices of the group individuals in the VCF sample names
    group1Id = []
    group2Id = []
    for i, n in enumerate(vcfInd):
        if n in group1Ind:
            group1Id.append(i)
        elif n in group2Ind:
            group2Id.append(i)

    return [group1Id, group2Id]


## Transfer the records in the window to a GenotypeArray object in scikit-allel package, ready for some divergence statistics calculation
def toGenotypeArray(records):
    recordList = []
    # will be filled with a three-dimentional list
    # 1st dimention: variants
    # 2nd dimention: individuals
    # 3rd dimention: ploidy
    # example: a window of 4 variants, 3 individuals, and ploidy of 2 will be
    # [[[0, 0], [0, 1], [0, 1]],
    #  [[0, 1], [0, 0], [1, 1]],
    #  [[0, 0], [1, 1], [0, 1]],
    #  [[0, 0], [0, 1], [0, 1]]]

    for record in records:
    # 1st loop: SNP
        snpGenos = []

        for s in record.samples:
        # 2nd loop: sample
            geno = []

            for g in re.split('/|\|', s['GT']):     # split with / or |
            # 3rd loop: ploidy
                if g == '.':
                # transfer missing genotype to -1
                    geno.append(-1)
                else:
                    geno.append(int(g))

            snpGenos.append(geno)

        recordList.append(snpGenos)

    return al.GenotypeArray(recordList)


## Return a subset of the genoArray containing only the groups interested, return a new groupID list
def subsetGenotypeArray(genoArray, groupIDs):
    newGenoArray = genoArray.subset(sel1=groupIDs[0] + groupIDs[1])
    newGroupIDs = [list(range(len(groupIDs[0]))), list(range(len(groupIDs[0]), len(groupIDs[0]) + len(groupIDs[1])))]
    return [newGenoArray, newGroupIDs]


## Check if genotypeArray contains invariants, remove them if specified
def genoArrayInvar(genoArray, remove=True):
    ac = genoArray.count_alleles()
    isNv = ac.is_non_segregating()
    if all(isNv):
        return [True, 0]

    if remove:
        removedArray = genoArray[[not x for x in isNv]]

        return [False, removedArray]
    else:
        return[False, genoArray]


## Filter the per-group missing rate of a genotypeArray and it's associated groupID
def genoArrayMiss(genoArray, groupIDs, threshold=1.0):
# default: remove variants that are all missings in either group
    removeVar = []  # list for the variant ID to be removed
    for groupID in groupIDs:
        groupGenoArray = genoArray.subset(sel1=groupID)
        for i, var in enumerate(groupGenoArray.is_missing()):
        # loop through an array of missing value boolean flags
            freqMiss = np.count_nonzero(var) / len(var)
            if freqMiss >= threshold:
                if i not in removeVar:
                    removeVar.append(i)

    if len(removeVar) == len(groupGenoArray):
    # if all of the variants are missing in either group
        return [True, 0]

    keptVar = []
    for i in range(len(groupGenoArray)):
        if i not in removeVar:
            keptVar.append(i)

    return [False, genoArray.subset(sel0=keptVar)]


## Test if the minor allele frequencies in both groups are below a certain threshold
def genoArrayMAF(genoArray, groupIDs, threshold):
# default: remove variants that are all missings in either group
    removeVar = []  # list for the variant ID to be removed
    for groupID in groupIDs:
        groupGenoArray = genoArray.subset(sel1=groupID)
        for i, ac in enumerate(groupGenoArray.count_alleles()):
        # loop through an array of allele counts for each variant
            freq = min(ac) / np.sum(ac)
            if freq < threshold:
                if i not in removeVar:
                    removeVar.append(i)

    if len(removeVar) == len(groupGenoArray):
    # if all of the variants are filtered out in either group
        return [True, 0]

    keptVar = []
    for i in range(len(groupGenoArray)):
        if i not in removeVar:
            keptVar.append(i)

    return [False, genoArray.subset(sel0=keptVar)]
        

## Return the mean ALT allele frequency difference between two groups in the window, calculated by allel package.
def meanAlleleFreqDiffAllel(genoArray, groupIDs, wSize=None):
    meanAlleleFreqs = []
    for groupID in groupIDs:
        alleleCount = genoArray.count_alleles(subpop=groupID)
        totalCount = np.sum(alleleCount)
        altCount = np.sum(alleleCount[:, 1])
        meanAlleleFreqs.append(altCount / totalCount)

    return abs(meanAlleleFreqs[0] - meanAlleleFreqs[1])


## Calculate FST score between two groups
def wcFst(genoArray, groupIDs, wSize=None):
    a, b, c = al.weir_cockerham_fst(genoArray, groupIDs)
    fst = np.sum(a) / (np.sum(a) + np.sum(b) + np.sum(c))
    return fst


## Calculate variace component between two groups as a absolute divergence statistics
def wcVarComp(genoArray, groupIDs, wSize=None):
    a, b, c = al.weir_cockerham_fst(genoArray, groupIDs)
    return np.sum(a)


## Calculate CSV score between two groups
## The Cluster Separation Score (CSS) measures genetic difference between all pairs of individuals as euclidean distance in two dimensions obtained from an MDS analysis of proportion sequence divergence.
## The CSS score is the average distance between pairs of individuals belonging to different populations minus the average distance between pairs belonging to the same populations
def css(genoArray, groupIDs, wSize=None):
# note that the genoArray here only contains the genotypes from the group interested, and groupIDs are adjusted.
    # Adjust the genoArray to 0, 0.5, 1 genotypes; -1 for missing values
    proGenoArray = genoArray.to_n_alt(fill=-2)/2

    # Fill missing genotype to the mean of the group
    for v in range(genoArray.n_variants):
        group0Geno = proGenoArray[v][groupIDs[0]]
        group0Mean = np.mean(group0Geno[group0Geno >= 0.0])
        group1Geno = proGenoArray[v][groupIDs[1]]
        group1Mean = np.mean(group1Geno[group1Geno >= 0.0])
        for s in groupIDs[0]:
            if proGenoArray[v, s] == -1.0:
            # if missing
                proGenoArray[v, s] = group0Mean
        for s in groupIDs[1]:
            if proGenoArray[v, s] == -1.0:
            # if missing
                proGenoArray[v, s] = group1Mean

    # Calculate pairwise sequence distances, scaled to proportion sequence divergence when we adjusted the genoArray to 0, 0.5, 1 on the above step
    psd = al.pairwise_distance(proGenoArray, metric='cityblock')

    # Perform pcoa, a.k.a. classical multi-dimensional scaling; select the first 2 axis
    mds = al.pcoa(psd)[0][:, [0, 1]]

    # Euclidean distance matrix calculation on the transposed MDS
    mdsG0 = mds[groupIDs[0], :]
    mdsG1 = mds[groupIDs[1], :]
    s00 = np.sum(pdist(mdsG0, metric='euclidean'))
    s11 = np.sum(pdist(mdsG1, metric='euclidean'))
    s01 = np.sum(cdist(mdsG0, mdsG1, metric='euclidean'))

    m = len(groupIDs[0])
    n = len(groupIDs[1])

    css = ( s01 / (m*n) ) - ( 1/(m+n) ) * ( s00/( (m-1)/2 ) + s11/( (n-1)/2 ) )

    return css


## Calculate dXY score between two groups
## dXY defines the pairwise nucleotide distance between the marine and freshwater group.
def dxy(genoArray, groupIDs, wSize):
    dxy = np.sum(al.mean_pairwise_difference_between(genoArray.count_alleles(subpop=groupIDs[0]), genoArray.count_alleles(subpop=groupIDs[1]))) / wSize
    return dxy

## Permutation test on selected statistics, report the p-values
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

    genoArray = kwargs['genoArray']
    decompGroupIDs = kwargs['groupIDs'][0] + kwargs['groupIDs'][1]
    g0Len = len(kwargs['groupIDs'][0])
    g1Len = len(kwargs['groupIDs'][1])

    nPerm = 0
    while True:
        if nPerm == maxPerm:
            break

        print(nPerm)
        print(extremeNum)
        print(returnP)
        print(time.time() - sTime)
        print()
        g0IDs = random.sample(decompGroupIDs, g0Len)
        g1IDs = list(set(decompGroupIDs) - set(g0IDs))
        newGroupIDs = [g0IDs, g1IDs]

        # Test if at least one of the group is completely missing
        isMiss, newGenoArray = genoArrayMiss(genoArray, newGroupIDs, threshold=1.0)
        if isMiss:
            continue
        else:
            kwargs['genoArray'] = newGenoArray
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
                        returnP[s] = maxExtm / (nPerm + 1)

                        if returnP.count(-1) == 0:
                        # If all test reach maximum extreme values
                            return returnP

        nPerm += 1

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
