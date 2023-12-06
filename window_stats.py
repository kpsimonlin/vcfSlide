## Module for statistics funcitons applied in a single sliding window
## Kung-Ping Lin, 20221110.
## All functions take a list of records (PyVCF object) in the sliding window as input

import numpy as np
import allel as al
import re
import scipy
import random
import time
from itertools import combinations
from scipy.spatial.distance import pdist, cdist

## Expand the iterator into lists of positions and genotypes, for the genotypes replace the missing values with -1 to get ready for the conversion to GenotypeArray
def expandIterator(windowVcf):
    try:
        poss = []
        gnts = []
        for record in windowVcf:
            poss.append(record.pos)
            newGts = []
            for s in record.samples.values():
                gt = s['GT']    # this has structure like (0, 1) or (None, 1) for which None is a missing value
                newGt = [-1 if g is None else g for g in gt]    # If there is None value, replace it with -1
                newGts.append(newGt)
            gnts.append(newGts)
            # will be filled with a three-dimentional list
            # 1st dimention: variants
            # 2nd dimention: individuals
            # 3rd dimention: ploidy
            # example: a window of 4 variants, 3 individuals, and ploidy of 2 will be
            # [[[0, 0], [0, 1], [0, 1]],
            #  [[0, 1], [0, 0], [1, 1]],
            #  [[0, 0], [1, 1], [0, 1]],
            #  [[0, 0], [0, 1], [0, 1]]]
        return [poss, gnts]

    except StopIteration:
        # if there is no record in this window
        # return [[poss], [genotypes]]
        return [[], []]


## Return the number of SNPs in the window
def snpNumber(gnts):
    num = len(gnts)
    return num


## Return the indices of the group individuals in the VCF sample names
def getGroupID(samples, g1samples, g2samples):
    # list of all samples, list of group1 samples, list of group2 samples
    # get the indices of the group individuals in the VCF sample names
    group1Id = []
    group2Id = []
    for i, n in enumerate(samples):
        if n in g1samples:
            group1Id.append(i)
        elif n in g2samples:
            group2Id.append(i)

    return [group1Id, group2Id]


## Transfer the records in the window to a GenotypeArray object in scikit-allel package, ready for some divergence statistics calculation
def toGenotypeArray(gnts):
    return al.GenotypeArray(gnts)


## Check if genotypeArray contains invariants, remove them if specified
def genoArrayInvar(genoArray, remove=True):
    ac = genoArray.count_alleles()
    isNv = ac.is_non_segregating()
    if all(isNv):
        return [True, 0, 0]

    if remove:
        removedArray = genoArray[[not x for x in isNv]]

        return [False, removedArray, [not x for x in isNv]]
    else:
        return[False, genoArray, [not x for x in isNv]]


## Filter the per-group missing rate of a genotypeArray and it's associated groupID
def genoArrayMiss(genoArray, groupIDs, threshold=1.0):
# default: remove variants that are all missings in either group
    keptList = [True]*len(genoArray)   # list of true and false for which true is kept
    for groupID in groupIDs:
        groupGenoArray = genoArray.subset(sel1=groupID)
        for i, var in enumerate(groupGenoArray.is_missing()):
        # loop through an array of missing value boolean flags
            freqMiss = np.count_nonzero(var) / len(var)
            if freqMiss >= threshold:
                if keptList[i]:
                    keptList[i] = False

    if sum(keptList) == 0:
    # if all of the variants are missing in either group
        return [True, 0, 0]
    keptVar = genoArray[keptList]

    return [False, keptVar, keptList]


## Test if the minor allele frequencies in both groups are below a certain threshold
## This function is abandoned.
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
    alleleFreqs = []
    for groupID in groupIDs:
        alleleCount = genoArray.count_alleles(subpop=groupID)
        refCount = alleleCount[:, 0]
        altCount = alleleCount[:, 1]
        alleleFreq = altCount / (refCount + altCount)
        alleleFreqs.append(alleleFreq)
    
    alleleDiffs = abs(alleleFreqs[0] - alleleFreqs[1])
    meanAlleleDiff = np.mean(alleleDiffs)

    return meanAlleleDiff


## Calculate FST score between two groups
def wcFst(genoArray, groupIDs, wSize=None):
    a, b, c = al.weir_cockerham_fst(genoArray, groupIDs)
    bgVar = np.sum(a)
    tVar = np.sum(a) + np.sum(b) + np.sum(c)
    fst = bgVar / tVar
    return [fst, bgVar, tVar]


## Calculate CSS score between two groups
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


## Calculate CSS-modified, which use the distance matrix itself without performing MDS.
def cssMod(genoArray, groupIDs, wSize=None):
# The only difference between this method and the above CSS is that MDS is skipped in this method, and the distance matrix is put directly into the calculation of the score, which is supposed to be faster if a lot of samples were put into the calculation.
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

    # Define the genotype matrices of the two groups, transpose them to fit in python pdist and cfist calculation
    g0GenoArray = np.transpose(proGenoArray[:, groupIDs[0]])
    g1GenoArray = np.transpose(proGenoArray[:, groupIDs[1]])

    # Calculate the sum of within- and between- group distances
    s00 = np.sum(pdist(g0GenoArray, metric='cityblock'))
    s11 = np.sum(pdist(g1GenoArray, metric='cityblock'))
    s01 = np.sum(cdist(g0GenoArray, g1GenoArray, metric='cityblock'))

    p1Time = time.time()

    # Calculate CSS
    m = len(groupIDs[0])
    n = len(groupIDs[1])

    css = ( s01 / (m*n) ) - ( 1/(m+n) ) * ( s00/( (m-1)/2 ) + s11/( (n-1)/2 ) )

    return css


## Calculate dXY score between two groups
## dXY defines the pairwise nucleotide distance between the marine and freshwater group.
def dxy(genoArray, groupIDs, wSize):
    dxy = np.nansum(al.mean_pairwise_difference_between(genoArray.count_alleles(subpop=groupIDs[0]), genoArray.count_alleles(subpop=groupIDs[1]))) / wSize
    return dxy


## Permutation test on selected statistics, report the p-values
def permuP(scores, stats, maxPerm, maxExtm, minExtm=0.5, *args, **kwargs):
# Return the p-value of the given score using a permutation test
# scores: a list of given divergence scores
# stats: a list of functions calculating the divergence score
# maxPerm: maximum number of permutations
# maxExtm: output the p-value if this number of score that is more extreme than the given score is observed
# minExtm: minimum number of extreme values if no extreme values were observed in permutations, minExtm / maxPerm = the given p-value of this situation
#    sTime = time.time()

    extremeNum = [0] * len(stats)   # number of extreme values observed in permutations
    returnP = [-1] * len(stats)     # the p-values that is going to be returned

    genoArray = kwargs['genoArray']
    decompGroupIDs = kwargs['groupIDs'][0] + kwargs['groupIDs'][1]
    g0Len = len(kwargs['groupIDs'][0])
    g1Len = len(kwargs['groupIDs'][1])

    nPerm = 0
    while True:
        if nPerm == maxPerm:
            break

        random.shuffle(decompGroupIDs)
#        g0IDs = decompGroupIDs[:g0Len]
#        g1IDs = decompGroupIDs[g0Len:]
        g0IDs = random.sample(decompGroupIDs, g0Len)
        g1IDs = list(set(decompGroupIDs) - set(g0IDs))
        newGroupIDs = [g0IDs, g1IDs]

        # Test if at least one of the group is completely missing
        isMiss, newGenoArray, trueList = genoArrayMiss(genoArray, newGroupIDs, threshold=1.0)
        if isMiss:
            continue
        else:
            kwargs['genoArray'] = newGenoArray
            kwargs['groupIDs'] = newGroupIDs

        # Call the divergence score function
        permScores = []
        for s, func in enumerate(stats):
            if func == 'afd':
                stat = meanAlleleFreqDiffAllel
            elif func == 'fst':
                stat = wcFst
            elif func == 'varcomp':
                stat = wcFst
            elif func == 'css':
                stat = css
            elif func == 'cssmod':
                stat = cssMod
            elif func == 'dxy':
                stat = dxy

            if extremeNum[s] < maxExtm:
                permScore = stat(*args, **kwargs)
                if func == 'fst':
                    permScore = permScore[0]
                elif func == 'varcomp':
                    permScore = permScore[1]

                permScores.append(permScore)

                if permScore >= scores[s]:
                # if the permutation score is as extreme as or more extreme than the observed score
                    extremeNum[s] += 1

                    if extremeNum[s] == maxExtm:
                    # If the number of extreme score is equal to the set threshold, return the p-value of the permutation test
                        returnP[s] = maxExtm / (nPerm + 1)

                        if returnP.count(-1) == 0:
                        # If all test reach maximum extreme values
                            return returnP

#        print(nPerm)
#        print(scores)
#        print(permScores)
#        print(extremeNum)
#        print(returnP)
#        print(time.time() - sTime)
#        print()

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


