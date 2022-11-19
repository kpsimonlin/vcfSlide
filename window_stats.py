## Module for statistics funcitons applied in a single sliding window
## Kung-Ping Lin, 20221110.
## All functions take a list of records (PyVCF object) in the sliding window as input

import numpy as np
import allel as al
import re
import scipy
import miscell

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
def toGenotypeArray(records, noRecord=False):
    if noRecord:
        # return 0 if no SNPs in window
        return 0

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
def subsetGenotypeArray(genoArray, groupIDs, noRecord=False):
    if noRecord:
        return [0, 0]

    newGenoArray = genoArray.subset(sel1=groupIDs[0] + groupIDs[1])
    newGroupIDs = [range(len(groupIDs[0])), range(len(groupIDs[0]), len(groupIDs[0]) + len(groupIDs[1]))]
    return [newGenoArray, newGroupIDs]


## Check if genotypeArray contains invariants, remove them if specified
def genoArrayInvar(genoArray, remove=True, noRecord=False):
    if noRecord:
        return [True, 0]

    ac = genoArray.count_alleles()
    isNv = ac.is_non_variant()
    if all(isNv):
        return [True, 0]

    if remove:
        removedArray = genoArray[[not x for x in isNv]]
        return [False, removedArray]
    else:
        return[False, genoArray]


## Return the mean ALT allele frequency difference between two groups in the window, calculated by allel package.
def meanAlleleFreqDiffAllel(genoArray, groupIDs, noRecord=False):
    if noRecord:
    # if no records in this window
        return 'NA'

    meanAlleleFreqs = []
    for groupID in groupIDs:
        alleleCount = genoArray.count_alleles(subpop=groupID)
        totalCount = np.sum(alleleCount)
        if totalCount == 0:
        # All missings for this subpopulation in this window
            return 'NA'

        try:
            altCount = np.sum(alleleCount[:, 1])
            meanAlleleFreqs.append(altCount / totalCount)
        except IndexError:
        # only invariants present in this window after filtering group individuals
            meanAlleleFreqs.append(0.0)
    if meanAlleleFreqs[0] == 0.0 and meanAlleleFreqs[1] == 0.0:
    # if no variants present after group individual filtering
        return 'NA'
    else:
        return abs(meanAlleleFreqs[0] - meanAlleleFreqs[1])


## Calculate FST score between two groups
def wcFst(genoArray, groupIDs, reportA=False, noRecord=False):
# reportA: whether to report the variance component a (among groups) as a absolute divergence statistics
    if noRecord:
    # if no records in this window
        if reportA:
            return ['NA', 'NA']
        else:
            return ['NA']

    a, b, c = al.weir_cockerham_fst(genoArray, groupIDs)

    fst = np.sum(a) / (np.sum(a) + np.sum(b) + np.sum(c))

    if reportA:
        return [fst, np.sum(a)]
    else:
        return [fst]


## Calculate CSV score between two groups
## The Cluster Separation Score (CSS) measures genetic difference between all pairs of individuals as euclidean distance in two dimensions obtained from an MDS analysis of proportion sequence divergence.
## The CSS score is the average distance between pairs of individuals belonging to different populations minus the average distance between pairs belonging to the same populations
def css(genoArray, groupIDs, noRecord=False):
# note that the genoArray here only contains the genotypes from the group interested, and groupIDs are adjusted.
    if noRecord:
        return 'NA'

    # Adjust the genoArray to 0, 0.5, 1 genotypes; -1 for missing values
    proGenoArray = genoArray.to_n_alt(fill=-1)/2

    # Fill missing genotype to the mean of the group
    for v in range(genoArray.n_variants):
        for s in range(genoArray.n_samples):
            if proGenoArray[v, s] == -1:
            # if missing
                if s < len(groupIDs[0]):
                # if in group0
                    proGenoArray[v, s] = np.mean(proGenoArray[v][groupIDs[0]])
                else:
                # if in group1
                    proGenoArray[v, s] = np.mean(proGenoArray[v][groupIDs[1]])

    # Calculate pairwise sequence distances, scaled to proportion sequence divergence when we adjusted the genoArray to 0, 0.5, 1 on the above step
    psd = al.pairwise_distance(proGenoArray, metric='cityblock')

    # Perform pcoa, a.k.a. classical multi-dimensional scaling
    mds = al.pcoa(psd)[0]
    # transpose and extract the first two axis of mds so that it's in the same format as genoArray
    tranMds = np.transpose(mds)[:2, ]
    
    # Euclidean distance matrix calculation on the transposed MDS
    eucd = al.pairwise_distance(tranMds, metric='euclidean')

    s00 = 0.0    # sum of distance for within-group0 pairs
    s01 = 0.0    # sum of distance for between-group pairs
    s11 = 0.0    # sum of distance for within-group1 pairs
    # Loop the matrix, define whether it's among or between populations
    for i, eu in enumerate(eucd):
        # get the individual IDs from each entry in the condensed distance matrix, here n = # samples
        x, y = miscell.condensed_to_square(i, n=len(tranMds[0]))
        if x < len(groupIDs) and y < len(groupIDs):
            s00 += eu
        elif (x < len(groupIDs[0]) and y >= len(groupIDs[0])) or (x >= len(groupIDs[0]) and y < len(groupIDs[0])):
            s01 += eu
        else:
            s11 += eu

    m = len(groupIDs[0])
    n = len(groupIDs[1])

    css = ( s01 / (m*n) ) - ( 1/(m+n) ) * ( s00/( (m-1)/2 ) + s11/( (n-1)/2 ) )

    return css
