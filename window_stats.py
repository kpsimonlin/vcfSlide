## Module for statistics funcitons applied in a single sliding window
## Kung-Ping Lin, 20221110.
## All functions take a list of records (PyVCF object) in the sliding window as input

import numpy as np
import allel as al
import re

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
    if len(records) == 0:
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

## OLD: Return the genotypes of the group individuals
def getGroupGenos(records, groupIds):
# group is a list containing two paths to two files recording two groups of individual IDs
    if len(records) == 0:
    # return 0 if no SNPs in window
        return 0

    group1Id = groupIds[0]
    group2Id = groupIds[1]

    group1Geno = []     # list of group1 genotypes [['0/0', '0/1', ... ], [ ... ], ... ]
    group2Geno = []     # list of group2 genotypes
    for record in records:
        # list for group1 genotypes in the current SNP ['0/0', '0/1', ... ]
        group1SnpGeno = [record.samples[i]['GT'] for i in group1Id]
        group1Geno.append(group1SnpGeno)

        # list for group2 genotypes in the current SNP
        group2SnpGeno = [record.samples[i]['GT'] for i in group2Id]
        group2Geno.append(group2SnpGeno)

    return [group1Geno, group2Geno]


## Return the mean ALT allele frequency difference between two groups in the window, calculated by allel package.
def meanAlleleFreqDiffAllel(genoArray, groupIDs):
    if type(genoArray) == int:
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

## OLD: Return the mean ALT allele frequency difference between two groups in the window
def meanAlleleFreqDiff(records, groupsGeno):
# groupsGeno is a list of two lists containing the genotypes of each group
# ignore missing genotypes - they are not included in the calculation
    if groupsGeno == 0:
    # return 'NA' if no SNP is in window
        return 'NA'

    groupAltFreqs = []
    for gGeno in groupsGeno:
    # gGeno = group genotypes, should be only two groups
        count0 = 0  # count of REF allele
        count1 = 0  # count of ALT allele
        for snpGeno in gGeno:
        # snpGeno = genotypes of individuals in a group
            for ind in snpGeno:
            # ind = individual genotype, e.g., '0/0'
                count0 += ind.count('0')
                count1 += ind.count('1')
        
        # if all missing, return NA (this happens because of group individual extraction)
        if count0 + count1 == 0:
            return 'NA'
        # calculate mean allele frequency
        groupAltFreq = count1/(count0 + count1)
        groupAltFreqs.append(groupAltFreq)

    return abs(groupAltFreqs[0] - groupAltFreqs[1])


## Calculate FST score between two groups
def wcFst(genoArray, groupIDs, reportA=False):
# reportA: whether to report the variance component a (among groups) as a absolute divergence statistics
    if type(genoArray) == int:
    # if no records in this window
        if reportA:
            return ['NA', 'NA']
        else:
            return 'NA'

    a, b, c = al.weir_cockerham_fst(genoArray, groupIDs)

    if np.sum(a) + np.sum(b) + np.sum(c) == 0.0:
    # if only invariants present in the window
        if reportA:
            return ['NA', 'NA']
        else:
            return 'NA'

    fst = np.sum(a) / (np.sum(a) + np.sum(b) + np.sum(c))

    if reportA:
        return [fst, np.sum(a)]
    else:
        return [fst]

