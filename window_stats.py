## Module for statistics funcitons applied in a single sliding window
## Kung-Ping Lin, 20221110.
## All functions take a list of records (PyVCF object) in the sliding window as input

import numpy as np

## Return the positions of SNPs in the window
def snpPositions(records):
    poss = [r.POS for r in records]
    return poss


## Return the number of SNPs in the window
def snpNumber(records):
    num = len(records)
    return num


## The following two functions return two lists of genotypes from two groups defined by two group files
# return the indices of the group individuals in the VCF sample names
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

# return the genotypes of the group individuals
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


## Return the mean ALT allele frequency difference between two groups in the window
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
    
