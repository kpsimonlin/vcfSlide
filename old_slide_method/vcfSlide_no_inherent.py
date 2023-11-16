## This program runs sliding windows on a given VCF file and applies statistics on them, using pysam package.
## Author: Kung-Ping Lin; last editted: 20231031.

import argparse
from pysam import VariantFile
import os
import sys
import random
import window_stats as ws
import time

#### Argument parsing
parser = argparse.ArgumentParser(
                description='Iterate through user-defined sliding windows on a VCF file and apply statistics on them.')

## Required arguments
parser.add_argument('-v', '--vcf', help='The path of the input VCF file', required=True)
parser.add_argument('-w', '--window-size', help='The size of the window', required=True, type=int)
parser.add_argument('-g', '--gap', help='The gap for every window slide', required=True, type=int)
parser.add_argument('-o', '--output-prefix', help='The path and prefix of output files', required=True)

## General SNP information
parser.add_argument('--snp-position', help='Output the positions of SNPs in every window', action='store_true', default=False)
#parser.add_argument('--snp-number', help='Output the number of SNPs in every window', action='store_true', default=False)

## Divergence score
parser.add_argument('-G', '--group', help='The path to two files defining two groups of individual to compared, required for divergence statistics', nargs=2, default=0)
parser.add_argument('--allele-freq-diff', help='Output the mean ALT allele frequency difference between groups in each window', action='store_true', default=False)
parser.add_argument('--fst', help='Output Weir-Cockerham FST score between groups in each window', action='store_true', default=False)
parser.add_argument('--var-comp', help='Output the numerator of the  Weir-Cockerham FST score between groups in each window', action='store_true', default=False)
parser.add_argument('--css', help='Output the Cluster Separation Score (CSS) by Jones et al. (2012) in each window', action='store_true', default=False)
parser.add_argument('--css-mod', help='Output the CSS-modified, which is similar to CSS but without MDS', action='store_true', default=False)
parser.add_argument('--dxy', help='Output the dXY score in each window', action='store_true', default=False)

## SNP filtering options
parser.add_argument('--missing', help='The missingness threshold for both groups; default=0.9', type=float, default=0.9)
#parser.add_argument('--maf', help='The minor allele frequency threshold for both groups; default=0.1', type=float, default=0.1)

## Permutation settings
parser.add_argument('-P', '--permutation', help='Whether to perform permutation and output p-values', action='store_true', default=False)
parser.add_argument('--max-permutation', help='The maximum permutation number to perform; default=100000', type=int, default=100000)
parser.add_argument('--max-extreme', help='The maximum number of permutation score equal to or more extreme than the observed value before calculating p-value; default=10', type=int, default=10)
parser.add_argument('--min-extreme', help='The number of extreme values when there is no extreme value after maximum number of permutations; default=0.5', type=float, default=0.5)

## Inherent settings
parser.add_argument('--inherent', help='Enable inheritance to accelerate calculation' action='store_true', default=False)

args = parser.parse_args()

## Check if --group present when group statistics are in arguments
if args.group == 0:
    if args.allele_freq_diff:
        parser.error('--allle-freq-diff requires --group')
    if args.fst:
        parser.error('--fst requires --group')
    if args.var_comp:
        parser.error('--var-comp requires --group')
    if args.css:
        parser.error('--css requires --group')
    if args.css_mod:
        parser.error('--css-mod requires --group')
    if args.dxy:
        parser.error('--dxy requires --group')
    if args.permutation:
        parser.error('--permutation requires --group')

## Check if --inherent has conficts
if args.inherent:
    if args.fst:
        parser.error('--fst cannot be applied with --inherent')
    if args.window_size % gap != 0:
        parser.error('The remainder of window size devided by gap must be zero when applying --inherent')

vcfPath = args.vcf
wSize = args.window_size
gap = args.gap
outPrefix = args.output_prefix
missThreshold = args.missing
#mafThreshold = args.maf

## Print argument inputs
print('Arguments:')
for arg in vars(args):
    print(arg, getattr(args, arg))
print('')


#### Hidden parameters
windowReportSize = 1      # the denominator of reporting window scanned


#### Function for writing output files for statistics
## Open output files
writeResult = open(outPrefix + '.result.tsv', 'w')
resultHeader = []

if args.snp_position:
    writeSP = open(outPrefix + '.snppos.tsv', 'w')

#if args.snp_number:
#    resultHeader.append('SNP_number')
resultHeader.append('SNP_number')
if args.allele_freq_diff:
    resultHeader.append('Allele_freq_diff')
if args.fst:
    resultHeader.append('FST')
if args.var_comp:
    resultHeader.append('Variance_comp')
if args.css:
    resultHeader.append('CSS')
if args.css_mod:
    resultHeader.append('CSS_mod')
if args.dxy:
    resultHeader.append('dXY')

if args.permutation:
    if args.allele_freq_diff:
        resultHeader.append('AFD_p')
    if args.fst:
        resultHeader.append('FST_p')
    if args.var_comp:
        resultHeader.append('VRC_p')
    if args.css:
        resultHeader.append('CSS_p')
    if args.css_mod:
        resultHeader.append('CSSm_p')
    if args.dxy:
        resultHeader.append('dXY_p')
    
writeResult.write('\t'.join(resultHeader) + '\n')


## When skipping the divergence score calculation part
def skipLoop(results, resultHeader):
    naNum = len(resultHeader) - 1
#    if args.snp_number:
#        naNum = len(resultHeader) - 1
#    else:
#        naNum = len(resultHeader)
    for i in range(naNum):
        results.append('NA')
    writeResult.write('\t'.join(results) + '\n') 
    return 0


## Calculate statistics for records in each window
def windowStats(windowVcf, groupIDs, wSize, resultHeader):
    results = {}    # dictionary for storing results
    # Try to expand the iterator into lists of pos and genotypes
    poss, gnts = ws.expandIterator(windowVcf)   
    # the returned list will be empty if there is no variants in this window

    ## If output SNP positions
    if args.snp_position:
        if len(poss) == 0:
        # if there is no record
            writeSP.write('NA\n')
        else:
            writeSP.write('\t'.join([str(x) for x in poss]) + '\n')
    
    ## if output SNP number
    #if args.snp_number:
    num = ws.snpNumber(gnts)
    results.append(str(num))
    if num == 0:
    # if there is no record
        skipLoop(results, resultHeader)
        return 0

    if groupIDs != 0:
        ## Convert records to genotypeArray object
        genoArray = ws.toGenotypeArray(gnts)

        ## Check if only invariants remain, if not, remove the invariant sites
        ifAllInvar, genoArray = ws.genoArrayInvar(genoArray, remove=True)
        if ifAllInvar:
            skipLoop(results, resultHeader)
            return 0

        ## Check if the variants are missed for a certain threshold in either group
        isMiss, genoArray = ws.genoArrayMiss(genoArray, groupIDs, threshold=missThreshold)
        if isMiss:
            skipLoop(results, resultHeader)
            return 0

        ## ABANDONED. Check if the MAF of the variants are lower than a specified threshold in either group
        #lowMAF, genoArray = ws.genoArrayMAF(genoArray, groupIDs, threshold=mafThreshold)
        #if lowMAF:
        #    skipLoop(results, resultHeader)
        #    return 0

        if args.permutation:
            stats = []      # a list of the function for permutation test
            scores = []     # a list of the divergence scores for permutation test

        ## If output mean allele frequency difference
        if args.allele_freq_diff:
            afd = ws.meanAlleleFreqDiffAllel(genoArray, groupIDs)
            results.append(str(afd))
            if args.permutation:
                stats.append(ws.meanAlleleFreqDiffAllel)
                scores.append(afd)

        ## If output FST
        if args.fst:
            fst = ws.wcFst(genoArray, groupIDs)
            results.append(str(fst))
            if args.permutation:
                stats.append(ws.wcFst)
                scores.append(fst)

        ## if output variance component
        if args.var_comp:
            varcomp = ws.wcVarComp(genoArray, groupIDs)
            results.append(str(varcomp))
            if args.permutation:
                stats.append(ws.wcVarComp)
                scores.append(varcomp)

        ## If output CSS
        if args.css:
            c = ws.css(genoArray, groupIDs)
            results.append(str(c))
            if args.permutation:
                stats.append(ws.css)
                scores.append(c)

        ## If output CSS-modified
        if args.css_mod:
            cm = ws.cssMod(genoArray, groupIDs)
            results.append(str(cm))
            if args.permutation:
                stats.append(ws.cssMod)
                scores.append(cm)

        ## If output dXY
        if args.dxy:
            dxy = ws.dxy(genoArray, groupIDs, wSize)
            results.append(str(dxy))
            if args.permutation:
                stats.append(ws.dxy)
                scores.append(dxy)

        ## If doing permuation test
        if args.permutation:
            ps = ws.permuP(scores, stats, args.max_permutation, args.max_extreme, totalWindow, args.min_extreme, genoArray=genoArray, groupIDs=groupIDs, wSize=wSize)
            results += [str(x) for x in ps]

        ## Write results
        writeResult.write('\t'.join(results) + '\n')

    return 0


## Write the results into output files
def writeResults(results, resultHeader):
        if len(poss) == 0:
        # if there is no record
            writeSP.write('NA\n')
        else:
            writeSP.write('\t'.join([str(x) for x in poss]) + '\n')

        skipLoop(results, resultHeader)
        ## Write results
        writeResult.write('\t'.join(results) + '\n')

## Write the start and end position of the window
writeWSE = open(outPrefix + '.winpos.tsv', 'w')
def writeWindowStartEnd(start, wSize):
# this function assumes the given start position is 1-based
    writeWSE.write(str(start) + '\t' + str(start + wSize - 1) + '\n')
    return 0

## Set timer
sTime = time.time()

#### Main loop for sliding window
## pysam reader
vcf = VariantFile(vcfPath, 'r')
#vcfReader = vcf.Reader(filename=vcfPath)

## Keep only the interested individuals defined in groups, if applicable
if args.group != 0:
    with open(args.group[0], 'r') as readG1:
        group1Ind = readG1.read().splitlines()      # group1 sample names
    with open(args.group[1], 'r') as readG2:
        group2Ind = readG2.read().splitlines()      # group2 sample names
    vcf.subset_samples(group1Ind + group2Ind)

## Construct the looping lists for sliding windows
print('Prepare window sliding ...')

# Find the base position of the first SNP
firstV = next(vcf)
posStart = firstV.pos    # note that this value is 1-based
# Find the current chromosome, assuming the vcf only contains one contig, for later fetching
chrom = firstV.chrom
# Find the indices for the two groups, if applicable
if args.group != 0:
    groupIDs = ws.getGroupID(list(firstV.samples), group1Ind, group2Ind)
# Find the base position of the last SNP, this is going to take a while
for r in vcf:
    pass
posEnd = r.pos              # note that this value is 1-based

# Reset the iterator to the first base
vcf.reset()

# Compose the list for each sub-sliding windows for inherent method, note that these values are converted to 0-based, half-opened coordinates
if args.inherent:
    fetchCoord = []
    i = posStart - 1
    while i + gap <= posEnd:
        fetchCoord.append((i, i + gap))
        i += gap
    

# Compose the list for fetching the sliding windows for normal method, note that these values are converted to 0-based, half-opened coordinates
# fetchCoord = [(s1, e1), (s2, e2), ... ]
else:
    fetchCoord = []
    i = posStart - 1
    while i + wSize <= posEnd:
        fetchCoord.append((i, i + wSize))
        i += gap



## Define storing values
totalWindow = 0     # store the total number of window scanned

## Loop through SNPs in the VCF
print('Time spent: ' + str(round(time.time() - sTime, 2)) + ' sec; window sliding began.')
if args.inherent:


else:
    for s, e in fetchCoord:
        # fetch window vcf   
        windowVcf = vcf.fetch(contig=chrom, start=s, end=e)
        # parse the record in writeRecords function
        writeRecords(windowVcf, groupIDs, wSize=wSize, resultHeader=resultHeader, totalWindow=totalWindow)
        writeWindowStartEnd(s + 1, wSize)

        # output the current number of window processed
        totalWindow += 1
        if totalWindow % windowReportSize == 0:
            print(str(totalWindow) + ' windows scanned; time spent: ' + str(round(time.time() - sTime, 2)) + ' sec.')


#### Close openned files
writeWSE.close()
writeResult.close()
if args.snp_position:
    writeSP.close()


#### Final output
print('Finished; total window scanned: ' + str(totalWindow) + '; total time spent: ' + str(round(time.time() - sTime, 2)) + ' sec.')
