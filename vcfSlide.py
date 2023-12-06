## This program runs sliding windows on a given VCF file and applies statistics on them, using pysam package.
## Author: Kung-Ping Lin; last editted: 20231031.

import argparse
from pysam import VariantFile
import numpy as np
import os
import sys
import random
import window_stats as ws
import window_merge as wm
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
parser.add_argument('--inherent', help='Enable inheritance to accelerate calculation', action='store_true', default=False)

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
    if args.css:
        parser.error('--css cannot be applied with --inherent, use --css-mod instead')
    if args.window_size % args.gap != 0:
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
writeWSE = open(outPrefix + '.winpos.tsv', 'w')

if args.snp_position:
    writeSP = open(outPrefix + '.snppos.tsv', 'w')

resultHeader = ['SNP_number']
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


## Calculate statistics for records in each window
def windowStats(windowVcf, groupIDs, wSize, resultHeader):
# return a list of [snp_number, results, hiddens]. snp_number default to [], remains [] if no record or not recording
    results = [0] + ['NA']*(len(resultHeader) - 1)
    # results = [0, 'NA', 'NA', ... ]. 0 for default SNP number, NAs for default statistics
    hiddens = {}
    # hiddens is a dictionary storing hidden statistics needed for the inherent method
    inhGArray = 0
    # the genoArray for inherent method, 0 if not using inherent method

    ## Try to expand the iterator into lists of pos and genotypes
    poss, gnts = ws.expandIterator(windowVcf)   
    # the returned list will be empty if there is no variants in this window

    ## Return default results if there is no record
    if len(poss) == 0:
        return [[], results, hiddens, inhGArray]

    else:
        ## Report the position if specified, otherwise []
        if args.snp_position:
            snpPoss = poss
        else:
            snpPoss = []

        ## If calculating divergence statistics
        if groupIDs != 0:
            ## Convert records to genotypeArray object
            genoArray = ws.toGenotypeArray(gnts)

            ## Check if only invariants remain, if not, remove the invariant sites
            ifAllInvar, genoArray, varList = ws.genoArrayInvar(genoArray, remove=True)
            if ifAllInvar:
                return [[], results, hiddens, inhGArray]
            else:
                # update SNP number and SNP positions
                if args.snp_position:
                    snpPoss = list(np.array(snpPoss)[varList])
                results[0] = len(genoArray)

            ## Check if the variants are missed for a certain threshold in either group
            isMiss, genoArray, varList = ws.genoArrayMiss(genoArray, groupIDs, threshold=missThreshold)
            if isMiss:
                return [[], results, hiddens, inhGArray]
            else:
                # update SNP number and SNP positions
                if args.snp_position:
                    snpPoss = list(np.array(snpPoss)[varList])
                results[0] = len(genoArray)

            ## ABANDONED. Check if the MAF of the variants are lower than a specified threshold in either group
            #lowMAF, genoArray = ws.genoArrayMAF(genoArray, groupIDs, threshold=mafThreshold)
            #if lowMAF:
            #    skipLoop(results, resultHeader)
            #    return 0

            if args.permutation:
                stats = []      # a list of the function names (string) for permutation test
                scores = []     # a list of the divergence scores for permutation test

            indicator = 1       # The indicator for each element in results, starting from 1 the first statistics (0 is snpNum)

            ## If output mean allele frequency difference
            if args.allele_freq_diff:
                afd = ws.meanAlleleFreqDiffAllel(genoArray, groupIDs)
                results[indicator] = afd
                indicator += 1
                if args.permutation:
                    stats.append('afd')
                    scores.append(afd)

            ## If output FST
            if args.fst or args.var_comp:
                fst, varcomp, totalvar = ws.wcFst(genoArray, groupIDs)
                if args.fst:
                    results[indicator] = fst
                    indicator += 1
                    if args.permutation:
                        stats.append('fst')
                        scores.append(fst)
                if args.var_comp:
                    results[indicator] = varcomp
                    indicator += 1
                    if args.permutation:
                        stats.append('varcomp')
                        scores.append(varcomp)
                if args.inherent:
                    hiddens['varcomp'] = varcomp
                    hiddens['totalvar'] = totalvar

            ## If output CSS
            if args.css:
                c = ws.css(genoArray, groupIDs)
                results[indicator] = c
                indicator += 1
                if args.permutation:
                    stats.append('css')
                    scores.append(c)

            ## If output CSS-modified
            if args.css_mod:
                cm = ws.cssMod(genoArray, groupIDs)
                results[indicator] = cm
                indicator += 1
                if args.permutation:
                    stats.append('cssmod')
                    scores.append(cm)

            ## If output dXY
            if args.dxy:
                dxy = ws.dxy(genoArray, groupIDs, wSize)
                results[indicator] = dxy
                indicator += 1
                if args.permutation:
                    stats.append('dxy')
                    scores.append(dxy)

            ## If doing permuation test
            if args.permutation:
                if args.inherent:
                    inhGArray = genoArray
                    # the permutation is allocated to another function in window_merge.py after all results are available
                else:
                    ps = ws.permuP(scores, stats, args.max_permutation, args.max_extreme, args.min_extreme, genoArray=genoArray, groupIDs=groupIDs, wSize=wSize)
                    for p in ps:
                        results[indicator] = p
                        indicator += 1

            return [snpPoss, results, hiddens, inhGArray]


## Write the results into output files
def writeResults(results):
    if args.snp_position:
        if len(results[0]) != 0:
            writeSP.write('\t'.join([str(x) for x in results[0]]) + '\n')
        else:
            writeSP.write('NA\n')
    
    writeResult.write('\t'.join([str(x) for x in results[1]]) + '\n')
    return 0


## Write the start and end position of the window
def writeWindowStartEnd(start, wSize):
# this function assumes the given start position is 1-based
    writeWSE.write(str(start) + '\t' + str(start + wSize - 1) + '\n')
    return 0



#### Main loop for sliding window
## Set timer
sTime = time.time()

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


## Compose the list for fetching the sliding windows, note that these values are converted to 0-based, half-opened coordinates
# fetchCoord = [(s1, e1), (s2, e2), ... ]
fetchCoord = []
i = posStart - 1
while i + wSize <= posEnd:
    fetchCoord.append((i, i + wSize))
    i += gap


## Define storing values
totalWindow = 0     # store the total number of window scanned

## Loop through SNPs in the VCF
print('Time spent: ' + str(round(time.time() - sTime, 2)) + ' sec; window sliding began.')
# If using inherent windows
if args.inherent:
    openResults = []    # the return of windowStats for each openned sub-window
    for S, E in fetchCoord:
        i = S           # the i is the flag that tracts the starting position of each sub-window
        if S == posStart - 1:
        # if this is the starting window
            # Construct the non-overlaping sub-windows until it reaches the end of the current window
            openWindowSE = []   # the (s, e) positions of the sub-windows that are openned and reserved
            while i + gap <= E:
                openWindowSE.append((i, i + gap))
                i += gap

            for s, e in openWindowSE:
                windowVcf = vcf.fetch(contig=chrom, start=s, end=e)
                results = windowStats(windowVcf, groupIDs, wSize=gap, resultHeader=resultHeader)
                openResults.append(results)

        else:
        # if this is not the starting window
            openResults.pop(0)      # remove the first element as it's not in the current window
            windowVcf = vcf.fetch(contig=chrom, start=lastEndPos, end=E)
            results = windowStats(windowVcf, groupIDs, wSize=gap, resultHeader=resultHeader)
            openResults.append(results)

        lastEndPos = E
        # the ending point of current window will serve as the starting point of the new sub-window in the next window

        mResults = wm.mergeResults(openResults, resultHeader)

        # Permutation
        if args.permutation:
            permuPs = wm.permuPMerge(openResults, resultHeader, mResults, args.max_permutation, args.max_extreme, minExtm=args.min_extreme, groupIDs=groupIDs, wSize=wSize)
            mResults += permuPs

        writeResults(mResults)
        writeWindowStartEnd(S + 1, wSize)

        # output the current number of window processed
        totalWindow += 1
        if totalWindow % windowReportSize == 0:
            print(str(totalWindow) + ' windows scanned; time spent: ' + str(round(time.time() - sTime, 2)) + ' sec.')


else:
    for s, e in fetchCoord:
        # fetch window vcf   
        windowVcf = vcf.fetch(contig=chrom, start=s, end=e)
        # parse the record and calculate statistics
        results = windowStats(windowVcf, groupIDs, wSize=wSize, resultHeader=resultHeader)
        writeResults(results)
        writeWindowStartEnd(s + 1, wSize)

        # output the current number of window processed
        totalWindow += 1
        if totalWindow % windowReportSize == 0:
            print(str(totalWindow) + ' windows scanned; time spent: ' + str(round(time.time() - sTime, 2)) + ' sec.')


## Close openned files
writeWSE.close()
writeResult.close()
if args.snp_position:
    writeSP.close()


#### Final output
print('Finished; total window scanned: ' + str(totalWindow) + '; total time spent: ' + str(round(time.time() - sTime, 2)) + ' sec.')
