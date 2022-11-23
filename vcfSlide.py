## This program runs sliding windows on a given VCF file and applies statistics on them.
## Author: Kung-Ping Lin; last editted: 20221108.

import argparse
import vcf
import window_stats as ws

#### Argument parsing
parser = argparse.ArgumentParser(
                description='Iterate through user-defined sliding windows on a VCF file and apply statistics on them.')

parser.add_argument('-v', '--vcf', help='The path of the input VCF file', required=True)
parser.add_argument('-w', '--window-size', help='The size of the window', required=True, type=int)
parser.add_argument('-g', '--gap', help='The gap for every window slide', required=True, type=int)
parser.add_argument('-o', '--output-prefix', help='The path and prefix of output files', required=True)
parser.add_argument('--snp-position', help='Output the positions of SNPs in every window', action='store_true', default=False)
parser.add_argument('--snp-number', help='Output the number of SNPs in every window', action='store_true', default=False)

parser.add_argument('-G', '--group', help='The path to two files defining two groups of individual to compared, required for divergence statistics', nargs=2, default=0)
parser.add_argument('--allele-freq-diff', help='Output the mean ALT allele frequency difference between groups in each window', action='store_true', default=False)
parser.add_argument('--fst', help='Output Weir-Cockerham FST score between groups in each window', action='store_true', default=False)
parser.add_argument('--var-comp', help='Output the numerator of the  Weir-Cockerham FST score between groups in each window', action='store_true', default=False)
parser.add_argument('--css', help='Output the Cluster Separation Score (CSS) by Jones et al. (2012) in each window', action='store_true', default=False)

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

vcfPath = args.vcf
wSize = args.window_size
gap = args.gap
outPrefix = args.output_prefix

## Print argument inputs
print('Arguments:')
for arg in vars(args):
    print(arg, getattr(args, arg))
print('')


#### Hidden parameters
windowReportSize = 100      # the denominator of reporting window scanned
missThreshold = 0.9         # the threshold for missing rate in either group


#### Function for writing output files for statistics
## Open output files
if args.snp_position:
    writeSP = open(outPrefix + '.snppos.tsv', 'w')
if args.snp_number:
    writeSN = open(outPrefix + '.snpnum.tsv', 'w')
if args.allele_freq_diff:
    writeAFD = open(outPrefix + '.afdiff.tsv', 'w')
if args.fst:
    writeFST = open(outPrefix + '.fst.tsv', 'w')
if args.var_comp:
    writeVCP = open(outPrefix + '.varcom.tsv', 'w')
if args.css:
    writeCSS = open(outPrefix + '.css.tsv', 'w')


## Write the statistics for the records in the window
def writeRecords(records, groupIDs):
# groupIDs was called when the first record encountered; it defines the indices of two groups in the VCF samples

    ## If output SNP positions
    if args.snp_position:
        poss = ws.snpPositions(records)
        writeSP.write('\t'.join([str(x) for x in poss]) + '\n')
    
    ## if output snp number
    if args.snp_number:
        num = ws.snpNumber(records)
        writeSN.write(str(num) + '\n')

    ## If has group information as input
    noRecord = False    # a flag if there's no record in this window or no variants after subsetting

    if groupIDs != 0:
        # if no records present in this window
        if len(records) == 0:
            noRecord = True

        ## Convert records to genotypeArray object and
        # subset the genoArray to contain only the two groups interested
        subGenoArray, subGroupIDs = ws.subsetGenotypeArray(ws.toGenotypeArray(records, noRecord=noRecord), groupIDs, noRecord=noRecord)

        ## Check if only invariants remain, if not, remove the invariant sites
        ifAllInvar, varSubGenoArray = ws.genoArrayInvar(subGenoArray, remove=True, noRecord=noRecord)
        if ifAllInvar:
            noRecord = True

        ## Check if the variants are missed for a certain threshold in either group
        isMiss, mfVarSubGenoArray = ws.genoArrayMiss(varSubGenoArray, subGroupIDs, threshold=missThreshold, noRecord=noRecord)
        if isMiss:
            noRecord = True

        ## If output mean allele frequency difference
        if args.allele_freq_diff:
            afd = ws.meanAlleleFreqDiffAllel(mfVarSubGenoArray, subGroupIDs, noRecord=noRecord)
            writeAFD.write(str(afd) + '\n')

        ## If output FST
        if args.fst:
            if args.var_comp:
                fst, varcomp = ws.wcFst(mfVarSubGenoArray, subGroupIDs, reportA=True, noRecord=noRecord)
                writeFST.write(str(fst) + '\n')
                writeVCP.write(str(varcomp) + '\n')

            else:
                fst = ws.wcFst(mfVarSubGenoArray, subGroupIDs, noRecord=noRecord)[0]
                writeFST.write(str(fst) + '\n')

        ## If output CSS
        if args.css:
            css = ws.css(mfVarSubGenoArray, subGroupIDs, noRecord=noRecord)
            writeCSS.write(str(css) + '\n')

    return 0


## Write the start and end position of the window
writeWSE = open(outPrefix + '.winpos.tsv', 'w')
def writeWindowStartEnd(start, wSize):
    writeWSE.write(str(start) + '\t' + str(start + wSize - 1) + '\n')
    return 0


#### Main loop for sliding window
# PyVCF reader
vcfReader = vcf.Reader(filename=vcfPath)

## Define storing lists
starts = []    # define the start position of the current opened windows; 1-based
openRecords = []    # store the openning SNP objects
totalWindow = 0     # store the total number of window scanned

## Loop through SNPs in the VCF
for record in vcfReader:
    pos = record.POS
    # add the first position into starts, also, load the indices of group ID if applicable
    if starts == []:
        starts.append(pos)
        totalWindow += 1

        # get the group indices if --group is defined, this is for statistics require group definition
        if args.group != 0:
            groupIDs = ws.getGroupID(record, args.group)
        else:
            groupIDs = 0

    openRecords.append(record)

    ## Determine whether new windows need to be opened
    # if the current SNP position exceeds where the next window should be, open new window
    if pos >= starts[-1] + gap:
        # set new end to current right border before looping
        newStart = starts[-1] + gap
        # add new window starts until all windows start before the SNP are added
        while newStart <= pos:
            starts.append(newStart)
            newStart += gap
            totalWindow += 1
            if totalWindow % windowReportSize == 0:
                print('Total window scanned: ' + str(totalWindow))

    ## Determine whether some old windows need to be closed and outputed
    firstRight = starts[0] + wSize - 1    # the end of the leftest openning window
    # if the current SNP position exceeds the first openning window, delete and output window
    if pos > firstRight:
        # set old end to the start of the leftest openning window
        oldEnd = firstRight
        # pull out all positions from all openned SNPs
        poss = [r.POS for r in openRecords]
        # set list for the indices of the SNPs to be deleted in openRecords
        indiTBD = []

        while oldEnd < pos:
            
            ## Extract the records in the window
            start = oldEnd - wSize + 1
            outRecords = [] # records in the window
            for i, p in enumerate(poss):
                if p >= start and p <= oldEnd:
                    outRecords.append(openRecords[i])
                    # if the current position only present in the leftest window
                    if p < start + gap:
                        indiTBD.append(i)

            ## Output statistics from the records in the window
            writeRecords(outRecords, groupIDs)
            writeWindowStartEnd(start, wSize)

            ## Remove the start of the leftest window from starts
            starts.pop(0)

            # new leftest window end before next loop
            oldEnd += gap

        ## Delete the SNPs that only present in the leftest openning window
        if len(indiTBD) != 0:
            del openRecords[indiTBD[0] : indiTBD[-1] + 1]


#### Output the remaining sliding windows
for start in starts:
    end = start + wSize - 1
    outRecords = []
    for r in openRecords:
        if r.POS >= start and r.POS <= end:
            outRecords.append(r)
    writeRecords(outRecords, groupIDs)
    writeWindowStartEnd(start, wSize)


#### Close openned files
if args.snp_position:
    writeSP.close()
if args.snp_number:
    writeSN.close()
if args.allele_freq_diff:
    writeAFD.close()
if args.fst:
    writeFST.close()
if args.var_comp:
    writeVCP.close()
if args.css:
    writeCSS.close()


#### Final output
print('Finished; total window scanned: ' + str(totalWindow))
