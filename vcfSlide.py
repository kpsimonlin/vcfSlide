## This program runs sliding windows on a given VCF file and applies statistics on them.
## Author: Kung-Ping Lin; last editted: 20221108.

import argparse
import vcf

#### Argument parsing
parser = argparse.ArgumentParser(
                description='Iterate through user-defined sliding windows on a VCF file and apply statistics on them.')

parser.add_argument('-v', '--vcf', help='The path of the input VCF file', required=True)
parser.add_argument('-w', '--window-size', help='The size of the window', required=True, type=int)
parser.add_argument('-g', '--gap', help='The gap for every window slide', required=True, type=int)
parser.add_argument('-o', '--output-prefix', help='The path and prefix of output files', required=True)
parser.add_argument('--snp-position', help='Write the SNP positions in every window', action='store_true', default=False)
args = parser.parse_args()

vcfPath = args.vcf
wSize = args.window_size
gap = args.gap
outPrefix = args.output_prefix

## Print argument inputs
print('Arguments:')
for arg in vars(args):
    print(arg, getattr(args, arg))
print()

#### Function for writing output files for statistics
## Open output files
if args.snp_position:
    writeSP = open(outPrefix + '.snppos.tsv', 'w')

## Write the statistics for the records in the window
def writeRecords(records):
    ## If output SNP positions in each window
    if args.snp_position:
        poss = [r.POS for r in records]
        writeSP.write('\t'.join([str(x) for x in poss]) + '\n')
    
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

## Loop through SNPs in the VCF
for record in vcfReader:
    pos = record.POS
    # add the first position into starts
    if starts == []:
        starts.append(pos)
        print('New window: ' + str(pos) + '-' + str(pos + wSize - 1))

    openRecords.append(record)

    ## Determine whether new windows need to be opened
    # if the current SNP position exceeds where the next window should be, open new window
    if pos >= starts[-1] + gap:
        # set new end to current right border before looping
        newStart = starts[-1] + gap
        # add new window starts until all windows start before the SNP are added
        while newStart <= pos:
            print('New window: ' + str(newStart) + '-' + str(newStart + wSize - 1))
            starts.append(newStart)
            newStart += gap

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
            writeRecords(outRecords)
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
    writeRecords(outRecords)
    writeWindowStartEnd(start, wSize)
