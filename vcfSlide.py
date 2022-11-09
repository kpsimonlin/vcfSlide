## This program runs sliding windows on a given VCF file and applies statistics on them.
## Author: Kung-Ping Lin; last editted: 20221108.

import argparse
import vcf
import gzip

## Argument parsing
parser = argparse.ArgumentParser(
                description='Iterate through user-defined sliding windows on a VCF file and apply statistics on them.')

parser.add_argument('-v', '--vcf', help='The path of the input VCF file.', required=True)
parser.add_argument('-w', '--windowSize', help='The size of the window.', required=True, type=int)
parser.add_argument('-g', '--gap', help='The gap for every window slide.', required=True, type=int)
args = parser.parse_args()

vcfPath = args.vcf
wSize = args.windowSize
gap = args.gap


## Main loop for sliding window
vcfReader = vcf.Reader(filename=vcfPath)

for record in vcfReader:
    
