## Module for methods combining statistics from several non-overlapping windows
## Kung-Ping Lin, 20231115.
## These methods take a list of "results" for input, the "results" are a two-element nested list structured for a single sliding window like this:
#### [[snp_pos1, snp_pos2 ... ], [snp_number, stats1, stats2, ... ]] where stats are reported statistics like CSS

import numpy as np

## Main function, of combining results from several sliding windows
def mergeResults(resultList, resultHeader):
# resultHeader is a list recording the statistics specified, starting with snp_number
    mergedResults = [[], []]                # the merged results list
    statsNum = len(resultList[0][1]) - 1    # the number of statistics ready to be combined w/o snp_number
    statsHeader = resultHeader[1:]          # statistics to combine

    ## Merge the number of SNPs and the statistics
    mergedNum = np.sum(r[1][0] for r in resultList)     # total number of SNPs in the combined window
    if mergedNum == 0:
    # if there is no SNPs in the combined window, return directly
        mergedResults[1].append(mergedNum)
        mergedResults[1] + ['NA']*statsNum
        return mergedResults

    if statsNum > 0:
        statsList = []
        for i in range(statsNum):
            statsList.append([])   # the list for storing each stats from each sub-window, NAs would not be stored
        snpNumList = []             # the list for storing SNP number from each sub-window
    ## Main loop for combining statistics and concatenate SNP positions
    for results in resultList:
        # concatenate SNP positions
        mergedResults[0] += results[0]

        if statsNum > 0:
        # if there are statsitics to combine, construct statsList and snpNumeList
            for i, s in enumerate(results[1][1:]):
                statsList[i].append(s)
                # note that NAs could present in this list, will deal with this issue when specific merging method is specified later
            snpNumList.append(results[1][0])

    # Constructing the second sublist of mergedResults
    mergedResults[1].append(mergedNum)
    for i, h in enumerate(statsHeader):
    # loop through each statistics, determine which merging method to use
        if h == 'Allele_freq_diff':
            afd = weightedAverage(statsList[i], snpNumList)
            mergedResults[1].append(afd)
        
        elif h == 'dXY':
            dxy = simpleAverage(statsList[i])
            mergedResults[1].append(dxy)
    
        elif h == 'CSS_mod':
            css = simpleSum(statsList[i])
            mergedResults[1].append(css)

    return mergedResults


##### Methods for merging statistics

## Average statistics weighted by SNP number
def weightedAverage(statsList, snpNumList):
    notNA = np.array(statsList) != 'NA'     # selection array indicating the element that is not NA
    fStatsList = np.array(np.array(statsList)[notNA], dtype=float)
    fSnpNumList = np.array(snpNumList)[notNA]

    weights = fSnpNumList / np.sum(fSnpNumList)
    weightedStats = fStatsList * weights

    return np.sum(weightedStats)


## Average statistics (simple mean) over all sub-windows, don't remove windows when no SNPs were present
def simpleAverage(statsList):
    nanStatsList = [np.nan if s=='NA' else s for s in statsList]
    return np.nansum(nanStatsList) / len(statsList)

## Sum up statistics over all sub-windows, ignore the windows with no SNPs
def simpleSum(statsList):
    nanStatsList = [np.nan if s=='NA' else s for s in statsList]
    return np.nansum(nanStatsList)



