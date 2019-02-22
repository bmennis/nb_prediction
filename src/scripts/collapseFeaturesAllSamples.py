"""Mk one feature matrix for sample and chrom.
   Use for all samples (not training).
"""
###Adjusting const import was import sfConst, but trying just include: const.py
include: "const.py"
import argparse, csv, sys, os, pandas
from functools import reduce
PWD = os.getcwd().split('code')[0]
sys.path.append(PWD + 'code/rules/')


# def loadData(dataFile):
#     """Load {} of chrom:pos:ref:alt to fields."""
    
def main(args):
    argDict = vars(args)
    dfs = {x:pandas.read_csv(argDict[x], dtype={'chrom': object}, delimiter='\t')
           for x in argDict if not x == 'outFile'}
    featureDfs = {x:dfs[x] for x in dfs if not x in ('wxsDepth', 'noCallCounts')}
    for dfName in featureDfs:
        if dfName != 'cgiFeats':
            featureDfs[dfName].rename(index=str, columns={'match':dfName + '_match'}, inplace=True)

    colKeyLs = ['chrom', 'pos', 'ref', 'alt']
    resultDf = reduce( lambda left, right: pandas.merge(left, right, on=colKeyLs), featureDfs.values())
    
    colKeyLs = ['chrom', 'pos']
    finalDf = pandas.merge(resultDf, dfs['noCallCounts'], on=colKeyLs, how='left')
    finalDf['noCallSampleCount'].fillna(0, inplace=True)

    finalDf.to_csv(args.outFile, sep='\t', index=False)

if __name__ == "__main__":
    desc = 'Make on feature matrix for sample and chrom.'
    parser = argparse.ArgumentParser(description=desc)
    argLs = ['noCallCounts', 'cgiFeats', 
             'outFile']

    for param in argLs:
        parser.add_argument(param)
    args = parser.parse_args()
    main(args)
