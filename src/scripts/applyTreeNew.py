"""Make decision tree for snv and indel cutoffs.
   limitFeatures means do not wxs depth
   limit data places cutoffs on len of variant and wxs depth
"""
from sklearn import preprocessing
import numpy
from collections import defaultdict
import csv, sys, pickle, pandas, argparse
from sklearn import tree, metrics
from sklearn.ensemble import ExtraTreesClassifier
from sklearn.externals.six import StringIO
from sklearn.cross_validation import train_test_split, StratifiedShuffleSplit
#from sklearn.metrics import precision_recall_curve

import mkTree

def loadData(mat, chrom, sex):
    df = pandas.read_csv(mat, delimiter='\t')

    # recode str cols
    df.loc[:, 'qualVal'] = df.apply(lambda row: mkTree.mkQualStr(row), axis=1)
    df.loc[:, 'mvarZygosityWrongVal'] = df.apply(lambda row: mkTree.mkStr(row, 'mvarZygosityWrong'), axis=1)
    if chrom in ('Y', 'X') and sex == 'M':
        df.loc[:, 'hasNoCallVal'] = 0
    else:
        df.loc[:, 'hasNoCallVal'] = df.apply(lambda row: mkTree.mkStr(row, 'hasNoCall'), axis=1)

    df.loc[:, 'varLen'] = df.apply(lambda row: mkTree.mkVarLen(row), axis=1)
    df.loc[:, 'tumorFrac'] = df.apply(lambda row: mkTree.mkTumorFrac(row), axis=1)
    df.loc[:, 'tumorTotFrac'] = df.apply(lambda row: mkTree.mkTumorTotFrac(row), axis=1)
    return df

def main(args):
    snvFeatures = ['varDepth', 'refDepth', 'totDepth', 'allDepth', 'varFrac',
                   'goodDepthFrac', 'hasNoCallVal', 'mvarZygosityWrongVal',
                   'tumorAlt', 'tumorRef', 'tumorEffDepth',
                   'tumorTotDepth', 'qualVal', 'tumorFrac',
                   'tumorTotFrac', 'noCallSampleCount']
    indelFeatures = snvFeatures + ['varLen']
    dfInit = loadData(args.mat, args.chrom, args.sex)

    indel_idx = dfInit.varLen != 0
    snvDatInit = dfInit[~indel_idx]
    indelDatInit = dfInit[indel_idx]

    snvDat = snvDatInit[snvFeatures]
    indelDat = indelDatInit[indelFeatures]

    treeSnv = pickle.load(open(args.snvTreePickle, 'rb'))
    treeIndel = pickle.load(open(args.indelTreePickle, 'rb'))
    resultSnv = treeSnv.predict_proba(snvDat)
    resultIndel = treeIndel.predict_proba(indelDat)

    dfInit.loc[~indel_idx, 'Prob'] = [x[1] for x in resultSnv]
    dfInit.loc[indel_idx, 'Prob'] = [x[1] for x in resultIndel]
    dfInit.to_csv(args.outFile, sep='\t', index=False)

if __name__ == "__main__":
    desc = 'Apply decision tree for indels and snvs.'
    parser = argparse.ArgumentParser(description=desc)
    argLs = ('chrom', 'sex', 'mat',
             'snvTreePickle', 'indelTreePickle',
             'outFile')
    for param in argLs:
        parser.add_argument(param)
    args = parser.parse_args()
    main(args)
