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
import pydotplus

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

def calcY(row):
    evalLs = ('D_match', 'TD_match', 'platypus_match',
              'trevorGatk_match', 'SI_match', 'nciGatk_match',
              'bam2mpg_match')
    return max([row[x] for x in evalLs])

def mkQualStr(row):
    if row['qual'] == 'PASS':
        return 1
    if row['qual'] == 'VQLOW':
        return 0
    print(row['qual'])
    i = 1/0

def mkStr(row, col):
    if row[col] == 'T':
        return 1
    if row[col] == 'F':
        return 0
    print(row[col])
    i = 1/0

def mkVarLen(row):
    return len(row['ref']) - len(row['alt'])

def mkTumorFrac(row):
    td = row['tumorTotDepth']
    if td:
        return float(row['tumorAlt']) / float(td)
    return 0

def mkTumorTotFrac(row):
    td = row['tumorTotDepth']
    if td:
        return float(row['tumorEffDepth']) / float(td)
    return 0

def loadData(mat, limitData, useChroms):
    """useChroms:
         M: only chrM
         male: X and Y from males
         other: all other
    """
    dfInit = pandas.read_csv(mat, delimiter='\t', dtype={'chrom':object})

    # select variant subset by chromosomes
    # M not treated differently b/c no M indel calls for exome
    # if useChroms == 'M':
    #     df = dfInit[dfInit.chrom=='M']

    if useChroms == 'male':
        # did not have enough data for male alone
        df = dfInit.query("(chrom=='Y' | chrom == 'X') & sex=='M'")
    elif useChroms == 'other':
        # exclude male X and Y b/c they have lower noCall sample counts
        # and they have true no call status b/c there is only one chrom
        df = dfInit.query("not ( (chrom=='Y' | chrom == 'X') & sex=='M')")
    else:
        i = 1/0

    # select one variant per position
    byPosRandom = df.groupby(['chrom', 'pos',]).agg(lambda x: x.iloc[numpy.random.randint(0, len(x))])
    byPosRandom['Y'] = byPosRandom.apply(lambda row: calcY(row), axis=1)
    # recode str cols
    byPosRandom['qualVal'] = byPosRandom.apply(lambda row: mkQualStr(row), axis=1)
    byPosRandom['mvarZygosityWrongVal'] = byPosRandom.apply(lambda row: mkStr(row, 'mvarZygosityWrong'), axis=1)
    byPosRandom['hasNoCallVal'] = byPosRandom.apply(lambda row: mkStr(row, 'hasNoCall'), axis=1)
    byPosRandom['varLen'] = byPosRandom.apply(lambda row: mkVarLen(row), axis=1)
    byPosRandom['tumorFrac'] = byPosRandom.apply(lambda row: mkTumorFrac(row), axis=1)
    byPosRandom['tumorTotFrac'] = byPosRandom.apply(lambda row: mkTumorTotFrac(row), axis=1)
    if limitData:
        limitWxsCovDf = byPosRandom[byPosRandom.wxsDepth > 10]
        limitVarLen = limitWxsCovDf[(limitWxsCovDf.varLen > -13) & (limitWxsCovDf.varLen < 20)]
        return limitVarLen
    return byPosRandom

def mkTree(dat, features, treePickle, treePdf):
#    print('debug', treePdf)
    clf = tree.DecisionTreeClassifier(max_depth=5)
    y = dat['Y']
    x = dat[features]
    clf = clf.fit(x, y)
    dot_data = StringIO()
    tree.export_graphviz(clf, feature_names=features, out_file=dot_data)
    graph = pydotplus.graph_from_dot_data( dot_data.getvalue() )
#    print('debug', treePdf)
    graph.write_pdf(treePdf)
    pickle.dump(clf, open(treePickle, 'wb'))

def findFeaturePerformance(dat, features, treeImportanceOut):
    y = dat['Y']
    x = dat[features]
    forest = ExtraTreesClassifier(n_estimators=300,
                                  random_state=13,
                                  bootstrap=True,
                                  max_features=len(features)-2,
                                  min_samples_split=2,
                                  max_depth=len(features)-1,
                                  min_samples_leaf=13,
                                  n_jobs=4)
    forest.fit(x, y)
    importances = forest.feature_importances_
    std = numpy.std([atree.feature_importances_ for atree in forest.estimators_],
                 axis=0)
    indices = numpy.argsort(importances)[::-1]
    
    # Print the feature ranking
    with open(treeImportanceOut, 'w') as fout:
        print("Feature ranking:", file=fout)
        for f in range(x.shape[1]):
            ls = (features[indices[f]],
                  f + 1, indices[f],
                  importances[indices[f]])
            print("%s, %d. feature %d (%f)" % ls, file=fout)

def crossVal(dat, features, treePerformanceOut, rocPng, roc_data_file):
    y = dat['Y']
    x = dat[features]
    sss = StratifiedShuffleSplit(y, 10, test_size=0.1, random_state=442)
    clf = tree.DecisionTreeClassifier(max_depth=5)
    plt.figure()
    with open(treePerformanceOut, 'w') as fout:
        for cVal, crossDat in enumerate(sss):
            train_index, test_index = crossDat
            X_train, X_test = x.iloc[train_index], x.iloc[test_index]
            y_train, y_test = y[train_index], y[test_index]
            clf = clf.fit(X_train, y_train)
            print('all data', file=fout)
            print(metrics.classification_report(y_test, clf.predict(X_test)), file=fout)
            probs = clf.predict_proba(X_test)[:,1]
            pre, rec, thresh_pr = metrics.precision_recall_curve(y_test, probs,
                                                                 pos_label=1)
            plt.plot(rec, pre, label='All_%d' % (cVal,))

            # find indels and snvs
            X_test_indel_idx = X_test.varLen != 0
            X_test_snv_idx = ~X_test_indel_idx

            x_test_indel, x_test_snv = X_test[X_test_indel_idx], X_test[X_test_snv_idx]
            y_test_indel, y_test_snv = y_test[X_test_indel_idx], y_test[X_test_snv_idx]

            # eval only test indels
            print('only indels', file=fout)
            print(metrics.classification_report(y_test_indel, clf.predict(x_test_indel)), file=fout)
            probs = clf.predict_proba(x_test_indel)[:,1]
            #print(probs)
            pre, rec, thresh_pr = metrics.precision_recall_curve(y_test_indel, probs,
                                                                 pos_label=1)
            plt.plot(rec, pre, label='Indel_%d' % (cVal,))

            # find snvs
            print('only snvs', file=fout)
            print(metrics.classification_report(y_test_snv, clf.predict(x_test_snv)), file=fout)
            probs = clf.predict_proba(x_test_snv)[:,1]
            pre, rec, thresh_pr = metrics.precision_recall_curve(y_test_snv, probs,
                                                                 pos_label=1)
            plt.plot(rec, pre, label='Snv_%d' % (cVal,))

    plt.xlim([0.0, 1.05])
    plt.ylim([0.0, 1.05])
    plt.xlabel('Recall')
    plt.ylabel('Precision')
    plt.title('Precision-Recall Curve for Decision Tree')
    plt.legend(loc="lower left")
    plt.savefig(rocPng)

def main(args):
    baseFeatures = ['varDepth', 'refDepth', 'totDepth', 'allDepth', 'varFrac',
                    'goodDepthFrac', 'hasNoCallVal', 'mvarZygosityWrongVal', 'tumorAlt', 'tumorRef', 'tumorEffDepth',
                    'tumorTotDepth', 'qualVal', 'tumorFrac', 'tumorTotFrac', 'varLen', 'noCallSampleCount']
    if args.limitFeatures == 'no':
        features = baseFeatures + ['wxsDepth'] #, 'varLen']
    elif args.limitFeatures == 'yes':
        features = baseFeatures
    else:
        i = 1/0

    if args.limitData == 'no':
        limitData = False
    if args.limitData == 'yes':
        limitData = True

    dat = loadData(args.mat, limitData, args.chroms)
    findFeaturePerformance(dat, features, args.treeImportance)
    crossVal(dat, features, args.treePerformance, args.rocPng)
    mkTree(dat, features, args.treePickle, args.treePdf)

if __name__ == "__main__":
    desc = 'Decision tree for indels and snvs.'
    parser = argparse.ArgumentParser(description=desc)
    argLs = ('limitFeatures', 'limitData', 'chroms', 'mat', 'treePickle', 'treePdf',
             'treePerformance', 'treeImportance', 'rocPng')

    for param in argLs:
        parser.add_argument(param)
    args = parser.parse_args()
    main(args)
