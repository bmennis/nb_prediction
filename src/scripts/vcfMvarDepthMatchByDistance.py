"""Match vcf ref and alt by edit distances to mvar ref and alt."""
import operator
import editdistance

def getEffFrac(mvarKey, mvar):
    # ignore ? b/c they are blank
    eff = sum([int(mvar[x]) for x in ('allele1ReadCount', 'allele2ReadCount') if mvar[x]])
    t = int(mvar['totalReadCount'])
    if eff > t:
        i = 1/0
    return float(eff)/t

def getTotDepth(mvarKey, mvar):
    return int(mvar['totalReadCount'])

def getLenDiff(mvarKeyMvar, alt):
    mvarKey, mvar = mvarKeyMvar
    return abs(len(mvar[mvarKey]) - len(alt))

def breakTies(bestLs, alt):
    """See 1:211260977 ACACG to A for an example of a real tie.
       I use the highest effective covergage fraction
       to break the tie.
    """
    if len(bestLs) == 1:
        return bestLs[0:1]

    ls = []
    for bb in bestLs:
        ls.append( (getEffFrac(*bb), bb) )
    ls.sort( key=operator.itemgetter(0) )
    ls.reverse()
#    print('shit', ls)
    bestEff = ls[0][0]
    bestLs2 = [(bb[0], bb[1]) for eff, bb
               in ls if eff == bestEff]
    if len(bestLs2) != 1:
        # choose by total depth
        ls = []
        for bb in bestLs2:
            ls.append( (getTotDepth(*bb), bb) )
        ls.sort( key=operator.itemgetter(0) )
        ls.reverse()
        bestEff = ls[0][0]
        bestLs3 = [(bb[0], bb[1]) for eff, bb
                   in ls if eff == bestEff]
        if len(bestLs3) != 1:
            ls = []
            for bb in bestLs3:
                ls.append( (getLenDiff(bb, alt), bb) )
            ls.sort( key=operator.itemgetter(0) )
            bestEff = ls[0][0]
            bestLs4 = [(bb[0], bb[1]) for eff, bb
                       in ls if eff == bestEff]
            if len(bestLs4) != 1:
                haveIt = False
                if len(bestLs4) == 2:
                    mk1, md1 = bestLs4[0]
                    mk2, md2 = bestLs4[1]
                    if md1[mk1] == md2[mk2]:
                        r1 = int(md1[mk1.replace('Seq','ReadCount')])
                        r2 = int(md2[mk2.replace('Seq','ReadCount')])
                        if r1>r2:
                            return bestLs4[0:1]
                        return bestLs4[1:]
                    else:
                        if int(md1['totalReadCount-T1']) > int(md2['totalReadCount-T1']):
                            return bestLs4[0:1]
                        elif int(md1['totalReadCount-T1']) < int(md2['totalReadCount-T1']):
                            return bestLs4[1:]

                        r1 = int(md1[mk1.replace('Seq','ReadCount')])
                        r2 = int(md2[mk2.replace('Seq','ReadCount')])
                        if r1>r2:
                            return bestLs4[0:1]
                        return bestLs4[1:]

                if not haveIt:
#                    print(ls)
                    i = 1/0
            return bestLs4[0:1]
        return bestLs3[0:1]
    return bestLs2[0:1]

def calcDis(ref, alt, mvarDict):
    mvarKeys = ('reference', 'allele1Seq', 'allele2Seq')
    mvarRef, mvarA1, mvarA2 = [mvarDict[x] for x in mvarKeys]
    if mvarA1 == mvarA2:
        d = editdistance.eval(ref, mvarRef) + editdistance.eval(alt, mvarA1)
        return [ (d, 'allele1Seq', mvarDict) ]

    ret = []
    for mvarKey in mvarKeys[1:]:
        if (mvarDict['reference'] != mvarDict[mvarKey] or (alt == mvarDict[mvarKey] and not '?' in mvarDict[mvarKey])) \
           and mvarDict[mvarKey.replace('Seq','ReadCount')]:
#            print('nuts', mvarKey, mvarDict[mvarKey.replace('Seq','ReadCount')])
            d = editdistance.eval(ref, mvarRef) + editdistance.eval(alt, mvarDict[mvarKey])
            ret.append( (d, mvarKey, mvarDict) )
    return ret

def findMatchesByDistance(ref, alt, mvarDictList):
    disLs = []
    for mvarDict in mvarDictList:
        disLs += calcDis(ref, alt, mvarDict)
    disLs.sort( key=operator.itemgetter(0) )
    # print('debug it')
    # for d, mvarKey, mvarDict in disLs:
    #     print('\\', d, mvarKey, mvarDict[mvarKey])
    bestDis = disLs[0][0]
    bestLs = [(mvarKey, mvarDict) for d, mvarKey, mvarDict
              in disLs if d == bestDis]
    if len(bestLs) == 0:
        i = 1/0
    if len(bestLs) == 1:
        return bestLs
    if len(bestLs) == 2:
        if bestLs[0] == bestLs[1]:
            return bestLs[0:1]
    # print( 'again!',  [(d, mvarKey, mvarDict) for d, mvarKey, mvarDict
    #                    in disLs if d == bestDis]) 
    #if len(bestLs) == 1:
    return breakTies(bestLs, alt)
    i = 1/0
