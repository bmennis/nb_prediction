import re
from collections import defaultdict
import mvarVcfFields, vcfMvarDepthDelFuncs, vcfMvarDepthInsFuncs, vcfMvarDepthSubsFuncs

def calcFracs(totDepth, altCount, effectiveDepth, pos):
    """Find alt frac for total reads.
       Find frac of both reads out of total.
    """
    totDepth = int(totDepth)
    altCount = int(altCount)
    # sometimes there are no reads
    if altCount == 0 and totDepth == 0:
        altFrac = 0
    else:
        altFrac = float(altCount)/totDepth
    if effectiveDepth == 0 and totDepth == 0:
        goodDepthFrac = 0
    else:
        goodDepthFrac = effectiveDepth/float(totDepth)

    if goodDepthFrac > float(1) or altFrac > float(1):
        # if goodDepthFrac > float(1) and pos == '9747513':
        #     # couldn't figure out why there are 45 tumor tot, but sum is 46
        #     goodDepthFrac = float(1)
        # else:
            print('bad calc', goodDepthFrac, altFrac, pos)
            i = 1/0

    return altFrac, goodDepthFrac

def evalResult(sp, totDepth, altCount, refDepth, effectiveDepth,
               hasNoCall, zygosityWrong,
               tumorTotDepth, tumorAltCount, tumorRefDepth, tumorEffectiveDepth,
               ):
    """Make allele fraction using total reads and
       fraction of reads counted towards allele1, 2, or ref.
    """
    altFracNormal, goodDepthFracNormal = calcFracs(totDepth, altCount, effectiveDepth, sp[1])
    altFracTumor, goodDepthFracTumor = calcFracs(tumorTotDepth, tumorAltCount, tumorEffectiveDepth, sp[1])

    # print( sp[0], sp[1], sp[3], sp[4], altCount, refDepth, str(effectiveDepth),
    #         str(totDepth), str(altFrac), str(goodDepthFrac),
    #         hasNoCall, zygosityWrong ) 

    return [sp[0], sp[1], sp[3], sp[4], altCount, refDepth, str(effectiveDepth),
            str(totDepth), str(altFracNormal), str(goodDepthFracNormal),
            hasNoCall, zygosityWrong, str(tumorAltCount), str(tumorRefDepth),
            str(tumorEffectiveDepth), str(tumorTotDepth)]

def guessVcfVar(ref, alt):
    vcfVar = 'snp'
    if len(alt) > len(ref):
        vcfVar = 'ins'
    elif len(alt) < len(ref):
        vcfVar = 'del'
    return vcfVar

def guessMvarType(mvarDict):
    if mvarDict['reference'] in (mvarDict['allele1Seq'], mvarDict['allele2Seq'] ):
        if mvarDict['reference'] == mvarDict['allele1Seq']:
            if len(mvarDict['reference']) < len(mvarDict['allele2Seq']):
                return 'ins', 'allele2Seq'
            elif len(mvarDict['reference']) > len(mvarDict['allele2Seq']):
                return 'del', 'allele2Seq'
            else:
                return 'snp', 'allele2Seq'

        if mvarDict['reference'] == mvarDict['allele2Seq']:
            if len(mvarDict['reference']) < len(mvarDict['allele1Seq']):
                return 'ins', 'allele1Seq'
            elif len(mvarDict['reference']) > len(mvarDict['allele1Seq']):
                return 'del', 'allele1Seq'
            else:
                return 'snp', 'allele1Seq'

    return '\\?', '\\?'

def guessMvarTypes(mvarDicts):
    """Loop over mvarDicts, and acc varType:(altAllele, mvarDict)
       One mvarDict can have two non-ref variants.
    """
    ret = defaultdict(list)
    for mvarDict in mvarDicts:
        if mvarDict['reference'] in (mvarDict['allele1Seq'], mvarDict['allele2Seq'] ):
            guessType, altAlleleKey = guessMvarType(mvarDict)
            ret[guessType].append( (altAlleleKey, mvarDict) )
        else: # two non-ref alleles
            for alleleSeq in ('allele1Seq', 'allele2Seq'):
                if len(mvarDict['reference']) > len(mvarDict[alleleSeq]):
                    ret['del'].append( (alleleSeq, mvarDict) )
                elif len(mvarDict['reference']) < len(mvarDict[alleleSeq]):
                    ret['ins'].append( (alleleSeq, mvarDict) )
                else:
                    ret['snp'].append( (alleleSeq, mvarDict) )
                    
    # # AGA to AG/nothing
    # print('shit', mvarDict['reference'][0:-1], mvarDict['allele1Seq'], mvarDict['allele2Seq'])
    # if mvarDict['reference'][0:-1] == mvarDict['allele1Seq'] and not mvarDict['allele2Seq'] \
    #         or mvarDict['reference'][0:-1] == mvarDict['allele2Seq'] and not mvarDict['allele1Seq']:
    #     return 'del'
                    
    if '\\?' in ret:
        print('cannot determine var type')
        print(mvarDicts)
        
    return ret

def evalMvarMatchForVcfPair(varType, mvarKey, mvarDat, ref, alt):
    # ins check needs more data from vcf T>TCAA vcf, CAA > blank/CAACAA
    otherAlleleKey = 'allele1Seq'
    if mvarKey == 'allele1Seq':
        otherAlleleKey = 'allele2Seq'
    print(ref, alt, mvarDat['reference'][2:], mvarDat[mvarKey][2:])
    v1 =  (varType == 'snp' and len(ref) > 1 and len(alt) > 1 and ref[1:] == mvarDat['reference']
            and alt[1:] == mvarDat[mvarKey] and alt[0] == ref[0]
            and not mvarDat[otherAlleleKey]) or \
           (varType == 'ins' and alt == mvarDat[mvarKey] and alt != mvarDat[otherAlleleKey]) or \
           (varType == 'ins' and ref[1:] == mvarDat['reference']
            and alt[1:] == mvarDat[mvarKey]
            and alt[1:] != mvarDat[otherAlleleKey] and alt != mvarDat[otherAlleleKey]) or \
           (varType == 'ins' and ref[1:] == mvarDat['reference']
            and alt[1:] == mvarDat[mvarKey]
            and alt[1:] == mvarDat[otherAlleleKey] and mvarDat['zygosity'] == 'hom') or \
           (varType != 'del' and ref == mvarDat['reference'] and alt == mvarDat[mvarKey]) or \
           (varType == 'del' and mvarDat['reference'] == ref[1:]
            and ref[0] == alt[0] and len(alt) == 1
            and not mvarDat[mvarKey] and mvarDat[otherAlleleKey]) or \
           (varType == 'del' and mvarDat['reference'] == ref[1:]
            and len(alt) == 1 and alt == ref[0] and 'hom' == mvarDat['zygosity']
            and not mvarDat[mvarKey] and mvarDat[otherAlleleKey]) or \
           (varType == 'del' and mvarDat['reference'] == ref[1:]
            and len(alt) == 1 and alt == ref[0] and 'hom' == mvarDat['zygosity']
            and not mvarDat[mvarKey] and not mvarDat[otherAlleleKey])
    v2 = False
    if len(mvarDat[mvarKey]) > 1 and mvarDat[otherAlleleKey] and len(mvarDat['reference']) > 1:
        v2 = (varType == 'snp' and mvarDat['reference'][0] == ref and \
              mvarDat[mvarKey][0] == alt and mvarDat[otherAlleleKey][0] != alt \
              and mvarDat['reference'][1:] == mvarDat[mvarKey][1:]) or \
             (varType == 'snp' and mvarDat['reference'][0] == ref and
              mvarDat[mvarKey][0] == alt and mvarDat[otherAlleleKey][0] != alt
              and ref == mvarDat[otherAlleleKey][0] and mvarDat['reference'][1] != mvarDat[otherAlleleKey][1]
              and len(mvarDat['reference']) == 2 and len(mvarDat[otherAlleleKey]) == 2 and 
              len(mvarDat[mvarKey]) == 2)
    v3 = False
    if not (v1 or v2):
        v3 = varType == 'snp' and ref[1:] == mvarDat['reference'] \
             and alt[1:] == mvarDat[mvarKey] and alt[1:] != mvarDat[otherAlleleKey] \
             and alt != mvarDat[otherAlleleKey] and alt[0] == ref[0]
    v4 = False
    if not (v1 or v2 or v3):
        print(ref, alt, mvarDat['reference'], mvarDat[mvarKey], mvarDat[otherAlleleKey])
        v4 = varType == 'del' and not mvarDat[mvarKey] and mvarDat[otherAlleleKey] and ref[1:] != mvarDat['reference']
        if len(alt) != 1 and varType == 'del':
            # this is a double deletion case
            print('yo', ref, alt)
            v4 |= mvarDat['reference'] == ref[1:] and mvarDat[mvarKey] == alt[1:] and not mvarDat[otherAlleleKey]

    return v1 or v2 or v3 or v4

def evalMvarMatch(useNewRef, varType, mvarKey, mvarDat, ref, alt, initRef, initAlt):
    if useNewRef:
        return evalMvarMatchForVcfPair(varType, mvarKey, mvarDat, ref, alt)
    return evalMvarMatchForVcfPair(varType, mvarKey, mvarDat, initRef, initAlt)

#    print('NUTS', varType, mvarDat['reference'], ref[1:], mvarDat[mvarKey], mvarDat[otherAlleleKey])

def handleSingleMvarGuessSingle(varType, mvarDatLs, ref, alt, initRef, initAlt):
    """mvaDat is (alleleSeqKey to mvarDict, mvarDict)"""
    print(mvarDatLs)
    mvarKey, mvarDat, useNewRef = mvarDatLs
    if not evalMvarMatch(useNewRef, varType, mvarKey, mvarDat, ref, alt, initRef, initAlt):
        print('YOOOOOO', varType, mvarKey, mvarDat, ref, alt, initRef, initAlt)
        i = 1/0
    return mvarKey, mvarDat, useNewRef

def mkNewMvarDict(mvarDict, ref, alt):
    """Make a new hom variant."""
    return {'allele1Seq':alt,
            'allele1ReadCount':str(int(mvarDict['allele1ReadCount']) + int(mvarDict['allele2ReadCount'])),
            'allele1ReadCount-T1':str(int(mvarDict['allele1ReadCount-T1']) + int(mvarDict['allele2ReadCount-T1'])),
            'allele2Seq':alt,
            'allele2ReadCount':str(int(mvarDict['allele1ReadCount']) + int(mvarDict['allele2ReadCount'])),
            'allele2ReadCount-T1':str(int(mvarDict['allele1ReadCount-T1']) + int(mvarDict['allele2ReadCount-T1'])),
            'reference':ref,
            'referenceAlleleReadCount':mvarDict['referenceAlleleReadCount'],
            'referenceAlleleReadCount-T1':mvarDict['referenceAlleleReadCount-T1'],
            'varType':'perryBlend',
            'zygosity':'hom',
            'totalReadCount':mvarDict['totalReadCount'],
            'totalReadCount-T1':mvarDict['totalReadCount-T1'],
            }

def bothMvarAllelesContributeToAlt(mvarDict, ref, alt, initRef, initAlt):
    """Ex AG > CA/CG in mvar file. Need counts for both alleles.
       Or C>T vcf and CGA/TGA/TG/.
    """
    #print('OUTCH', ref, alt, mvarDict['reference'][0], mvarDict['allele1Seq'][0], mvarDict['allele2Seq'][0] )
#    print(mvarDict['allele1Seq'], mvarDict['allele2Seq'])

    # refHits = len(re.findall(ref, mvarDict['reference']))
    # altHits1 = len(re.findall(alt, mvarDict['allele1Seq']))
    # altHits2 = len(re.findall(alt, mvarDict['allele2Seq']))
    # matchInside = False
    # if 1 == refHits and 1 == altHits1 and 1 == altHits2:
    #     refIdx = mvarDict['reference'].index(ref)
    #     altIdx1 = mvarDict['allele1Seq'].index(alt)
    #     altIdx2 = mvarDict['allele2Seq'].index(alt)
    #     if refIdx == altIdx1 and refIdx == altIdx2:
    #         matchInside = True        

    return mvarDict['allele1Seq'] != mvarDict['allele2Seq'] and mvarDict['allele1Seq'] and mvarDict['allele2Seq'] and \
           ( (len(mvarDict['reference']) == 2 and len(mvarDict['allele1Seq']) == 2 and len(mvarDict['allele2Seq']) == 2 and \
              len(ref) == 1 and len(alt) == 1 and mvarDict['reference'][0] == ref and \
              mvarDict['allele1Seq'][0] == alt and mvarDict['allele2Seq'][0] == alt)
             or \
             (ref == mvarDict['reference'][0]
              and alt == mvarDict['allele1Seq'][0] and alt == mvarDict['allele2Seq'][0]) )
             # or \
             #  matchInside )

def expandWhenMatchesBoth(mvarGuessTypes, varType, ref, alt, initRef, initAlt):
    """It is very complicated now. When a single mar line with 2 alt alleles
       that can match the vcf variant, collapse the counts into a new mvarDict.
       Can't use mvardict as blend and initial.
    """
    #print(ref, alt, varType)
    if varType != 'snp':
        return [], []
    mvarDicts = [x[1] for x in mvarGuessTypes]
    r = []
    seen = []
    initialToToss = []
    for mvarDict in mvarDicts:
        if not mvarDict in seen:
            #print('BEF TEST IT', ref, alt)
            if bothMvarAllelesContributeToAlt(mvarDict, ref, alt, initRef, initAlt):
                #print('TEST IT', ref, alt)
                r.append( ('allele1Seq', mkNewMvarDict(mvarDict, ref, alt), True) )
                initialToToss.append(mvarDict)
        seen.append(mvarDict)
    return r, initialToToss
    
def handleSingleMvarGuessMulti(guessType, mvarGuessTypes, ref, alt, initRef, initAlt):
#    print('AAAAAAH')
    mergedMatches, initialToToss = expandWhenMatchesBoth(mvarGuessTypes, guessType, ref, alt, initRef, initAlt)
    delMatches = []
    # if guessType == 'del':
    #     delMatches = vcfMvarDepthDelFuncs.findDelMatches(initRef, initAlt, 
    matches = []
    for m in mergedMatches:
        matches.append(m)
#    print('HERE', ref, alt, initRef, initAlt, mvarGuessTypes, guessType)
    for mvarKey, mvarDat in mvarGuessTypes:
        if not mvarDat in initialToToss:
            if evalMvarMatch(False, guessType, mvarKey, mvarDat, ref, alt, initRef, initAlt):
                matches.append( (mvarKey, mvarDat, False) )
    if len(matches) != 1:
        # check hom case where alleles are the same, and both will match the allele
        homCase = False
        if len(matches) == 2:
            alleleKeys = [x[0] for x in matches]
            if 'allele1Seq' in alleleKeys and 'allele2Seq' in alleleKeys and matches[0][1] == matches[1][1]:
                homCase = True
        if not homCase:
            print('these', matches)
            i = 1/0
    return matches[0]

def handleSingleMvarGuess(mvarGuessTypes, guessType, vcfVar, ref, alt, initRef, initAlt):
    """There is one mvar guess type (snp, del, etc)"""
    if len(mvarGuessTypes[guessType]) == 1:
#        print('OK', ref, alt, guessType)
        if guessType != vcfVar:
            print('SHIT', guessType, vcfVar, ref, alt, mvarGuessTypes)
            i = 1/0
        initMatches, initialToToss = expandWhenMatchesBoth(mvarGuessTypes[guessType], guessType, ref, alt, initRef, initAlt)
        # print('ON NO', initialToToss)
        # print( mvarGuessTypes[guessType][0] )
        if mvarGuessTypes[guessType][0][1] in initialToToss:
 #           print("HERE I MA")
            return handleSingleMvarGuessSingle(guessType, initMatches[0], ref, alt, initRef, initAlt)
        else:
  #          print('HELP ME')
            return handleSingleMvarGuessSingle(guessType, list(mvarGuessTypes[guessType][0]) + [False], ref, alt, initRef, initAlt)
    else: # multilpe del, snp, ins guesses
        if guessType != vcfVar:
            i = 1/0
        #print('SHOOOOO', mvarGuessTypes[guessType], ref, alt, initRef, initAlt)
        return handleSingleMvarGuessMulti(guessType, mvarGuessTypes[guessType], ref, alt, initRef, initAlt)

def handleMultiMvarGuess(mvarGuessTypes, vcfVar, ref, alt, initRef, initAlt):
    if not vcfVar in mvarGuessTypes:
        i = 1/0
    return handleSingleMvarGuess(mvarGuessTypes, vcfVar, vcfVar, ref, alt, initRef, initAlt)

def matchMvarData(ref, alt, mvarData, initRef, initAlt):
    """Match mvar data to this alt.
       Use with all mvarData instances,
       but is only pertinent to cases
       where the vcf line matches multiple
       mvar lines
    """
    mvarDicts = [ buildMvarData(d) for d in mvarData.split(',') ]
    vcfVar = guessVcfVar(ref, alt)
    if vcfVar == 'del':
        return vcfMvarDepthDelFuncs.findDelMatches(initRef, initAlt, mvarDicts, ref, alt)
    if vcfVar == 'ins':
        return vcfMvarDepthInsFuncs.findInsMatches(initRef, initAlt, mvarDicts, ref, alt)

    return vcfMvarDepthSubsFuncs.findSubsMatches(ref, alt, initRef, initAlt, mvarDicts)

    # else snp/subs
    # get's complicated b/c I might need to make a new mvar dataset in some cases
    initMatches, initialToToss = expandWhenMatchesBoth(mvarGuessTypes[guessType], guessType, ref, alt, initRef, initAlt)
    i = 1/0

    # mvarGuessTypes = guessMvarTypes(mvarDicts)
    # if len(mvarGuessTypes) == 1:
    #     mvarGuessType = list(mvarGuessTypes.keys())[0]
    #     return handleSingleMvarGuess(mvarGuessTypes, mvarGuessType, vcfVar, ref, alt, initRef, initAlt)
    # else:
    #     return handleMultiMvarGuess(mvarGuessTypes, vcfVar, ref, alt, initRef, initAlt)



    # # mvarLineCount = getMvarLen(mvarData)
    # # if mvarLineCount == 1:
    # #     return buildMvarData(mvarData)
    
    # # else there is more than one mvar line
    # # for this vcf line
    # vcfVar = guessVcfVar(ref, alt)
    # mvarDicts = [ buildMvarData(d) for d in mvarData.split(',') ]
    # mvarGuessTypes = guessMvarTypes(mvarDicts)
    # #print(mvarGuessTypes)
    # if len(mvarGuessTypes) == 1:
    #     mvarGuessType = list(mvarGuessTypes.keys())[0]
    #     return handleSingleMvarGuess(mvarGuessTypes, mvarGuessType, vcfVar, ref, alt, initRef, initAlt)
    # else:
    #     return handleMultiMvarGuess(mvarGuessTypes, vcfVar, ref, alt, initRef, initAlt)

    # if len( set(mvarTypes) ) != len(mvarTypes):
    #     if mvarTypes[0] == 'del' and len(mvarTypes) == 2:
    #         if mvarDicts[0]['reference'] == ref[1:]:
    #             return mvarDicts[0]
    #         elif mvarDicts[1]['reference'] == ref[1:]:
    #             return mvarDicts[1]
    #         else:
    #             print('cannot decide which mvar deletion line')
    #             print(ref, alt)
    #             print(mvarTypes)
    #             i = 1/0
    #     elif mvarTypes[0] == 'sub' and len(mvarTypes) == 2:
    #         if alt in [mvarDicts[0][x] for x in ('allele1Seq', 'allele2Seq')] and not alt in [mvarDicts[1][x] for x in ('allele1Seq', 'allele2Seq')]:
    #             return mvarDicts[0]
    #         elif alt in [mvarDicts[1][x] for x in ('allele1Seq', 'allele2Seq')] and not alt in [mvarDicts[0][x] for x in ('allele1Seq', 'allele2Seq')]:
    #             return mvarDicts[1]
    #         else:
    #             i = 1/0
    #     else:
    #         print('cannot decide which mvar line')
    #         print(ref, alt)
    #         print(mvarTypes)
    #         i = 1/0
    # mvarLs = []
    # for d, mvarGuess in zip(mvarDicts, mvarGuessTypes):
    #     print('debug2', mvarGuess, d, vcfVar)
    #     if vcfVar in d['varType'] or mvarGuess == vcfVar:
    #         mvarLs.append(d)
    # print('debug new', mvarLs)
    # print('debug new', len(mvarLs))
    # if len(mvarLs) != 1:
    #     print( list(set(mvarTypes)), ['del'] )
    #     if len(set(mvarTypes)) == 2 and 'del' in mvarTypes and 'sub' in mvarTypes:
    #         if mvarTypes[0] == 'del' and len(ref) == 2 and len(alt) == 1 and ref[0] == alt[0] and ref[1] == mvarDicts[0]['reference'] and (not mvarDicts[0]['allele1Seq'] or not mvarDicts[0]['allele2Seq']):
    #             return mvarDicts[0]
    #         elif mvarTypes[1] == 'del' and len(ref) == 2 and len(alt) == 1 and ref[0] == alt[0] and ref[1] == mvarDicts[1]['reference'] and (not mvarDicts[1]['allele1Seq'] or not mvarDicts[1]['allele2Seq']):
    #             return mvarDicts[1]
    #         elif mvarTypes[0] == 'sub' and andmvarDicts[0]['reference'] == ref and alt in (mvarDicts[0]['allele1Seq'], mvarDicts[0]['allele2Seq']):
    #             return mvarDicts[0]
    #         elif mvarTypes[1] == 'sub' and mvarDicts[1]['reference'] == ref and alt in (mvarDicts[1]['allele1Seq'], mvarDicts[1]['allele2Seq']):
    #             return mvarDicts[1]
    #         elif mvarTypes[0] == 'del' and mvarDicts[0]['reference'] == ref[1:] and ((mvarDicts[0]['allele1Seq'] == ref[1:] and not mvarDicts[0]['allele2Seq']) or (mvarDicts[0]['allele2Seq'] == ref[1:] and not mvarDicts[0]['allele1Seq'])):
    #             return mvarDicts[0]
    #         elif mvarTypes[1] == 'del' and mvarDicts[1]['reference'] == ref[1:] and ((mvarDicts[1]['allele1Seq'] == ref[1:] and not mvarDicts[1]['allele2Seq']) or (mvarDicts[1]['allele2Seq'] == ref[1:] and not mvarDicts[1]['allele1Seq'])):
    #             return mvarDicts[1]            
    #         else:
    #             i = 1/0
    #     elif mvarDicts[0]['reference'] == ref[1:] and ((mvarDicts[0]['allele1Seq'] == alt[1:] and not mvarDicts[0]['allele2Seq']) or (mvarDicts[0]['allele2Seq'] == alt[1:] and not mvarDicts[0]['allele1Seq'])):
    #         return mvarDicts[0]
    #     elif mvarDicts[1]['reference'] == ref[1:] and ((mvarDicts[1]['allele1Seq'] == alt[1:] and not mvarDicts[1]['allele2Seq']) or (mvarDicts[1]['allele2Seq'] == alt[1:] and not mvarDicts[1]['allele1Seq'])):
    #         return mvarDicts[1]
    #     else:
    #         print(ref, alt)
    #         print(len(d))
    #         print('debug2', vcfVar)
    #         print(vcfVar, mvarDicts)
    #         print('cannot match mvar line')
    #         print(ref, alt, mvarData)
    #         i = 1/0
    # return mvarLs[0]

# def matchMvarData(ref, alt, mvarData):
#     """Match mvar data to this alt.
#        Use with all mvarData instances,
#        but is only pertinent to cases
#        where the vcf line matches multiple
#        mvar lines
#     """
#     mvarLineCount = getMvarLen(mvarData)
#     if mvarLineCount == 1:
#         return buildMvarData(mvarData)
    
#     # else there is more than one mvar line
#     # for this vcf line
#     vcfVar = guessVcfVar(ref, alt)
#     print('debug', vcfVar)
#     mvarDicts = [ buildMvarData(d) for d in mvarData.split(',') ]
#     mvarGuessTypes = [ guessMvarType(d) for d in mvarDicts ]
#     mvarTypes = [ d['varType'] for d in mvarDicts ]

#     print(mvarData)

#     if len( set(mvarTypes) ) != len(mvarTypes):
#         if mvarTypes[0] == 'del' and len(mvarTypes) == 2:
#             if mvarDicts[0]['reference'] == ref[1:]:
#                 return mvarDicts[0]
#             elif mvarDicts[1]['reference'] == ref[1:]:
#                 return mvarDicts[1]
#             else:
#                 print('cannot decide which mvar deletion line')
#                 print(ref, alt)
#                 print(mvarTypes)
#                 i = 1/0
#         else:
#             print('cannot decide which mvar line')
#             print(ref, alt)
#             print(mvarTypes)
#             i = 1/0
#     mvarLs = []
#     for d, mvarGuess in zip(mvarDicts, mvarGuessTypes):
#         print('debug2', mvarGuess, d)
#         if vcfVar in d['varType'] or mvarGuess == vcfVar:
#             mvarLs.append(d)
#     if len(mvarLs) != 1:
#         print(len(d))
#         print('debug2', vcfVar)
#         print(vcfVar, mvarDicts)
#         print('cannot match mvar line')
#         print(ref, alt, mvarData)
#         i = 1/0
#     return mvarLs[0]

def buildMvarData(mvarData):
    return {f:v for f,v in zip(mvarVcfFields.fields, mvarData.split('/')) }

def assertMvarSingle(mvarData):
    if len( mvarData.split(',') ) != 1:
        print('should have one mvar entry')
        print(mvarData)
        i = 1/0

def getMvarLen(mvarData):
    return len( mvarData.split(',') )

def printHomAlt(sp, totDepth, altCounts, refDepth, mvarData):
    """Might not really be hom."""
    #print('here', sp)
    #assertMvarSingle(mvarData)
    # a1,a2 = altCounts.split(',')
    # altCount = a1
    # if a1 != a2:
    #     print('should be hom alt')
    #     print(sp)
    #     i = 1/0
    vcfRef, vcfAlt = sp[3:5]
    initRef, initAlt = getInitRefAndAlt(sp)
    mvarAlleleKey, mvarDict = matchMvarData(vcfRef, vcfAlt, mvarData, initRef, initAlt)
    refDepth = mvarDict['referenceAlleleReadCount']
    altCount = mvarDict[mvarAlleleKey.replace('Seq', 'ReadCount')]

#    print(mvarDict)

    # if mvarDict['allele1ReadCount'] != mvarDict['allele2ReadCount']:
    #     i = 1/0
    # if mvarDict['allele1ReadCount'] != altCount:
    #     print(mvarDict['allele1ReadCount'], altCount)
    #     i = 1/0
    # if mvarDict['allele1ReadCount-T1'] != mvarDict['allele2ReadCount-T1']:
    #     i = 1/0    

    homIsWrong = 'F'
    totDepth = mvarDict['totalReadCount']
    if mvarDict['allele1ReadCount'] == mvarDict['allele2ReadCount'] and mvarDict['allele1Seq'] == mvarDict['allele2Seq']:
        effectiveDepth = sum([int(x) for x in (refDepth, altCount)])
        tumorEffectiveDepth = int(mvarDict['allele1ReadCount-T1']) + int(mvarDict['referenceAlleleReadCount-T1'])

    else:
        effectiveDepth = sum([int(mvarDict[x]) for x in ('allele1ReadCount', 'allele2ReadCount')])
        tumorEffectiveDepth = int(mvarDict['allele1ReadCount-T1']) + int(mvarDict['allele2ReadCount-T1'])
        homIsWrong = 'T'
    outSp = evalResult(sp, totDepth, altCount, refDepth, effectiveDepth, 'F', homIsWrong,
                       int(mvarDict['totalReadCount-T1']),
                       int(mvarDict['allele1ReadCount-T1']),
                       int(mvarDict['referenceAlleleReadCount-T1']), 
                       tumorEffectiveDepth, )
    return outSp

def hasAlleleSubset(mvarDict, totDepth, tumorTotDepth):
    """ex TG/CGC/CG
          ref/allele1/allele2
       not all nested result in alt sums greater than the total
       not all un-nested are safe from having the sum be greater
       This is now just a test of the sum of allele counts being greater than the total 
       read count.
    """
    s1 = mvarDict['allele1Seq']
    s2 = mvarDict['allele2Seq']
    lenS1 = len(s1)
    lenS2 = len(s2)
    acc, accTumor = 0, 0
    if mvarDict['allele1ReadCount']:
        acc += int(mvarDict['allele1ReadCount'])
    if mvarDict['allele2ReadCount']:
        acc += int(mvarDict['allele2ReadCount'])

    if mvarDict['allele1ReadCount-T1']:
        accTumor += int(mvarDict['allele1ReadCount-T1'])
    if mvarDict['allele2ReadCount-T1']:
        accTumor += int(mvarDict['allele2ReadCount-T1'])

    if acc > int(totDepth) or accTumor > int(tumorTotDepth):
        if lenS1 < lenS2 and lenS1 != 0:
            if s1 in s2:
                return True
        elif lenS1 > lenS2 and lenS2 != 0:
            if s2 in s1:
                return True

    return False 

def printRefAndAltForNoSubset(sp, mvarDict):
    """I cannot trust the vcf allele depths here."""
    #print('deubg3', mvarDict)
    if mvarDict['reference'] == mvarDict['allele1Seq']:
        altCount = mvarDict['allele2ReadCount']
        altCountTumor = mvarDict['allele2ReadCount-T1']
    elif mvarDict['reference'] == mvarDict['allele2Seq']:
        altCount = mvarDict['allele1ReadCount']
        altCountTumor = mvarDict['allele1ReadCount-T1']
    else:
        print(sp, 'should be one ref, one alt')
        i = 1/0

    if mvarDict['referenceAlleleReadCount'] == '0':
        zygosityWrong = 'T'

    effectiveDepth = sum([int(x) for x in (mvarDict['referenceAlleleReadCount'], altCount)])
    effectiveDepthTumor = sum([int(x) for x in (mvarDict['referenceAlleleReadCount-T1'], altCountTumor)])
    #print('yo', effectiveDepth, altCount )
    return effectiveDepth, altCount, mvarDict['referenceAlleleReadCount'], effectiveDepthTumor, altCountTumor, mvarDict['referenceAlleleReadCount-T1']

def guessMvarVars(mvarData):
    """For each allele seq, guess ins/del.snp."""
    refSeqLen = len(mvarData['reference'])
    status1, status2 = 'snp', 'snp'

    if len(mvarData['allele1Seq']) > refSeqLen:
        status1 = 'ins'
    elif len(mvarData['allele1Seq']) < refSeqLen:
        status1 = 'del'

    if len(mvarData['allele2Seq']) > refSeqLen:
        status2 = 'ins'
    elif len(mvarData['allele2Seq']) < refSeqLen:
        status2 = 'del'

    if status1 == status2:
        print('cant distinguish mvar alleles')
        print(mvarData)
        i = 1/0

    return status1, status2

def printRefAndAltForSubset(sp, mvarKey, mvarData, totDepth):
    """Here the sum of a1 and a2 depths is greater than the total depth.
       This is because one allele is always present b/c it is a
       subset of the other.
       You cannot use the call (ex 0/1) to determine the order of the 
       vcf allele counts, so you can't use the vcf allele counts at all.
       You must get the count from the mVar line.
    """
    # what is the alt count?
    # if one is not ref?
    refVcf, altVcf = sp[3:5]
    vcfVar = guessVcfVar(refVcf, altVcf)
    
    #mvar1Var, mvar2Var = guessMvarVars(mvarData)
    if mvarKey == 'allele1Seq':
        altCount = mvarData['allele1ReadCount']
        altCountTumor = mvarData['allele1ReadCount-T1']
    elif mvarKey == 'allele2Seq':
        altCount = mvarData['allele2ReadCount']
        altCountTumor = mvarData['allele2ReadCount-T1']
    else:
        i = 1/0
    effectiveDepth = max([int(mvarData[x]) for x in ('allele1ReadCount', 'allele2ReadCount')])#int(totDepth)
    effectiveDepthTumor = max([int(mvarData[x]) for x in ('allele1ReadCount-T1', 'allele2ReadCount-T1')]) #int(altCountTumor)
    return effectiveDepth, altCount, effectiveDepthTumor, altCountTumor

def printRefAndAlt(sp, totDepth, altCounts, refDepth, mvarData):
    vcfRef, vcfAlt = sp[3:5]
    initRef, initAlt = getInitRefAndAlt(sp)
    mvarAlleleSeqKey, mvarDict = matchMvarData(vcfRef, vcfAlt, mvarData, initRef, initAlt)

    # bad
    # can't always trust the vcf call!
    # a1,a2 = altCounts.split(',')

    # if one alt is a subset of another,
    # the total reads considered changes
    if hasAlleleSubset(mvarDict, totDepth, int(mvarDict['totalReadCount-T1'])):
        (effectiveDepth, altCount,
         tumorEffectiveDepth, tumorAltCount) = printRefAndAltForSubset(sp, mvarAlleleSeqKey,
                                                                       mvarDict, totDepth)
#        print('oh no', altCount)
    else:
        # can't trust ref depth
        effectiveDepth, altCount, refDepth, tumorEffectiveDepth, tumorAltCount, tumorRefDepth = printRefAndAltForNoSubset(sp, mvarDict)
        #print('here', altCount)
    
    zygosityWrong = 'F'

    tumorTotDepth = int(mvarDict['totalReadCount-T1'])
    tumorRefDepth = int(mvarDict['referenceAlleleReadCount-T1'])
    return evalResult(sp, totDepth, altCount, refDepth, effectiveDepth, 'F', zygosityWrong,
                      tumorTotDepth, tumorAltCount, tumorRefDepth, tumorEffectiveDepth, )

def printMultiAlt(sp, totDepth, altCounts, refDepth, mvarData, fout):
    # OLD_MULTIALLELIC=6:29910558:T/C/A
    assertMvarSingle(mvarData)
    ref, alt1, alt2 = sp[-3].split('OLD_MULTIALLELIC=')[1].split(';')[0].split(':')[-1].split('/')
    a1,a2 = altCounts.split(',')
    if sp[4] == alt1:
        altCount = a1
    elif sp[4] == alt2:
        altCount = a2
    else:
         print('bad alt ls', sp)
         i = 1/0

    mvarDict = buildMvarData(mvarData)
    if not altCount in (mvarDict['allele1ReadCount'], mvarDict['allele2ReadCount']):
        i = 1/0

    effectiveDepth = sum( (int(a1), int(a2), int(refDepth)) )
    effectiveDepthTumor = sum( [int(mvarDict[x]) for x in ('allele1ReadCount-T1', 'allele2ReadCount-T1', 'referenceAlleleReadCount-T1')] )
    tumorRefDepth = int(mvarDict['totalReadCount-T1'])
    tumorTotDepth = int(mvarDict['referenceAlleleReadCount-T1']) 
    if altCount == int(mvarDict['allele1ReadCount']) and altCount != int(mvarDict['allele2ReadCount']):
        tumorAltCount = mvarDict['allele1ReadCount']
    elif altCount == int(mvarDict['allele2ReadCount']) and altCount != int(mvarDict['allele1ReadCount']):
        tumorAltCount = mvarDict['allele2ReadCount']
    else:
        i = 1/0
    evalAndPrint(sp, totDepth, altCount, refDepth, effectiveDepth, 'F', 'F', 
                 tumorTotDepth, tumorAltCount, tumorRefDepth, tumorEffectiveDepth,
                 fout)

def printMissingCall(sp, totDepth, altCounts, refDepth, mvarData):
    #assertMvarSingle(mvarData)
    # one alt count is .
    # a1,a2 = altCounts.split(',')
    # if '.' == a1:
    #     altCount = a2
    # elif '.' == a2:
    #     altCount = a1
    # else:
    #      print('bad alt ls', sp)
    #      i = 1/0

    vcfRef, vcfAlt = sp[3:5]
    initRef, initAlt = getInitRefAndAlt(sp)
    mvarAlleleKey, mvarDict = matchMvarData(vcfRef, vcfAlt, mvarData, initRef, initAlt)
    #mvarDict = buildMvarData(mvarData)
    
    if mvarAlleleKey == 'allele1Seq':
        altCount = mvarDict['allele1ReadCount']
        tumorAltCount = mvarDict['allele1ReadCount-T1']
    elif mvarAlleleKey == 'allele2Seq':
        altCount = mvarDict['allele2ReadCount']
        tumorAltCount = mvarDict['allele2ReadCount-T1']
    else:
        i = 1/0

    refDepth = mvarDict['referenceAlleleReadCount']

    # if not altCount in (mvarDict['allele1ReadCount'], mvarDict['allele2ReadCount']):
    #     i = 1/0
    # if not '' in (mvarDict['allele1ReadCount'], mvarDict['allele2ReadCount']):
    #     i = 1/0

    # if str(altCount) == mvarDict['allele1ReadCount']:
        
    # elif str(altCount) == mvarDict['allele2ReadCount']:
        
    # else:
    #     i = 1/0

    effectiveDepth = int(altCount) + int(refDepth)
    tumorRefDepth = int(mvarDict['referenceAlleleReadCount-T1'])
    tumorEffectiveDepth = int(tumorAltCount) + int(tumorRefDepth)
    tumorTotDepth = int(mvarDict['totalReadCount-T1'])

    isReallyHet = 'F'
    if int(refDepth) > 0:
        isReallyHet = 'T'

    return evalResult(sp, totDepth, altCount, refDepth, effectiveDepth, 'T', isReallyHet,
                      tumorTotDepth, tumorAltCount, tumorRefDepth, tumorEffectiveDepth, )

# is it really fake missing?
# def doMultiWithFakeMissing(call, sp, totDepth, altCounts, refDepth, mvarData, fout):
#     # no call is caused by vt when variants are split
#     # ref is near zero
#     # use call to see which alt allele to use for this pos
#     a1,a2 = altCounts.split(',')
#     # mvarLen = getMvarLen(mvarData)
#     # if mvarLen > 1:
#     #     i = 1/0
#     mvarDict = buildMvarData(mvarData)
#     if hasAlleleSubset(mvarDict, totDepth):
#         effectiveDepth, altCount = printRefAndAltForSubset(sp, mvarData, totDepth)
#     else:
#         if '|' in call:
#             c1,c2 = call.split('|')
#         elif '/' in call:
#             c1,c2 = call.split('/')
#         else:
#             i = 1/0

#         if not c1 in ('.','1') or not c2 in ('.','1'):
#             i = 1/0
#         if c1 == '1':
#             altCount = a1
#         elif c2 == '1':
#             altCount = a2
#         else:
#             i = 1/0

#         effectiveDepth = int(a1) + int(a2) + int(refDepth)
#     evalAndPrint(sp, totDepth, altCount, refDepth, effectiveDepth, 'F', 'F', fout)    

def getInitRefAndAlt(sp):
    if 'OLD_MULTIALLELIC' in sp[-3]:
        initRef, multAllele1, multAllele2 = sp[-3].split('OLD_MULTIALLELIC=')[1].split(';')[0].split(':')[-1].split('/')
    if 'OLD_CLUMPED=' in sp[-3]:
        initRef, useAllele = sp[-3].split('OLD_CLUMPED=')[1].split(';')[0].split(':')[-1].split('/')
    elif 'OLD_VARIANT' in sp[-3]:
        initRef, useAllele = sp[-3].split('OLD_VARIANT=')[1].split(';')[0].split(':')[-1].split('/')
    else:
        initRef, useAllele = sp[3], sp[4]

    if 'OLD_MULTIALLELIC' in sp[-3]:
        if not useAllele in (multAllele1, multAllele2):
            print(sp)
            print(useAllele)
            i = 1/0

    return initRef, useAllele

def findAltCountFromOldVcf(sp, mvarDict):
    """
    TARGET-30-PASNZU-10A-01D.10
    10      8007559 .       G       C
    1/.:.:PASS
    OLD_MULTIALLELIC=10:8007559:GT/CC/GC;OLD_CLUMPED=10:8007559:GT/CC
    66:27,33:0:GT/CC/GC/het-alt/complex/27/33/0/66/44/26/0/80
    You cannot trust the order of call here: 1|.
    The order in the multi allele vcf is not the order in mvar.
    You must check alleles.
    """
    initRef, useAllele = getInitRefAndAlt(sp)

    # initRef, multAllele1, multAllele2 = sp[-3].split('OLD_MULTIALLELIC=')[1].split(';')[0].split(':')[-1].split('/')
    # if 'OLD_CLUMPED=' in sp[-3]:
    #     initRef, useAllele = sp[-3].split('OLD_CLUMPED=')[1].split(';')[0].split(':')[-1].split('/')
    # elif 'OLD_VARIANT' in sp[-3]:
    #     initRef, useAllele = sp[-3].split('OLD_VARIANT=')[1].split(';')[0].split(':')[-1].split('/')
    # else:
    #     useAllele = sp[4]
    

    if useAllele == mvarDict['allele1Seq']:
        return mvarDict['allele1ReadCount'], mvarDict['allele1ReadCount-T1']
    elif useAllele == mvarDict['allele2Seq']:
        return mvarDict['allele2ReadCount'], mvarDict['allele2ReadCount-T1']
    else: # I see a case where the blank allele is not the deletion?
        if len(useAllele) == 1 and (not mvarDict['allele2Seq'] or not mvarDict['allele1Seq']) and len(useAllele) < len(sp[3]):
            # have a deletion that is coded differently in the mvar file
            # use the empty mvar allele
            if not (not mvarDict['allele2Seq'] and not mvarDict['allele1Seq']):
                if not mvarDict['allele2Seq']:
                    return mvarDict['allele2ReadCount'], mvarDict['allele2ReadCount-T1']
                if not mvarDict['allele1Seq']:
                    return mvarDict['allele1ReadCount'], mvarDict['allele1ReadCount-T1']
                i = 1/0
            else:
                i = 1/0
        elif len(sp[3]) == 1 and sp[3] == useAllele[0] and sp[3] == initRef[0] and useAllele[1:] in ( mvarDict['allele2Seq'], mvarDict['allele1Seq'] ) and initRef[1:] == mvarDict['reference']:
            if useAllele[1:] == mvarDict['allele1Seq']:
                return mvarDict['allele1ReadCount'], mvarDict['allele1ReadCount-T1']
            elif useAllele[1:] == mvarDict['allele2Seq']:
                return mvarDict['allele2ReadCount'], mvarDict['allele2ReadCount-T1']
        elif initRef[0] == useAllele[0] and \
             len(initRef) > 1 and len(useAllele) > 1 and \
             initRef[1:] == mvarDict['reference'] and \
             useAllele[1:] in ( mvarDict['allele2Seq'], mvarDict['allele1Seq'] ):
            if useAllele[1:] == mvarDict['allele1Seq']:
                return mvarDict['allele1ReadCount'], mvarDict['allele1ReadCount-T1']
            elif useAllele[1:] == mvarDict['allele2Seq']:
                return mvarDict['allele2ReadCount'], mvarDict['allele2ReadCount-T1']
        elif sp[1] == '49334018' and sp[0] == '22' and sp[3] == 'GAGA' and sp[4] == 'AAG' and mvarDict['reference'] == 'AGA' and mvarDict['allele1Seq'] == 'AG' and not mvarDict['allele2Seq']:
            # too hard to figure out
            return mvarDict['allele1ReadCount'], mvarDict['allele1ReadCount-T1']
        elif sp[1] == '97291797' and sp[0] == '12' and sp[3] == 'GGGC' and sp[4] == 'CGG' and mvarDict['reference'] == 'GGC' and mvarDict['allele1Seq'] == 'GG' and not mvarDict['allele2Seq']:
            # too hard to figure out
            return mvarDict['allele1ReadCount'], mvarDict['allele1ReadCount-T1']
        elif sp[4][1:] == mvarDict['allele1Seq'] and sp[3][1:] == mvarDict['reference'] and not mvarDict['allele2Seq']:
            return mvarDict['allele1ReadCount'], mvarDict['allele1ReadCount-T1']
        elif sp[4][1:] == mvarDict['allele2Seq'] and sp[3][1:] == mvarDict['reference'] and not mvarDict['allele1Seq']:
            return mvarDict['allele2ReadCount'], mvarDict['allele2ReadCount-T1']            
        else:
            print(sp[4][1:], mvarDict['allele1Seq'], initRef[1:], mvarDict['reference'])
            print(sp)
            print(initRef)
            print(  sp[3] == initRef[0], useAllele[1:] in ( mvarDict['allele2Seq'], mvarDict['allele1Seq'] ),  )
            print(useAllele, multAllele1, multAllele2)
            i = 1/0

def printMultiWithFakeMissing(sp, totDepth, altCounts, refDepth, mvarData):
    # there can be multiple mvar sites here
    # subset problem here?
    vcfRef, vcfAlt = sp[3:5]
    #print('debug', sp)
    
    initRef, initAlt = getInitRefAndAlt(sp)
    mvarAlleleSeqKey, mvarDict = matchMvarData(vcfRef, vcfAlt, mvarData, initRef, initAlt)
#    print(mvarDict)
    if hasAlleleSubset(mvarDict, totDepth, int(mvarDict['totalReadCount-T1'])):
        (effectiveDepth, altCount,
         tumorEffectiveDepth, tumorAltCount) = printRefAndAltForSubset(sp, mvarAlleleSeqKey,
                                                                       mvarDict, totDepth)
    else:
        # alt sum not screwed up
        totDepth = mvarDict['totalReadCount']
        if mvarDict['zygosity'] == 'hom':
            refDepth = mvarDict['referenceAlleleReadCount']
            tumorRefDepth = int(mvarDict['referenceAlleleReadCount-T1'])
            altCount = mvarDict['allele1ReadCount']
            tumorAltCount = mvarDict['allele1ReadCount-T1']
            effectiveDepth = int(mvarDict['allele1ReadCount'])
            tumorEffectiveDepth = int(mvarDict['allele1ReadCount-T1'])
        elif mvarDict['zygosity'] in ('het-ref',):
            # I have the right allele
            refDepth = mvarDict['referenceAlleleReadCount']
            tumorRefDepth = int(mvarDict['referenceAlleleReadCount-T1'])

            if mvarAlleleSeqKey == 'allele1Seq':
                altCount = mvarDict['allele1ReadCount']
                tumorAltCount = mvarDict['allele1ReadCount-T1']
            elif mvarAlleleSeqKey == 'allele2Seq':
                altCount = mvarDict['allele2ReadCount']
                tumorAltCount = mvarDict['allele2ReadCount-T1']
            else: i = 1/0
        
            effectiveDepth = int(mvarDict['allele1ReadCount']) + int(mvarDict['allele2ReadCount'])
            tumorEffectiveDepth = int(mvarDict['allele1ReadCount-T1']) + int(mvarDict['allele2ReadCount-T1'])
        else:
            #print('yo2')
            a1 = mvarDict['allele1ReadCount']
            a2 = mvarDict['allele2ReadCount']
            refDepth = mvarDict['referenceAlleleReadCount']
            tumorRefDepth = int(mvarDict['referenceAlleleReadCount-T1'])
            if mvarAlleleSeqKey == 'allele1Seq':
                altCount = a1
                tumorAltCount = mvarDict['allele1ReadCount-T1']
            elif mvarAlleleSeqKey == 'allele2Seq':
                altCount = a2
                tumorAltCount = mvarDict['allele2ReadCount-T1']
            else:
                # figure it out from the initial vcf
                #altCount, tumorAltCount = findAltCountFromOldVcf(sp, mvarDict)
                # print(sp)
                # print(mvarDict['zygosity'])
                # print(mvarDict)
                i = 1/0
            #print(mvarDict)
            effectiveDepth = 0
            if a1:
                effectiveDepth += int(a1)
            if a2:
                effectiveDepth += int(a2) # + int(refDepth)

            tumorEffectiveDepth = 0
            for x in ('allele1ReadCount-T1', 'allele2ReadCount-T1'):
                if mvarDict[x]:
                    tumorEffectiveDepth += int(mvarDict[x])

            #sum([int(mvarDict[x]) for x in ('allele1ReadCount-T1', 'allele2ReadCount-T1',)])
            
    hasNoCall = 'T'
    if mvarDict['allele1ReadCount'].strip() and mvarDict['allele2ReadCount'].strip():
        hasNoCall = 'F'
    tumorTotDepth = mvarDict['totalReadCount-T1']
    tumorRefDepth = mvarDict['referenceAlleleReadCount-T1']
        
    printLs = evalResult(sp, totDepth, altCount, refDepth, effectiveDepth, hasNoCall, 'F',
                         tumorTotDepth, tumorAltCount, tumorRefDepth, tumorEffectiveDepth)
    return printLs
