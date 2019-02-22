"""These functions are used to match vcf vars with mvar depths for substitutions.
   See the README.md for a full description.
"""
import vcfMvarDepthMatchByDistance

def mvarDictHasHomSub(mvarDict):
    return len(mvarDict['reference']) == len(mvarDict['allele1Seq']) \
           and mvarDict['allele1Seq'] == mvarDict['allele2Seq']

def matchHomSub(ref, alt, mvarDict):
    if ref == mvarDict['reference'] and alt == mvarDict['allele1Seq']:
        return [ ('allele1Seq', mvarDict) ]

# def bothMvarAllelesContributeToAlt(mvarDict, ref, alt, initRef, initAlt):
#     """Ex AG > CA/CG in mvar file. Need counts for both alleles.
#        Or C>T vcf and CGA/TGA/TG/.
#     """
#     return mvarDict['allele1Seq'] != mvarDict['allele2Seq'] and mvarDict['allele1Seq'] and mvarDict['allele2Seq'] and \
#            ( (len(mvarDict['reference']) == 2 and len(mvarDict['allele1Seq']) == 2 and len(mvarDict['allele2Seq']) == 2 and \
#               len(ref) == 1 and len(alt) == 1 and mvarDict['reference'][0] == ref and \
#               mvarDict['allele1Seq'][0] == alt and mvarDict['allele2Seq'][0] == alt)
#              or \
#              (ref == mvarDict['reference'][0]
#               and alt == mvarDict['allele1Seq'][0] and alt == mvarDict['allele2Seq'][0]) )

# def expandWhenMatchesBoth(mvarGuessTypes, varType, ref, alt, initRef, initAlt):
#     """It is very complicated now. When a single mar line with 2 alt alleles
#        that can match the vcf variant, collapse the counts into a new mvarDict.
#        Can't use mvardict as blend and initial.
#     """
#     #print(ref, alt, varType)
#     if varType != 'snp':
#         return [], []
#     mvarDicts = [x[1] for x in mvarGuessTypes]
#     r = []
#     seen = []
#     initialToToss = []
#     for mvarDict in mvarDicts:
#         if not mvarDict in seen:
#             #print('BEF TEST IT', ref, alt)
#             if bothMvarAllelesContributeToAlt(mvarDict, ref, alt, initRef, initAlt):
#                 #print('TEST IT', ref, alt)
#                 r.append( ('allele1Seq', mkNewMvarDict(mvarDict, ref, alt), True) )
#                 initialToToss.append(mvarDict)
#         seen.append(mvarDict)
#     return r, initialToToss

#def mvarDictHasSingleSub(mvarDict):
    
def matchSingleSub(ref, alt, mvarDict):
    if ref == mvarDict['reference'] and alt == mvarDict['allele1Seq']:
        return [ ('allele1Seq', mvarDict) ]
    if ref == mvarDict['reference'] and alt == mvarDict['allele2Seq']:
        return [ ('allele2Seq', mvarDict) ]

    matches = []
    # issue when paired with something else
    if ref[0] == alt[0]:
        for a1,a2 in (('allele1Seq', 'allele2Seq'), ('allele2Seq', 'allele1Seq')):
            if len(mvarDict['reference']) == len(mvarDict[a1]) \
               and len(mvarDict['reference']) != len(mvarDict[a2]) \
               and ref[1:] == mvarDict['reference'] and alt[1:] == mvarDict[a1]:
                matches.append( (a1, mvarDict) )
    elif ref[1:] == alt[1:]:
        matchRef, matchAlt = ref[0], alt[0]
        for a1,a2 in (('allele1Seq', 'allele2Seq'), ('allele2Seq', 'allele1Seq')):
            if mvarDict['reference'] == matchRef \
               and mvarDict[a1] == matchAlt and mvarDict[a2] != matchAlt:
                matches.append( (a1, mvarDict) )
    
    return matches

def matchBySubset(ref, alt, mvarDictList):
    """Try subset of vcf ref and alt."""
    matchLs = []
    if len(ref) != len(alt):
        i = 1/0
    subsetMatches = []
    for idx in range(len(ref)):
        if ref[idx:] == alt[idx:]:
            subsetMatches.append( (len(ref[idx:]), idx) )
    if subsetMatches:
        subsetMatches.sort()
        _, idx = subsetMatches[-1]
        for mvarDict in mvarDictList:
            if mvarDict['allele1Seq'] == mvarDict['allele2Seq']:
                i = 1/0
            if mvarDict['reference'] == ref[:idx] \
               and mvarDict['allele1Seq'] == alt[:idx]:
                matchLs.append( ('allele1Seq', mvarDict) )
            if mvarDict['reference'] == ref[:idx] \
               and mvarDict['allele2Seq'] == alt[:idx]:
                matchLs.append( ('allele2Seq', mvarDict) )

    return matchLs

def findSubsMatches(updatedRef, updatedAlt, initRef, initAlt, mvarDictList):
    """Sometimes I need to merge mvar depths into one allele.
       That's why I have updatedRef and updatedAlt.
    """
    #initMatches, initialToToss = expandWhenMatchesBoth(mvarGuessTypes[guessType], guessType, ref, alt, initRef, initAlt)
    matches = []
    for mvarDict in mvarDictList:
        if mvarDictHasHomSub(mvarDict):
            matches += matchHomSub(initRef, initAlt, mvarDict)
        else: #if mvarDictHasSingleSub(mvarDict):
            matches += matchSingleSub(initRef, initAlt, mvarDict)

    if not matches:
        matches += matchBySubset(initRef, initAlt, mvarDictList)
    if len(matches) != 1:
        print('here?')
        matches = vcfMvarDepthMatchByDistance.findMatchesByDistance(initRef, initAlt, mvarDictList)
    if len(matches) != 1:
        print('matches', matches)
        i = 1/0
    return matches[0]

