"""These functions are used to match vcf vars with mvar depths for ins.
   See the README.md for a full description.
"""
import vcfMvarDepthMatchByDistance

def mvarDictHasHomIns(mvarDict):
    return len(mvarDict['reference']) < len(mvarDict['allele1Seq']) and mvarDict['allele1Seq'] == mvarDict['allele2Seq']

def matchHomIns(ref, alt, mvarDict):
    if mvarDict['reference'] == ref[1:] and mvarDict['allele1Seq'] == alt[1:]:
        return [ ('allele1Seq', mvarDict) ]
    return []

def hasOneIns(mvarDict):
    oneIns = len(mvarDict['reference']) < len(mvarDict['allele1Seq'])
    twoIns = len(mvarDict['reference']) < len(mvarDict['allele2Seq'])
    return not (oneIns and twoIns) and (oneIns or twoIns)

def hasTwoIns(mvarDict):
    """Not the same ins"""
    oneIns = len(mvarDict['reference']) < len(mvarDict['allele1Seq'])
    twoIns = len(mvarDict['reference']) < len(mvarDict['allele2Seq'])
    return oneIns and twoIns

def matchOneIns(ref, alt, mvarDict):
    if mvarDict['reference'] == ref[1:] and mvarDict['allele1Seq'] == alt[1:]:
        return [ ('allele1Seq', mvarDict) ]
    if mvarDict['reference'] == ref[1:] and mvarDict['allele2Seq'] == alt[1:]:
        return [ ('allele2Seq', mvarDict) ]
    return []

def matchTwoIns(ref, alt, mvarDict):
    match1 = mvarDict['reference'] == ref[1:] and mvarDict['allele1Seq'] == alt[1:]
    match2 = mvarDict['reference'] == ref[1:] and mvarDict['allele2Seq'] == alt[1:]
    if match1 and match2:
        i = 1/0
    if match1:
        return [ ('allele1Seq', mvarDict) ]
    if match2:
        return [ ('allele2Seq', mvarDict) ]
    return []

# def buildInsMatches(d):
#     match1 = len(d['reference']) < len(d['allele1Seq'])
#     match2 = len(d['reference']) < len(d['allele2Seq'])
#     if match1 and match2:
#         return vcfMvarDepthMatchByDistance.findMatchesByDistance(ref, alt, mvarDictList)
#     if match1:
#         return ('allele1Seq', d)
#     if match2:
#         return ('allele2Seq', d)
#     i = 1/0

# def locateByIns(mvarDictList):
#     """Last resort for ins"""
#     return [buildInsMatches(d) for d in mvarDictList if len(d['reference']) < max([len(d[x]) for x in ('allele1Seq', 'allele2Seq')])]

def findMatchesByUpdatedVcf(matches, updatedRef, updatedAlt):
    newMatches = []
    for m in matches:
        mKey, mvar = m
        if mvar['reference'] == updatedRef and mvar[mKey] == updatedAlt:
            newMatches.append(m)
    return newMatches

def findInsMatches(ref, alt, mvarDictList, updatedRef, updatedAlt):
    """Use cgatools ref and alt to locate possible deletions."""
    matches = []
    for mvarDict in mvarDictList:
        if mvarDictHasHomIns(mvarDict):
            matches += matchHomIns(ref, alt, mvarDict)
        elif hasOneIns(mvarDict):
            matches += matchOneIns(ref, alt, mvarDict)
        elif hasTwoIns(mvarDict):
            matches += matchTwoIns(ref, alt, mvarDict)
    if not matches:
        matches = vcfMvarDepthMatchByDistance.findMatchesByDistance(ref, alt, mvarDictList)
        #locateByIns(mvarDictList)
    if len(matches) > 1:
        matches = findMatchesByUpdatedVcf(matches, updatedRef, updatedAlt)
    if len(matches) != 1:        
        matches = vcfMvarDepthMatchByDistance.findMatchesByDistance(ref, alt, mvarDictList)
    if len(matches) != 1:
        print('matches', matches)
        i = 1/0
    return matches[0]
