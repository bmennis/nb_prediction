"""These functions are used to match vcf vars with mvar depths for del.
   See the README.md for a full description.
"""
import vcfMvarDepthMatchByDistance

def determineWhatIsDeleted(ref, alt):
    """Must take the deletion from the ends.
    """
    altIdx = ref.find(alt)
    if altIdx == 0:
        # the tail is deleted
        # AGGTCTG to A
        # TCCTACACCTACT to TCCTACA so TCCTACA(CCTACT) (want)
        endDelIdx = altIdx + len(alt)
        return ref[endDelIdx:]

    # try dropping the first nuc from both
    # and recompute the match
    # the head is deleted after chopping off first T
    # TCCTACACCTACT to TCCTACT (alt)
    # T(CCTACA)CCTACT (want)
    ref, alt = ref[1:], alt[1:]
    altIdx = ref.find(alt)
    if -1 == altIdx or altIdx+len(alt) != len(ref):
        return '-1'
    return ref[:altIdx]

def mvarDictHasHomBlank(mvarDict):
    return not mvarDict['allele1Seq'] and not mvarDict['allele2Seq']

def matchHomBlank(ref, alt, mvarDict):
    """For hom del, both alt counts are the same, so use one."""
    whatIsDeleted = determineWhatIsDeleted(ref, alt)
    if whatIsDeleted == mvarDict['reference']:
        return [ ('allele1Seq', mvarDict) ]
    return []

def mvarDictHasHomDel(mvarDict):
    """Will not be have a blank for del."""
    return mvarDict['allele1Seq'] == mvarDict['allele2Seq'] \
           and len(mvarDict['reference']) > len(mvarDict['allele1Seq'])

def mvarDictHasBlank(mvarDict):
    return not mvarDict['allele1Seq'] or not mvarDict['allele2Seq']

def matchSingleBlank(ref, alt, mvarDict):
    """One allele is blank, and the other is not."""
    whatIsDeleted = determineWhatIsDeleted(ref, alt)
    matches = []
    if whatIsDeleted == mvarDict['reference']:
        if not mvarDict['allele1Seq']:
            return [ ('allele1Seq', mvarDict) ]
        if not mvarDict['allele2Seq']:
            return [ ('allele2Seq', mvarDict) ]
    else:
        # vcf TGGGGGTGCAAGGTGAG/TG
        # GGGGGTGCAAGGTGAG/G// (I want the G)
        for a1, a2 in ( ('allele1Seq', 'allele2Seq'),
                        ('allele2Seq', 'allele1Seq') ):
            if len(mvarDict['reference']) > len(mvarDict[a2]) \
               and not mvarDict[a1] \
               and mvarDict['reference'].find(mvarDict[a2]) == 0 \
               and mvarDict[a2] == alt[1:] and ref.find(alt) == 0:
                matches.append( (a2, mvarDict) )
    if not matches:
        if not mvarDict['allele1Seq'] and mvarDict['allele2Seq'] == '?':
            matches.append( ('allele1Seq', mvarDict) )
        if not mvarDict['allele2Seq'] and mvarDict['allele1Seq'] == '?':
            matches.append( ('allele2Seq', mvarDict) )

    return matches

def doesVcfAltMatchMvarAltForDelCase(vcfAlt, mvarAlleleKey, mvarDict):
    return vcfAlt[1:] == mvarDict[mvarAlleleKey] \
           and len(mvarDict['reference']) > len(mvarDict[mvarAlleleKey])

def matchNonBlank(ref, alt, mvarDict):
    """No alleles are blank, but this is still a del."""
    matches1 = doesVcfAltMatchMvarAltForDelCase(alt, 'allele1Seq', mvarDict)
    matches2 = doesVcfAltMatchMvarAltForDelCase(alt, 'allele2Seq', mvarDict)
    if matches1 and matches2:
        i = 1/0
    if matches1:
        return [ ('allele1Seq', mvarDict) ]
    if matches2:
        return [ ('allele2Seq', mvarDict) ]
    return []

def matchRefAndAlt(ref, alt, mvarDict):
    if mvarDict['allele1Seq'] == mvarDict['allele2Seq']:
        i = 1/0
    if ref == mvarDict['reference'] and mvarDict['allele1Seq'] == alt:
        return [ ('allele1Seq', mvarDict) ]
    if ref == mvarDict['reference'] and mvarDict['allele2Seq'] == alt:
        return [ ('allele2Seq', mvarDict) ]
    return []

def matchHomNonBlank(ref, alt, mvarDict):
    matches1 = doesVcfAltMatchMvarAltForDelCase(alt, 'allele1Seq', mvarDict)
    if matches1:
        return [ ('allele1Seq', mvarDict) ]
    return []

def buildDelMatches(d, alt):
    if not d['allele1Seq']:
        return ('allele1Seq', d)
    if not d['allele2Seq']:
        return ('allele2Seq', d)

    if alt == d['allele2Seq'] and alt == d['allele1Seq']:
        i = 1/0

    if alt == d['allele1Seq']:
        return ('allele1Seq', d)
    if alt == d['allele2Seq']:
        return ('allele2Seq', d)

    print(d['allele1Seq'])
    print(d['allele2Seq'])
    i = 1/0

def buildDelMatchesFromAlt(d, alt):
    if d['allele1Seq'] == d['allele2Seq']:
        i = 1/0

    if alt == d['allele1Seq']:
        return ('allele1Seq', d)
    if alt == d['allele2Seq']:
        return ('allele2Seq', d)

    print(d['allele1Seq'])
    print(d['allele2Seq'])
    i = 1/0

def locateByDelAndRef(mvarDictList, alt):
    """Last resort. Find dels basd on varType and ref in one allele.
       When it get's here I cannot match the vcf alleles to what
       I have in the mvar fields.
    """
    # d['varType'] == 'del' and 
    m = [buildDelMatches(d, alt) for d in mvarDictList
         if ('' in (d['allele1Seq'], d['allele2Seq'])
             or alt in (d['allele1Seq'], d['allele2Seq'])) and
         d['reference'] in (d['allele1Seq'], d['allele2Seq'])]
    if m:
        return m
    # m = [buildDelMatches(d, alt) for d in mvarDictList
    #         if d['varType'] == 'del' or 
    #            '' in (d['allele1Seq'], d['allele2Seq']) ]
    # if m:
    #     return m
    return [buildDelMatchesFromAlt(d, alt) for d in mvarDictList
            if alt in (d['allele1Seq'], d['allele2Seq']) ]

def lastChance(mvarDictList):
    newMatches = []
    for m in mvarDictList:
        if m['reference'] in (m['allele1Seq'], m['allele2Seq']):
 #           print('ah')
            if m['reference'] == m['allele1Seq']:
                newMatches.append( ('allele2Seq',m) )
            if m['reference'] == m['allele2Seq']:
                newMatches.append( ('allele1Seq',m) )
    return newMatches

def findMatchesByUpdatedVcf(matches, updatedRef, updatedAlt, initRef):
#    print('found it', len(matches))
    if len(matches) == 2:
        if matches[0] == matches[1]:
            return matches[0:1]
    # print('what?')
    # print(matches)
    newMatches = []
    for m in matches:
        mKey, mvar = m
        if mvar['reference'] == updatedRef and mvar[mKey] == updatedAlt:
            newMatches.append(m)
        elif mvar['reference'] == updatedRef[1:] and mvar[mKey] == updatedAlt[1:]:
            newMatches.append(m)

    if len(newMatches) == 2:
        if newMatches[0] == newMatches[1]:
            return newMatches[0:1]

    if len(newMatches) == 2:
        b1 = mvarDictHasBlank(newMatches[0][1])
        b2 = mvarDictHasBlank(newMatches[1][1])
        if b1 and b2:
            # print('bbbbb', newMatches[0][1])
            # print('bbb', newMatches[1][1])
            initRefChars = set()
            for c in initRef:
                initRefChars.add(c)

            oneC, twoC = set(), set()
            for c in newMatches[0][1]['reference']:
                oneC.add(c)
 #           print('aaaa', oneC, newMatches[0][1]['reference'])
            for c in newMatches[1][1]['reference']:
                twoC.add(c)
#            print('aaaa', twoC, newMatches[1][1]['reference'])
            u1 = len(initRefChars & oneC) == len(initRefChars) and len(initRefChars) == len(oneC)
            u2 = len(initRefChars & twoC) == len(initRefChars) and len(initRefChars) == len(twoC)
            if u1 and u2:
                print('help', oneC, twoC, initRefChars)
                i = 1/0
            if u1:
                return newMatches[0:1]
            if u2:
                return newMatches[1:]
            i = 1/0
        if b1:
            return newMatches[0:1]
        if b2:
            return newMatches[1:]

    if not newMatches:
        b1 = mvarDictHasBlank(matches[0][1])
        b2 = mvarDictHasBlank(matches[1][1])
        if b1 and b2:
            # print(matches[0][1])
            # print(matches[1][1])
            # i = 1/0
            return []
        if b1:
            return matches[0:1]
        if b2:
            return matches[1:]

    return newMatches

def findDelMatches(ref, alt, mvarDictList, updatedRef, updatedAlt):
    """Use cgatools ref and alt to locate possible deletions."""
    matches = []
    for mvarDict in mvarDictList:
        if mvarDictHasHomBlank(mvarDict):
            matches += matchHomBlank(ref, alt, mvarDict)
        elif mvarDictHasBlank(mvarDict):
            matches += matchSingleBlank(ref, alt, mvarDict)
        elif mvarDictHasHomDel(mvarDict):
            # 12:849054
            matches += matchHomNonBlank(ref, alt, mvarDict)
        else:
            matches += matchNonBlank(ref, alt, mvarDict)
    if not matches:
#        print('here-debug1', ref, alt, updatedRef, updatedAlt)
        matches += locateByDelAndRef(mvarDictList, alt)
    #if not matches:
        # print('wth')
        # print(matches)
        # print('here-debug2', ref, alt, updatedRef, updatedAlt)
        matches += locateByDelAndRef(mvarDictList, updatedAlt)
    if not matches:
        #print('ok')
        matches += lastChance(mvarDictList)
    if len(matches) > 1:
        #print('here-debug3')
        matches = findMatchesByUpdatedVcf(matches, updatedRef, updatedAlt, ref)
        if len(matches) > 1:
            matches = vcfMvarDepthMatchByDistance.findMatchesByDistance(ref, alt, mvarDictList)
    if len(matches) != 1:        
        print('how is this', mvarDictList)
        matches = vcfMvarDepthMatchByDistance.findMatchesByDistance(ref, alt, mvarDictList)
    if len(matches) != 1:
        print('matches', matches)
        i = 1/0
    return matches[0]
