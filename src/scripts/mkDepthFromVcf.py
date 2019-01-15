"""Use only vcf file, make depths."""
import argparse, sys
from vcfDepthHelpers import *

def unPhasedAndNoMissing(altCounts, call):
    return not '.' in call and not '.' in altCounts and '/' in call

def missingCall(altCounts):
    """Missing data in allele1 and allele2 counts."""
    return '.' in altCounts

def mkDepthLine(sp):
    depthData = sp[-1]
    totDepth, altCounts, refDepth, mvarData = depthData.split(':')[-4:]
    qual = depthData.split(':')[2]
    call = depthData.split(':')[0]
    printLs = ['junk']
    if not 'OLD_MULTIALLELIC' in sp[-3] and (call == '1/1' or call == '1|1'):
        printLs = printHomAlt(sp, totDepth, altCounts, refDepth, mvarData)
    elif unPhasedAndNoMissing(altCounts, call):
        if not 'OLD_MULTIALLELIC' in sp[-3]:
            #print('\t'.join(sp), file=sys.stderr)
            printLs = printRefAndAlt(sp, totDepth, altCounts, refDepth, mvarData)
        elif 'OLD_MULTIALLELIC' in sp[-3]:
            #print('\t'.join(sp), file=sys.stderr)
            i = 1/0 #printMultiAlt(sp, totDepth, altCounts, refDepth, mvarData, fout)
        else:
            i = 1/0
    elif missingCall(altCounts):
        printLs = printMissingCall(sp, totDepth, altCounts, refDepth, mvarData)
    elif '|' in call: # phased
        if call in ('1|0', '0|1'):
            # one is ref
            printLs = printRefAndAlt(sp, totDepth, altCounts, refDepth, mvarData)
        else:
            printLs = printMultiWithFakeMissing(sp, totDepth, altCounts, refDepth, mvarData)
    elif '/' in call and '.' in call:
        #print('\t'.join(sp), file=sys.stderr)
        printLs = printMultiWithFakeMissing(sp, totDepth, altCounts, refDepth, mvarData)
    else: 
         print('debug', sp)
         i = 1/0
    intervalId = sp[-1].split('/')[-1]
    return printLs + [qual, intervalId]

    #print('\t'.join(printLs) + '\t' + qual, file=fout)

def printDepths(f, fout):
    print('chrom\tpos\tref\talt\tvarDepth\trefDepth\ttotDepth\tallDepth\tvarFrac\tgoodDepthFrac\thasNoCall\tmvarZygosityWrong\ttumorAlt\ttumorRef\ttumorEffDepth\ttumorTotDepth\tqual\tintervalId', file=fout)
    for line in f:
        if line[0] != '#':
            print('\t'.join( mkDepthLine(line.split()) ), file=fout)

def main(args):
    with open(args.vcfIn) as f, open(args.outFile, 'w') as fout:
        printDepths(f, fout)

if __name__ == '__main__':
    desc = 'Make depths for decomposed vcf file.'
    parser = argparse.ArgumentParser(description=desc)
    for param in ('vcfIn', 'outFile'):
        parser.add_argument(param)
    args = parser.parse_args()
    main(args)

