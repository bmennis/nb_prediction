"""Use predictions to make a vcf for this sample and chrom.
   Fill in no-calls.
"""
import argparse, csv

def applyThreshold(homCutoff, varCutoff, predsFile, fout):
    refUpperBound = float(1)-varCutoff
    with open(predsFile) as f:
        reader = csv.DictReader(f, delimiter='\t')
        head = reader.fieldnames
        print('\t'.join(head + ['isHomVar', 'isVar', 'isNoCall']), file=fout)
        for row in reader:
            isHomVar, isVar, isNoCall = 0, 0, 1
            if float(row['Prob']) < refUpperBound:
                isNoCall = 0
            elif float(row['Prob']) > varCutoff:
                isNoCall = 0
                isVar = 1
                if int(row['totDepth']) > 10:
                    if float(row['varFrac']) > homCutoff:
                        isHomVar = 1
                elif int(row['refDepth']) < 2:
                    isHomVar = 1
            vcfLine = '\t'.join([row[x] for x in head] + [str(x) for x in (isHomVar, isVar, isNoCall)])
            print(vcfLine, file=fout)

def main(args):
    with open(args.out, 'w') as fout:
        applyThreshold(float(args.homCutoff), float(args.varCutoff),
                       args.predsFile, fout)

if __name__ == "__main__":
    desc = 'Mk vcf file using predictions.'
    parser = argparse.ArgumentParser(description=desc)
    argLs = ('homCutoff', 'varCutoff', 'chrom',
             'predsFile', 'out')
    for param in argLs:
        parser.add_argument(param)
    args = parser.parse_args()
    main(args)

