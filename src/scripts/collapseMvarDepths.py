"""Matching mvar and vcf depths sometimes has the same varinat twice.
   I collapse these into one variant.
"""
import argparse, pandas, numpy

def getUniq(x):
    return x.unique()[0]

def main(args):
    df = pandas.read_csv(args.depthFile, delimiter='\t')

    dfg = df.groupby(['chrom', 'pos', 'ref', 'alt']).agg({'chrom':getUniq,
                                                          'pos':getUniq,
                                                          'ref' :getUniq,
                                                          'alt':getUniq,
                                                          'varDepth':numpy.sum,
                                                          'refDepth':max,
                                                          'totDepth':max,
                                                          'allDepth':max,
                                                          'varFrac':numpy.sum,
                                                          'goodDepthFrac':max,
                                                          'hasNoCall':getUniq,
                                                          'mvarZygosityWrong':getUniq,
                                                          'tumorAlt':numpy.sum,
                                                          'tumorRef':max,
                                                          'tumorEffDepth':max,
                                                          'tumorTotDepth':max,
                                                          'qual': getUniq,
                                                          'intervalId':getUniq})

    dfg.to_csv(args.outFile, sep='\t', index=False)

if __name__ == "__main__":
    desc = 'Add var depths for same variants.'
    parser = argparse.ArgumentParser(description=desc)
    argLs = ('depthFile', 'outFile',)
    for param in argLs:
        parser.add_argument(param)
    args = parser.parse_args()
    main(args)
