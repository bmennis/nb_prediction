"""Make hdf5 varinant db for one sample/chrom combo.
   The order of chroms in the calls file matters: 1..22 X M Y
   X:23
   M:24
   Y:25
   Entries must go into the db sorted for the index to work!
"""
import tables, argparse, csv

TRANS = {'X':23, 'M':24, 'Y':25}

def mkNewChrom(chrom):
    if chrom in TRANS:
        return TRANS[chrom]
    return int(chrom)

class Variant(tables.IsDescription):
    chrom = tables.Int32Col()
    pos = tables.Int64Col()
    refDepth = tables.Int32Col()
    altDepth = tables.Int32Col()
    alt = tables.StringCol(30)
    tumorRef = tables.StringCol(30)
    ref = tables.StringCol(30)
    tumorTotDepth = tables.Int32Col()
    totDepth = tables.Int32Col()
    varDepth = tables.Int32Col()
    varFrac = tables.Float64Col()
    tumorAlt = tables.Int32Col()
    noCallSampleCount = tables.Int32Col()
    qualVal = tables.Int32Col()
    varLen = tables.Int32Col()
    tumorFrac = tables.Float64Col()
    Prob = tables.Float64Col()
    isHomVar = tables.Int32Col()
    isVar = tables.Int32Col()
    isNoCall = tables.Int32Col()

def fixBadDel(ref, alt):
    """Less than 10 times, the alt is . and should
       be the first ref char."""
    if alt == '.':
        if len(ref) == 1:
            i = 1/0
        return ref[0]
    return alt

def updateDb(afile, chrpos):
    sample = afile.split('/')[-1].split('.')[0]
    with open(afile) as f:
        reader = csv.DictReader(f, delimiter='\t')
        for row in reader:
            chrpos['chrom'] = mkNewChrom(row['chrom'])
            chrpos['pos'] = int(row['pos'])
            chrpos['refDepth'] = int(row['refDepth'])
            chrpos['altDepth'] = int(row['allDepth'])
            chrpos['alt'] = fixBadDel(row['ref'], row['alt'])
            chrpos['tumorRef'] = row['tumorRef']
            chrpos['ref'] = row['ref']
            chrpos['tumorTotDepth'] = int(row['tumorTotDepth'])
            chrpos['totDepth'] = int(row['totDepth'])
            chrpos['varDepth'] = int(row['varDepth'])
            chrpos['varFrac'] = float(row['varFrac'])
            chrpos['tumorAlt'] = int(row['tumorAlt'])
            chrpos['noCallSampleCount'] = int(float(row['noCallSampleCount']))
            chrpos['qualVal'] = int(row['qualVal'])
            chrpos['varLen'] = int(row['varLen'])
            chrpos['tumorFrac'] = float(row['tumorFrac'])
            chrpos['Prob'] = float(row['Prob'])
            chrpos['isHomVar'] = int(row['isHomVar'])
            chrpos['isVar'] = int(row['isVar'])
            chrpos['isNoCall'] = int(row['isNoCall'])

            chrpos.append()

def main(args):
    h5file = tables.open_file(args.hdfDb, mode = "w", title = "Vars")
    group = h5file.create_group("/", 'posCollection', 'genomic information')
    # expectedrows is from wc on all files. takes a while!
    table = h5file.create_table(group, 'posLs', Variant, "Genomic coord list",
                                expectedrows=int(args.rows))
    chrpos = table.row
    
    with open(args.fileOfVarFiles) as f:
        for afile in f:
            updateDb(afile.strip(), chrpos)

    table.flush()
    table.cols.chrom.create_index()
    table.cols.pos.create_index()
    table.flush()
    h5file.close()

if __name__ == "__main__":
    desc = 'Mk hdf5 db.'
    parser = argparse.ArgumentParser(description=desc)
    argLs = ('rows', 'fileOfVarFiles', 'hdfDb',)
    for param in argLs:
        parser.add_argument(param)
    args = parser.parse_args()
    main(args)
