"""Add exac freqs to single positions."""
import sys, csv, os
import tabix, queryExAC 

def ann(filteredFile, tb, outFile):
    with open(filteredFile) as f, open(outFile, 'w') as fout:
        reader = csv.DictReader(f, delimiter='\t')
        fields = reader.fieldnames
        print('\t'.join(fields + ['ExACfilter','ExACfilterMatchRefAlt']), file=fout)
        for row in reader:
            chrom = row['chrom']
            pos = int(row['pos'])
            ref, alt = row['ref'], row['alt']
            ref_alt_exacFilter_ls = queryExAC.getPosFilterWithRefAlt(chrom, pos, tb)
            exacFilters = [ r[2] for r in ref_alt_exacFilter_ls ]
            exacFiltersMatchRefAlt = [ r[2] for r in ref_alt_exacFilter_ls if r[0] == ref and r[1] == alt ]
            line = ''
            for field in fields:
                line += row[field] + '\t'
            print(line + ';'.join(exacFilters) + '\t' + ';'.join(exacFiltersMatchRefAlt), file=fout)

if __name__ == '__main__':
    filteredFile, exacFile, outFile = sys.argv[1:]
    tb = tabix.open(exacFile)
    ann(filteredFile, tb, outFile)
