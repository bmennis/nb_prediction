import sys, tables, csv, os
from multiprocessing import Process, Pool
import mkHdfVarsForSample

useFields = ['isHomVar', 'refDepth', 'altDepth', 'totDepth',
             'varFrac', 'tumorRef',
             'tumorAlt', 'tumorTotDepth', 'tumorFrac',
             'noCallSampleCount', 'Prob']

def filterRow(row, ref, alt):
    return row['ref'].decode("utf-8") == ref \
        and row['alt'].decode("utf-8") == alt \
        and row['isVar'] == 1

def fix_field(field):
    if "'" in str(field):
        return field.decode("utf-8")
    return str(field)

def query_hdf(table, useFields, ref, alt, chrom, start):
    return [ [x[uf] for uf in useFields] for x in 
             table.where('(chrom == %d) & (pos == %d)'
                         % (chrom, start))
             if filterRow(x, ref, alt) ]

def process_row(sample, row, table, useFields, fieldNames, fout):
    start = int(row['start']) + 1
    ref = row['ref']
    alt = row['alt']
    rawChrom = row['chrom']
    chrom = mkHdfVarsForSample.mkNewChrom(rawChrom)
    resLs = query_hdf(table, useFields, ref, alt, chrom, start)
    for r in resLs:
        print('\t'.join([ row[x] for x in fieldNames ] + [fix_field(ff) for ff in r] + [sample]),
              file=fout)    

async def load_vars(args, gemini_file, file_handle, fout, field2idx, fieldNames):
#    print(field2idx)
    h5file = tables.open_file(args.hdf_file, mode="r")
    table = h5file.root.posCollection.posLs
    # i = 0
    # fout = 
    # print('here0')
    # with aiofiles.open(gemini_file, mode='r') as f:
    #     print('here1')
        
        
#        async for line in f:
    async for line in file_handle.readlines():
        print('here0')
        sp = line.strip().split('\t')
        if not sp[0] in field2idx:
            row = {field:sp[field2idx[field]] for field in fieldNames}
            await trio.sleep(1) #process_row(row, table, useFields, fout)
        print(gemini_file)
    # i += 1
    # if i == 100:
    #     break

async def main(args):
    gemini_ls = (args.gemini_file_1, args.gemini_file_2, args.gemini_file_3,
                 args.gemini_file_4, args.gemini_file_5)
    files = [open(x) for x in gemini_ls]
    out_files = [open(gemini_file + '.' + args.sample + '.tmp', 'w')
                 for gemini_file in gemini_ls]
    for afile in files:
        fieldNames = afile.readline().strip().split('\t')
        field2idx = {field:x for x,field in enumerate(fieldNames)}
    for afile in out_files:
        print('\t'.join(fieldNames + useFields + ['sample']), file=afile)

    for afile in files:
        afile.close()

    async with trio.open_nursery() as nursery:
        print('one')
        async with aiofiles.open(args.gemini_file_1, mode='r') as f1:
            print('two')
            async with aiofiles.open(args.gemini_file_2, mode='r') as f2:
                async with aiofiles.open(args.gemini_file_3, mode='r') as f3:
                    print('spawning')
                    nursery.spawn(load_vars, args, gemini_file_1, f1, out_files[0], field2idx, fieldNames)
                    nursery.spawn(load_vars, args, gemini_file_2, f2, out_files[1], field2idx, fieldNames)
                    nursery.spawn(load_vars, args, gemini_file_3, f3, out_files[2], field2idx, fieldNames)
        # for out_file, file_handle, gemini_file in zip(out_files, files, gemini_ls):
        #     nursery.spawn(load_vars, args, gemini_file, file_handle, out_file, field2idx, fieldNames)


    for afile in out_files:
        afile.close()

    os.system('head -1 {} > {}'.format(args.gemini_file_1 + '.' + args.sample + '.tmp',
                                       args.outFile))
    for gemini_file in (args.gemini_file_1, args.gemini_file_2, args.gemini_file_3,
                        args.gemini_file_4, args.gemini_file_5):
        os.system('tail -n +2 {} >> {}'.format(args.gemini_file_1 + '.' + args.sample + '.tmp',
                                               args.outFile))

def read_file(params):
    afile_name, out_file_name, name, hdf_file, sample = params
    afile = open(afile_name)
    out_file = open(out_file_name, 'w')
    chrom = name.split('/')[-1].split('.')[0]
    # with open('shit.' + chrom + '.' + sample, 'w') as fout:
    #     print('shit', file=fout)
    # print(out_file)
    # print('test2', file=out_file)
    h5file = tables.open_file(hdf_file, mode="r")
    table = h5file.root.posCollection.posLs

    i = 0
#    stream = curio.io.FileStream(afile)
    sp = afile.readline().strip().split('\t')
    fieldNames = sp
    field2idx = {field:x for x,field in enumerate(fieldNames)}
    print('\t'.join(fieldNames), file=out_file)
#    print('here')
    for line in afile:
        sp = line.strip().split('\t')
        i += 1
        row = {field:sp[field2idx[field]] for field in fieldNames}
        process_row(sample, row, table, useFields, fieldNames, out_file)
        # print(i, name)
        # if i == 1000:
        #     break
    afile.close()
    out_file.close()

# class myThread(threading.Thread):
#     def __init__(self, threadID, open_file, out_file, hdf_file, sample):
#         threading.Thread.__init__(self)
#         self.threadID = threadID
#         self.open_file = open_file
#         self.out_file = out_file
#         self.hdf_file = hdf_file
#         self.sample = sample
#     def run(self):
#         read_file(self.open_file, self.out_file, self.threadID, self.hdf_file, self.sample)

def test(sample, hdf_file, files, out_file):
#    open_files = [open(afile, 'r') for afile in files]
    out_files = [afile + '.' + sample + '.tmp'
                 for afile in files]
    #print('test', file=out_files[0])

    i = 0
    pool = Pool(processes=12)
    input_ls = []
    for open_file, out_afile in zip(files, out_files):
        args = [open_file, out_afile, 't' + str(i), hdf_file, sample]
        input_ls.append(args)
        i += 1
    result = pool.map(read_file, input_ls)
        # i += 1
        # procs.append(proc)
        # proc.start()
    # for t in procs:
    #     t.join()

    # for afile in out_files:
    #     afile.close()

    os.system('head -1 {} > {}'.format(files[0] + '.' + sample + '.tmp',
                                       out_file))
    for gemini_file in files:
        os.system('tail -n +2 {} >> {}'.format(gemini_file + '.' + sample + '.tmp',
                                               out_file))
        # clean up tmp file
        os.system('rm ' + gemini_file + '.' + sample + '.tmp')

if __name__ == "__main__":
    sample, hdf_file = sys.argv[1:3]
    gemini_files = sys.argv[3:-1]
    out_file = sys.argv[-1]

    # desc = 'Pull data for report.'
    # parser = argparse.ArgumentParser(description=desc)
    # argLs = ('sample', 'hdf_file', 'gemini_file_1', 
    #          'gemini_file_2', 'gemini_file_3',
    #          'gemini_file_4', 'gemini_file_5',
    #          'out_file',)
    # for param in argLs:
    #     parser.add_argument(param)
    # args = parser.parse_args()
    test(sample, hdf_file, gemini_files, out_file)
