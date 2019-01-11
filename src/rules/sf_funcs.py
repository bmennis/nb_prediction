


def fixVcfDel(inVcf, outVcf):
    with open(inVcf) as f, open(outVcf, 'w') as fout:
        for line in f:
            if line[0] == '#':
                print(line.strip(), file=fout)
            else:
                sp = line.strip().split('\t')
                if sp[4] == '.':
                    if len(sp[3]) == 1:
                        i = 1/0
                    print('\t'.join(sp[0:4] + [ sp[3][0] ] + sp[5:]),
                          file=fout)
                else:
                    print(line.strip(), file=fout)

def getSamples():
    samples = {}
    pindelList = os.listdir(PINDEL_DIR)
    with open(WG_WX_FILE) as f:
        reader = csv.DictReader(f, delimiter='\t')
        for row in reader:
            ss = '-'.join(row['exomeSample'].split('-')[0:3])
            if not ss in ('TARGET-30-PARKGJ', 'TARGET-30-PAPTAN'): # this is a good sample, but has two normals
                samples[ss] = True
    return set( [x.split('.')[0] for x in pindelList if 'TARGET' == x[0:6] and '-'.join(x.split('.')[0].split('-')[0:3]) in samples] )

def getSampleToName():
    names = defaultdict(dict)
    with open(PATHS) as f:
        for line in f:
            sample, chrom, path = line.strip('\n').split('\t')
            shortSample = '-'.join(sample.split('-')[0:3])
            names[shortSample][sample] = True
    return names

def getAllSamples():
    """Return tumor and normal samples"""
    samples = set()
    with open(FULL_SAMPLE_FILE) as f:
        reader = csv.DictReader(f, delimiter='\t')
        for row in reader:
            if row['platform'] == 'CGI':
                samples.add( row['sample'] )
    return samples

def loadMaleSamples():
    """Return normal and tumor samples with Y. (males)"""
    samples= {}
    with open('/home/evansj/me/projects/diskin/target_meta/working/yStatus') as f:
        for line in f:
            sample, status = line.strip('\n').split('\t')
            if status == 'True':
                samples[sample] = True
    return samples

def loadSamplesByHist(hist):
    samples = set()
    with open(FULL_SAMPLE_FILE) as f:
        reader = csv.DictReader(f, delimiter='\t')
        for row in reader:
            if row['cancer'] == hist:
                samples.add( row['sample'] )
    return samples

def getAllNormalSamples():
    """Return tumor and normal samples"""
    ignoreSamples = ('TARGET-30-PARKGJ-10A-02D',
                     'TARGET-50-PAKJGM-11A-01D',
                     'TARGET-10-PANYYV-10A-01D')
    sample2chrom = []
    with open('/home/evansj/me/projects/diskin/target_meta/working/ntStatus') as f:
        for line in f:
            sample, status = line.strip('\n').split('\t')
            if not sample in ignoreSamples and status == 'normal':
                sample2chrom.append(sample)
    return sample2chrom

def getAllNormalSamplesNbl():
    """Return tumor and normal samples"""
    ignoreSamples = ('PARJXH', 'PASCHP', 'PASYXM', 'PANYYV')
    sample2chrom = []
    nblSamples = loadSamplesByHist('NBL')
    with open('/home/evansj/me/projects/diskin/target_meta/working/ntStatus') as f:
        for line in f:
            sample, status = line.strip('\n').split('\t')
            if not sample in ignoreSamples and status == 'normal' and sample in nblSamples:
                sample2chrom.append(sample)
    return sample2chrom

def rmNoCallAndRefFromVcf(vcfIn, vcfOut):
    """Remove GTs of ./0 or 0/.
       Couldn't get bcftools to do this.
    """
    with open(vcfIn) as f, open(vcfOut, 'w') as fout:
        for line in f:
            if line[0] == '#':
                print(line.strip(), file=fout)
            else:
                gt = line.strip().split('\t')[-1].split(':')[0]
                if not gt in ('./0', '0/.'):
                    print(line.strip(), file=fout)

def getHg19Genes():
    genes = set()
    with open(REF_GENE_HG19_BED) as f:
        for line in f:
            genes.add( line.strip().split('\t')[-1] )
    return genes

def getNazGenes():
    genes = set()
    with open('/home/evansj/me/projects/diskin/target_gene_ls/geneLs/naz.ls') as f:
        for line in f:
            genes.add( line.strip() )
    return genes

def getFavGenes():
    genes = set()
    with open('/home/evansj/me/projects/diskin/target_gene_ls/geneLs/forEnrich') as f:
        for line in f:
            genes.add( line.strip() )
    return genes

def mkGenesForRegression():
    genes = {}
    with open('/home/evansj/me/projects/diskin/nb_convergence/data/queryAccAll_synForBurdenTest') as f:
        reader = csv.DictReader(f, delimiter='\t')
        for row in reader:
            genes[ row['gene'] ] = True
    return set(genes)
