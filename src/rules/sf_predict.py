include: "const.py"

rule predict:
    input:  DATA_LOCAL + 'features/all_mat_sample_chrom/{sample}.{chrom}.mat',
            MODEL + 'snv/limitFeaturesyes.limitDatayes/tree.other.pickle',
            MODEL + 'indel/limitFeaturesyes.limitDatayes/tree.other.pickle'
    output: QUICK_DIR + 'scores/{sample}.{chrom}'
    run:  
        sex = 'F'
        if wildcards.sample in ALL_MALE_SAMPLES:
            sex = 'M'
        shell('python {SCRIPTS}applyTreeNew.py {wildcards.chrom} {sex} {input} {output}')

rule mkCalls:
    input:  QUICK_DIR + 'scores/{sample}.{chrom}'
    output: DATA_LOCAL + 'tmpCalls/{sample}.{chrom}'
    shell:  'python {SCRIPTS}mkSampleChromTab.py {HOM_CUT} {VAR_CUT} {wildcards.chrom} {input} {output}'

def mkAllFinal():
    ls = []
    for sample in ALL_SAMPLES:
        for chrom in mkChroms(sample):
            ls += [DATA_LOCAL + 'tmpCalls/%s.%s' % (sample, chrom)]
    return ls

def mkAllDbsFinal():
    ls = []
    for sample in ALL_SAMPLES:
        for chrom in mkChroms(sample):
            ls += [DATA_LOCAL + 'hdf5/%s.%s.hdf' % (sample, chrom)]
    return ls

rule allCgiScores:
    input:  mkAllFinal()
    output: DATA_LOCAL + 'varsWgsCgi.hdf'
    run: 
        varLs = 'varLs'
        with open(varLs, 'w') as fout:
            for afile in list(input):
                print(afile, file=fout)
        shell('python {SCRIPTS}mkHdfVars.py {varLs} {output}')        

rule mkSampleChromDb:
    input:  DATA_LOCAL + 'tmpCalls/{sample}.{chrom}'
    output: DATA_LOCAL + 'hdf5/{sample}.{chrom}.hdf'
    run:  
        shell('wc -l {input} > {input}.tmp')
        with open(list(input)[0] + '.tmp') as f:
            rows = f.readline().split()[0]
        shell('python {SCRIPTS}mkHdfVarsForSampleChrom.py {rows} {input} {output}')
        shell('rm {input}.tmp')

rule tmp:
    input: mkAllDbsFinal() #DATA_LOCAL + 'hdf5/TARGET-30-PAPKXS-10A-01D.11.hdf'

#rule mkVarDb
