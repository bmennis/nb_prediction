include: "const.py"


rule ls_missing_samples:
"""Rule to print and list the missing samples from analysis"""
    output: o = WORK + 'missing_cgi_no_vcf.tab'
    run:
        with open(output.o, 'w') as fout:
            for sample in ALL_SAMPLES:
                p = DATA_LOCAL + 'cgiMerge/allSampleVcf/{}.vcf'.format(sample)
#p = '/nas/nbl3/masterVarBeta/{}.tsv.bz2'.format(sample)
                if not os.path.exists(p):
                    print(sample, file=fout)


rule gunzipNormVcf:
"""Rule to unzip vcf files"""
    input:  DATA_LOCAL + 'cgi/normalize/{sample}.{chrom}.vcf.gz'
    output: DATA_LOCAL + 'cgi/normalize/{sample}.{chrom}.vcf'
    shell:  'gunzip -c {input} > {output}'

rule mkDepthVcfFull:
"""List normal and tumor depths."""
input:  DATA_LOCAL + 'cgi/normalize/{sample}.{chrom}.vcf'
    output: DATA_LOCAL + 'featuresDepthAllSamples/{sample}.{chrom}'
    log:    LOG + 'featuresFull/{sample}.{chrom}.depths'
    shell:  'python {SCRIPTS}mkDepthFromVcf.py {input} {output} &> {log}'

rule collapseDepthVcf:
    input:  DATA_LOCAL + 'featuresDepthAllSamples/{sample}.{chrom}'
    output: DATA_LOCAL + 'featuresDepthAllSamplesCollapse/{sample}.{chrom}'
    shell:  'python {SCRIPTS}collapseMvarDepths.py {input} {output}'

def getSampleChroms(wc):
# chrom has no chr
    return [ DATA_LOCAL + 'cgi/normalize/%s.%s.vcf' % (wc.sample, chrom) for chrom in mkChroms(wc.sample) ]

# for all pos, get no-calls using list of samples from getSampleChroms
rule catVarsFull:
    input:  getSampleChroms
    output: DATA_LOCAL + 'cgi/allSampleVarBed/{sample}.bed'
    shell:  """cut -f 1,2 {input} | grep -v ^# | uniq | awk '{{print $1 "\\t" $2-1 "\\t" $2}}'| {SORT_BED} > {output}"""

rule collapseFeatsAllSamples:
"""Collapse samples down"""
    input:  DATA_LOCAL + 'cgi/noCallAllCounts',
            DATA_LOCAL + 'featuresDepthAllSamplesCollapse/{sample}.{chrom}',
    output: DATA_LOCAL + 'features/all_mat_sample_chrom/{sample}.{chrom}.mat'
    shell:  'python {SCRIPTS}collapseFeaturesAllSamples.py {input} {output}'

def mkAllFinal():
"""Make list of samples"""
    ls = []
    for sample in ALL_SAMPLES:
        for chrom in mkChroms(sample):
            ls += [DATA_LOCAL + 'features/all_mat_sample_chrom/%s.%s.mat' % (sample, chrom)]
    return ls # [x for x in ls if '.Y' in x]

rule tmp:
    input: DATA_LOCAL + 'features/all_mat_sample_chrom/TARGET-30-PAPKXS-10A-01D.11.mat'

rule allCgi:
    input: mkAllFinal() #expand(, sample=ALL_SAMPLES, chrom=CHROMS  )
