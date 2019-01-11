include: "const.py"

"""Make one vcf file of variants from all samples.
   This is not a multi sample vcf.
"""

def getSampleChroms(wc):
# chrom has no chr
    return [ DATA_LOCAL + 'cgi/normalize/%s.%s.vcf' % (wc.sample, chrom) for chrom in mkChroms(wc.sample) ]

def mkChroms(sample):
    if sample in ALL_MALE_SAMPLES:
        return ['Y'] + CHROMS
    return CHROMS

rule catSampleChromVcfFull:
    input:  getSampleChroms
    output: DATA_LOCAL + 'cgiMerge/allSampleVcf/{sample}.vcf'
    shell:  """cut -f 1-5 {input} | grep -v ^# | uniq | awk '{{print $0 "\\t.\\t.\\t.\\tGT\\t0/1" }}' | {SORT_VCF} -T /home/evansj/me/tmp/ > {output}"""

rule collapseSampleVcfs:
    input:  expand(DATA_LOCAL + 'cgiMerge/allSampleVcf/{sample}.vcf', sample=ALL_SAMPLES)
    output: DATA_LOCAL + 'cgiMerge/allSampleCollapseVars.vcf'
    run:  
        shell('cat {ANN_VCF_HEADER} > {output}')
        shell('{SORT_BED} -T /home/evansj/me/tmp/ -m {input} | uniq >> {output}')
