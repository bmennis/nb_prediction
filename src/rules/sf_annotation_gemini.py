include: "const.py"

"""Split chroms
   vt and normalize
   Build wgs gemini database (one per chrom).
   Do not store samples.
"""

DBNSFP = 'Interpro_domain,SIFT_score,Polyphen2_HVAR_pred,RadialSVM_pred,LR_pred,Reliability_index,FATHMM_pred,MutationAssessor_pred,MutationTaster_pred,phyloP100way_vertebrate,phastCons100way_vertebrate'

rule fixVcfDelRule:
"""11 entries have . instead of the ref allele for dels"""
    input:  DATA_LOCAL + 'cgiMerge/allSampleCollapseVars.vcf'
    output: DATA_LOCAL + 'cgiMerge/allSampleCollapseVarsDelFix.vcf'
    run:
        fixVcfDel( list(input)[0], list(output)[0] )

rule snpeff:
    input:  DATA_LOCAL + 'cgiMerge/allSampleCollapseVarsDelFix.vcf'
    output: DATA_LOCAL + 'cgiMerge/allSampleCollapseVars.eff.vcf'
    threads: 20
    shell:  """{JAVA} -Xmx32g -Xms16g -jar {EFF} eff \
               -strict -noStats hg19 -c {EFF_CONFIG} \
               {input} > {output}"""

rule annotateDbnsfp:
    input:  DATA_LOCAL + 'cgiMerge/allSampleCollapseVars.eff.vcf'
    output: DATA_LOCAL + 'cgiMerge/allSampleCollapseVars.eff.dbnsfp.vcf'
    threads: 20
    shell:  """{JAVA} -Xmx32g -Xms16g -jar {SIFT} dbnsfp -v \
               -db {SIFT_DBNSFP} -f {DBNSFP} {input} > {output}"""

rule vcfanno:
    input:   vcf = DATA_LOCAL + 'cgiMerge/allSampleCollapseVars.eff.dbnsfp.vcf',
             conf = CONFIG + 'vcfanno.conf',
             lua = VCFANNO_LUA_FILE
    output:  DATA_LOCAL + 'cgiMerge/allSampleCollapseVars.eff.dbnsfp.anno.new.vcf'
    threads: 15
    shell:   """{VCFANNO} -p {threads} \
                -base-path {GEMINI_ANNO} -lua {input.lua} \
                {input.conf} {input.vcf} > {output}"""

HEADER_FIX = 'eff_indel_splice,1,Flag AC,1,Integer AF,1,Float dbNSFP_FATHMM_pred,.,String dbNSFP_Interpro_domain,.,String dbNSFP_LR_pred,.,String dbNSFP_phyloP100way_vertebrate,.,String dbNSFP_phastCons100way_vertebrate,.,String dbNSFP_SIFT_score,.,String dbNSFP_Reliability_index,.,String dbNSFP_RadialSVM_pred,.,String dbNSFP_RadialSVM_pred,.,String dbNSFP_Polyphen2_HVAR_pred,.,String dbNSFP_MutationTaster_pred,.,String dbNSFP_MutationAssessor_pred,.,String'

rule fixHeader:
    input:  DATA_LOCAL + 'cgiMerge/allSampleCollapseVars.eff.dbnsfp.anno.new.vcf'
    output: DATA_LOCAL + 'cgiMerge/allSampleCollapseVars.eff.dbnsfp.anno.new.hHack.vcf'
    shell:  'python {HEADER_HCKR} {input} {output} {HEADER_FIX}'

rule limitByChrom:
    input:  DATA_LOCAL + 'cgiMerge/allSampleCollapseVars.eff.dbnsfp.anno.new.hHack.vcf'
    output: DATA_LOCAL + 'cgiMerge/allSampleCollapseVars.eff.dbnsfp.anno.new.hHack.{chrom}.vcf.gz'
    shell:  "grep '^#\\|^{wildcards.chrom}'$'\t' {input} | bgzip -c > {output}"

rule gemini:
    input:  DATA_LOCAL + 'cgiMerge/allSampleCollapseVars.eff.dbnsfp.anno.new.hHack.{chrom}.vcf.gz'
    output: DATA_LOCAL + 'geminiDb/{chrom}.tgtCgiWgs.new.db'
    shell:  '{GPY} {VCFTODB} --legacy-compression {input} {JUNK_PED} {output}'

rule all_gemini:
    input: expand(DATA_LOCAL + 'geminiDb/{chrom}.tgtCgiWgs.new.db', \
                  chrom=CHROMS_WITH_Y)
