include: "const.py"


###Pipeline begins with taking output from snpeff and preparing and building gemini from that output 
###from original code /nb_convergence/sf_build_wgs_genimi.py:
"""Split chroms
   vt and normalize
   Build wgs gemini database (one per chrom).
   Do not store samples.
"""

###Begin with annotating snpeff output with dbsift
rule annotateDbnsfp:
    input:  DATA_LOCAL + 'cgiMerge/allSampleCollapseVars.eff.vcf'
    output: DATA_LOCAL + 'cgiMerge/allSampleCollapseVars.eff.dbnsfp.vcf'
    threads: 20
    shell:  """{JAVA} -Xmx32g -Xms16g -jar {SIFT} dbnsfp -v \
               -db {SIFT_DBNSFP} -f {DBNSFP} {input} > {output}"""

###Continue annotation of snpeff output with vcfanno
rule vcfanno:
    input:   vcf = DATA_LOCAL + 'cgiMerge/allSampleCollapseVars.eff.dbnsfp.vcf',
             conf = CONFIG + 'vcfanno.conf',
             lua = VCFANNO_LUA_FILE
    output:  DATA_LOCAL + 'cgiMerge/allSampleCollapseVars.eff.dbnsfp.anno.new.vcf'
    threads: 15
    shell:   """{VCFANNO} -p {threads} \
                -base-path {GEMINI_ANNO} -lua {input.lua} \
                {input.conf} {input.vcf} > {output}"""

###Fix header, need more info on step
HEADER_FIX = 'eff_indel_splice,1,Flag AC,1,Integer AF,1,Float dbNSFP_FATHMM_pred,.,String dbNSFP_Interpro_domain,.,String dbNSFP_LR_pred,.,String dbNSFP_phyloP100way_vertebrate,.,String dbNSFP_phastCons100way_vertebrate,.,String dbNSFP_SIFT_score,.,String dbNSFP_Reliability_index,.,String dbNSFP_RadialSVM_pred,.,String dbNSFP_RadialSVM_pred,.,String dbNSFP_Polyphen2_HVAR_pred,.,String dbNSFP_MutationTaster_pred,.,String dbNSFP_MutationAssessor_pred,.,String'

rule fixHeader:
    input:  DATA_LOCAL + 'cgiMerge/allSampleCollapseVars.eff.dbnsfp.anno.new.vcf'
    output: DATA_LOCAL + 'cgiMerge/allSampleCollapseVars.eff.dbnsfp.anno.new.hHack.vcf'
    shell:  'python {HEADER_HCKR} {input} {output} {HEADER_FIX}'

###Split vcf file by chr for run in gemini
rule limitByChrom:
    input:  DATA_LOCAL + 'cgiMerge/allSampleCollapseVars.eff.dbnsfp.anno.new.hHack.vcf'
    output: DATA_LOCAL + 'cgiMerge/allSampleCollapseVars.eff.dbnsfp.anno.new.hHack.{chrom}.vcf.gz'
    shell:  "grep '^#\\|^{wildcards.chrom}'$'\t' {input} | bgzip -c > {output}"

###Rule for gemini
rule gemini:
    input:  DATA_LOCAL + 'cgiMerge/allSampleCollapseVars.eff.dbnsfp.anno.new.hHack.{chrom}.vcf.gz'
    output: DATA_LOCAL + 'geminiDb/{chrom}.tgtCgiWgs.new.db'
    shell:  '{GPY} {VCFTODB} --legacy-compression {input} {JUNK_PED} {output}'

###Run gemini across all samples and chromosomes
rule all_gemini:
    input: expand(DATA_LOCAL + 'geminiDb/{chrom}.tgtCgiWgs.new.db', \
                  chrom=CHROMS_WITH_Y)

###Next part of pipeline is to query the gemini database
###from original code /nb_convergence/sf_query_cgi_1kg.py:
"""CGI WGS.
   Query gemini and variant db before enrichment.
   Use for splicing variants.
   Just CGI, filtered by global 1kg.
"""

###Query gemini
# kv_is_cgi=0, right? and 1kg needs a none, and exac cov?
# cg54 and esp
rule query_gemini:
    input:  DATA_LOCAL + 'geminiDb/{chrom}.tgtCgiWgs.new.db'
    output: DATA_DISKIN + 'wgsQueryNoPop/geminiQuery/{chrom}'
    run:  
        QB = """query --header -q "select * from variants where cg54_ac is NULL and (kv_is_cgi is NULL or kv_is_cgi=0) """
        Q = QB + "and af_1kg_all<0.01 and af_esp_all<0.01"
        GEMINI_Q = Q + """" """
        shell('{GEMINI} {GEMINI_Q}{input} > {output}')

rule add_exac_vqsr:
    input:  DATA_DISKIN + 'wgsQueryNoPop/geminiQuery/{chrom}',
            '/mnt/isilon/cbmi/variome/bin/gemini/data/gemini_data/ExAC.r0.3.sites.vep.tidy.vcf.gz'
    output: DATA_DISKIN + 'wgsQueryNoPop/geminiQueryFlagExac/{chrom}'
    shell:  'python {OTHER_SCRIPTS}annExacFilterForModel.py {input} {output}'

rule point_one_percent:
    input:  i = DATA_DISKIN + 'wgsQueryNoPop/geminiQueryFlagExac/{chrom}'
    output: o = DATA_DISKIN + 'wgsQueryNoPop/geminiQueryFlagExac1per/{chrom}'
    run:
        apply_point_one_percent(input.i, output.o)

rule rm_regions:
    input:  i = DATA_DISKIN + 'wgsQueryNoPop/geminiQueryFlagExac1per/{chrom}'
    output: o = DATA_DISKIN + 'wgsQueryNoPop/geminiQueryFlagExac1perRmRegion/{chrom}'
    run:
        clean_cgi(input.i, output.o)

rule add_lp:
    """Likely pathogenic. Splice moderate."""
    input:  i = DATA_DISKIN + 'wgsQueryNoPop/geminiQueryFlagExac1perRmRegion/{chrom}'
    output: o = DATA_DISKIN + 'wgsQueryNoPop/geminiQueryFlagExac1perRmRegion_lp/{chrom}'
    run:
        apply_p_or_lp(input.i, output.o, is_lp, 'lp')

# frameshift, stop gain, start lost, canonical splice site, and known hgmd
# examine splice. what is canonical?
# splice high  
rule add_p:
    """Pathogenic"""
    input:  i = DATA_DISKIN + 'wgsQueryNoPop/geminiQueryFlagExac1perRmRegion/{chrom}'
    output: o = DATA_DISKIN + 'wgsQueryNoPop/geminiQueryFlagExac1perRmRegion_p/{chrom}'
    run:
        apply_p_or_lp(input.i, output.o, is_p, 'p')

rule merge_p_lp:
    input:  expand(DATA_DISKIN + 'wgsQueryNoPop/geminiQueryFlagExac1perRmRegion_{x}/{{chrom}}', \
                   x = ('lp', 'p'))
    output: DATA_DISKIN + 'wgsQueryNoPop/geminiQueryFlagExac1perRmRegion_plp/{chrom}'
    run:
        lp, p = list(input)
        shell('cat {p} | rev | cut -f 1 | rev > {output}.tmp')
        shell('paste {lp} {output}.tmp > {output}')

rule update_pathogenic:
    input:  i = DATA_DISKIN + 'wgsQueryNoPop/geminiQueryFlagExac1perRmRegion_plp/{chrom}',
            h = '/home/evansj/me/projects/me/format_hgmd/docs/hgmd.collapse.tab'
    output: o = DATA_DISKIN + 'wgsQueryNoPop/geminiQueryFlagExac1perRmRegion_plp_morep/{chrom}'
    run:
        update_p_and_lp(input.i, input.h, output.o)

rule all_p:
    input: expand(DATA_DISKIN + 'wgsQueryNoPop/geminiQueryFlagExac1perRmRegion_plp_morep/{chrom}', chrom=CHROMS)

# run on resplica
RES_PY = '~/me/respublica/miniconda3/envs/cgi/bin/python'

rule queryHdf_trio:
    input:  DATA_LOCAL + 'hdf5bySample/{sample}.varsWgsCgi.hdf',
            expand(DATA_DISKIN + 'wgsQueryNoPop/geminiQueryFlagExac1perRmRegion_plp_morep/{chrom}', chrom=CHROMS)
    output: DATA_DISKIN + 'wgsQueryNoPop/hdf5QueryAll/{sample}'
    threads: 12
    shell:  '{RES_PY} {SCRIPTS}pull_hdf_mp.py {wildcards.sample} {input} {output}'

rule cat_hdf5:
    input: expand(DATA_DISKIN + 'wgsQueryNoPop/hdf5QueryAll/{sample}.use', \
                  sample=getAllNormalSamplesHaveDb())
    output: DATA_DISKIN + 'wgsQueryNoPop/acc'
    run:
        f = list(input)[0]
        shell('head -1 {f} > {output}')
        for afile in list(input):
            shell('tail -n +2 {afile} >> {output}')

rule finalFilter:
    input:  i = DATA_DISKIN + 'wgsQueryNoPop/acc'
    output: o = DATA_DISKIN + 'wgsQueryNoPop/acc.tab'
    run:
        with open(input.i) as f, open(output.o, 'w') as fout:
            reader = csv.DictReader(f, delimiter='\t')
            fields = reader.fieldnames
            print('\t'.join(fields),
                  file=fout)
            for row in reader:
                # filter any vqsr or this pos
                # correct position for hg19 checks
                row['start'] = str( int(row['start']) + 1 )
                ls = [row[x] for x in fields]
                print('\t'.join(ls), file=fout)

rule applyPop:
    input:  i=DATA_DISKIN + 'wgsQueryNoPop/acc.tab',
            pop='/mnt/isilon/diskin_lab/target_pe/WGS_SAMPLES.tab',
    output: o=DATA_DISKIN + 'wgsQueryNoPop/wgs.tab'
    run:
        pops = pandas.read_csv(input.pop, sep='\t')
        vars = pandas.read_csv(input.i, sep='\t')
        df = pandas.merge(vars, pops, on='sample', how='left')
        df.to_csv(output.o, index=False, sep='\t')

rule uploadNaz:
    input:  DATA_LOCAL + 'wgsQueryByPop/endRmBadExacNazByPop/wgs.queryAccAll_{varType}.{pop}.naz.tab'
    output: LOG + 'boxUpload/{pop}.{varType}.n'
    run: 
        shell('curl -u {USR}:{PASS} -T {input} "https://dav.box.com/dav/Diskin_LAB/5_Projects/aacr_late_2017_data/"')
        shell('touch {output}')


###

