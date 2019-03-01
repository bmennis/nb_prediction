include: "const.py"
include: "sfFuncs.py"
include: "sf_tools.py"
sys.path.append(SCRIPTS)
from var_classes import *

###Pipeline begins with taking output from snpeff and preparing and building gemini from that output 
###from original code /nb_convergence/sf_build_wgs_genimi.py:
"""Split chroms
   vt and normalize
   Build wgs gemini database (one per chrom).
   Do not store samples.
"""

DBNSFP = 'Interpro_domain,SIFT_score,Polyphen2_HVAR_pred,RadialSVM_pred,LR_pred,Reliability_index,FATHMM_pred,MutationAssessor_pred,MutationTaster_pred,phyloP100way_vertebrate,phastCons100way_vertebrate'

###Begin with annotating snpeff output with dbsift
rule annotateDbnsfp:
    input:  DATA + 'cgiMerge/allSampleCollapseVars.eff.vcf'
    output: DATA + 'cgiMerge/allSampleCollapseVars.eff.dbnsfp.vcf'
    threads: 20
    shell:  """{JAVA} -Xmx32g -Xms16g -jar {SIFT} dbnsfp -v \
               -db {SIFT_DBNSFP} -f {DBNSFP} {input} > {output}"""

###Continue annotation of snpeff output with vcfanno
rule vcfanno:
    input:   vcf = DATA + 'cgiMerge/allSampleCollapseVars.eff.dbnsfp.vcf',
             conf = CONFIG + 'vcfanno.conf',
             lua = VCFANNO_LUA_FILE
    output:  DATA + 'cgiMerge/allSampleCollapseVars.eff.dbnsfp.anno.new.vcf'
    threads: 15
    shell:   """{VCFANNO} -p {threads} \
                -base-path {GEMINI_ANNO} -lua {input.lua} \
                {input.conf} {input.vcf} > {output}"""

###Fix header, need more info on step what does {HEADER_HCKR} refer to?
HEADER_FIX = 'eff_indel_splice,1,Flag AC,1,Integer AF,1,Float dbNSFP_FATHMM_pred,.,String dbNSFP_Interpro_domain,.,String dbNSFP_LR_pred,.,String dbNSFP_phyloP100way_vertebrate,.,String dbNSFP_phastCons100way_vertebrate,.,String dbNSFP_SIFT_score,.,String dbNSFP_Reliability_index,.,String dbNSFP_RadialSVM_pred,.,String dbNSFP_RadialSVM_pred,.,String dbNSFP_Polyphen2_HVAR_pred,.,String dbNSFP_MutationTaster_pred,.,String dbNSFP_MutationAssessor_pred,.,String'

rule fixHeader:
    input:  DATA + 'cgiMerge/allSampleCollapseVars.eff.dbnsfp.anno.new.vcf'
    output: DATA + 'cgiMerge/allSampleCollapseVars.eff.dbnsfp.anno.new.hHack.vcf'
    shell:  'python {HEADER_HCKR} {input} {output} {HEADER_FIX}'

###Split vcf file by chr for run in gemini
rule limitByChrom:
    input:  DATA + 'cgiMerge/allSampleCollapseVars.eff.dbnsfp.anno.new.hHack.vcf'
    output: DATA + 'cgiMerge/allSampleCollapseVars.eff.dbnsfp.anno.new.hHack.{chrom}.vcf.gz'
    shell:  "grep '^#\\|^{wildcards.chrom}'$'\t' {input} | bgzip -c > {output}"

###Rule for gemini
rule gemini:
    input:  DATA + 'cgiMerge/ZVmod.allSampleCollapseVars.eff.dbnsfp.anno.new.hHack.{chrom}.vcf.gz'
    output: DATA + 'geminiDb/{chrom}.ZVmod.tgtCgiWgs.new.db'
    conda: 'envs/vcf2db'
    shell:  '{GPY} {VCFTODB} --legacy-compression {input} {JUNK_PED} {output}'

###Run gemini across all samples and chromosomes
rule all_gemini:
    input: expand(DATA + 'geminiDb/{chrom}.ZVmod.tgtCgiWgs.new.db', \
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
    input:  DATA + 'geminiDb/{chrom}.ZVmod.tgtCgiWgs.new.db'
    output: DATA + 'wgsQueryNoPop/geminiQuery/{chrom}.ZVmod'
    run:  
        QB = """query --header -q "select * from variants where cg54_ac is NULL and (kv_is_cgi is NULL or kv_is_cgi=0) """
        Q = QB + "and af_1kg_all<0.01 and af_esp_all<0.01"
        GEMINI_Q = Q + """" """
        shell('{GEMINI} {GEMINI_Q}{input} > {output}')

rule add_exac_vqsr:
    input:  DATA + 'wgsQueryNoPop/geminiQuery/{chrom}.ZVmod',
            '/mnt/isilon/cbmi/variome/bin/gemini/data/gemini_data/ExAC.r0.3.sites.vep.tidy.vcf.gz'
    output: DATA + 'wgsQueryNoPop/geminiQueryFlagExac/{chrom}.ZVmod'
    shell:  'python {SCRIPTS}annExacFilter.py {input} {output}'

###This rule is in reference to frequency and filters on that frequency
rule point_one_percent:
    input:  i = DATA + 'wgsQueryNoPop/geminiQueryFlagExac/{chrom}.ZVmod'
    output: o = DATA + 'wgsQueryNoPop/geminiQueryFlagExac1per/{chrom}.ZVmod'
    run:
        apply_point_one_percent(input.i, output.o)

###Cross reference vqsr and point one percent outputs
rule rm_regions:
    input:  i = DATA + 'wgsQueryNoPop/geminiQueryFlagExac1per/{chrom}.ZVmod'
    output: o = DATA + 'wgsQueryNoPop/geminiQueryFlagExac1perRmRegion/{chrom}.ZVmod'
    run:
        clean_cgi(input.i, output.o)

rule add_lp:
    """Likely pathogenic. Splice moderate."""
    input:  i = DATA + 'wgsQueryNoPop/geminiQueryFlagExac1perRmRegion/{chrom}.ZVmod'
    output: o = DATA + 'wgsQueryNoPop/geminiQueryFlagExac1perRmRegion_lp/{chrom}.ZVmod'
    run:
        apply_p_or_lp(input.i, output.o, is_lp, 'lp')

# frameshift, stop gain, start lost, canonical splice site, and known hgmd
# examine splice. what is canonical?
# splice high  
rule add_p:
    """Pathogenic"""
    input:  i = DATA + 'wgsQueryNoPop/geminiQueryFlagExac1perRmRegion/{chrom}.ZVmod'
    output: o = DATA + 'wgsQueryNoPop/geminiQueryFlagExac1perRmRegion_p/{chrom}.ZVmod'
    run:
        apply_p_or_lp(input.i, output.o, is_p, 'p')

###Rule to merge pathogenic and likely pathogenic
rule merge_p_lp:
    input:  expand(DATA + 'wgsQueryNoPop/geminiQueryFlagExac1perRmRegion_{x}/{{chrom}}.ZVmod', \
                   x = ('lp', 'p'))
    output: DATA + 'wgsQueryNoPop/geminiQueryFlagExac1perRmRegion_plp/{chrom}.ZVmod'
    run:
        lp, p = list(input)
        shell('cat {p} | rev | cut -f 1 | rev > {output}.tmp')
        shell('paste {lp} {output}.tmp > {output}')

rule update_pathogenic:
    input:  i = DATA + 'wgsQueryNoPop/geminiQueryFlagExac1perRmRegion_plp/{chrom}.ZVmod',
            h = '/mnt/isilon/cbmi/variome/perry/projects/me/format_hgmd/docs/hgmd.collapse.tab'
    output: o = DATA + 'wgsQueryNoPop/geminiQueryFlagExac1perRmRegion_plp_morep/{chrom}.ZVmod'
    run:
        update_p_and_lp(input.i, input.h, output.o)

rule all_p:
    input: expand(DATA + 'wgsQueryNoPop/geminiQueryFlagExac1perRmRegion_plp_morep/{chrom}.ZVmod', chrom=CHROMS)

# run on resplica
RES_PY = '~/me/respublica/miniconda3/envs/cgi/bin/python'

rule queryHdf_trio:
    input:  DATA + 'hdf5bySample/{sample}.varsWgsCgi.hdf',
            expand(DATA + 'wgsQueryNoPop/geminiQueryFlagExac1perRmRegion_plp_morep/{chrom}.ZVmod', chrom=CHROMS)
    output: DATA_LOCAL + 'wgsQueryNoPop/hdf5QueryAll/{sample}'
    threads: 12
    shell:  '{RES_PY} {SCRIPTS}pull_hdf_mp.py {wildcards.sample} {input} {output}'

rule cat_hdf5:
    input: expand(DATA + 'wgsQueryNoPop/hdf5QueryAll/{sample}.use', \
                  sample=getAllNormalSamplesHaveDb())
    output: DATA + 'wgsQueryNoPop/acc'
    run:
        f = list(input)[0]
        shell('head -1 {f} > {output}')
        for afile in list(input):
            shell('tail -n +2 {afile} >> {output}')

rule finalFilter:
    input:  i = DATA + 'wgsQueryNoPop/acc'
    output: o = DATA + 'wgsQueryNoPop/acc.tab'
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
    input:  i= DATA + 'wgsQueryNoPop/acc.tab',
            pop='/mnt/isilon/diskin_lab/target_pe/WGS_SAMPLES.tab',
    output: o= DATA + 'wgsQueryNoPop/wgs.tab'
    run:
        pops = pandas.read_csv(input.pop, sep='\t')
        vars = pandas.read_csv(input.i, sep='\t')
        df = pandas.merge(vars, pops, on='sample', how='left')
        df.to_csv(output.o, index=False, sep='\t')

rule uploadNaz:
    input:  DATA + 'wgsQueryByPop/endRmBadExacNazByPop/wgs.queryAccAll_{varType}.{pop}.naz.tab'
    output: LOG + 'boxUpload/{pop}.{varType}.n'
    run: 
        shell('curl -u {USR}:{PASS} -T {input} "https://dav.box.com/dav/Diskin_LAB/5_Projects/aacr_late_2017_data/"')
        shell('touch {output}')


###The query actually is the final step of pipeline to get results of the predictions


###Including the matrix generation snakefile for pipeline, not entirely sure where this fits in the order
###from original code file nb_convergence/sfCgiFull.py

rule ls_missing_samples:
    output: o = DATA + 'missing_cgi_no_vcf.tab'
    run:
        with open(output.o, 'w') as fout:
            for sample in ALL_SAMPLES:
                p = DATA + 'cgiMerge/allSampleVcf/{}.vcf'.format(sample)
                #p = '/nas/nbl3/masterVarBeta/{}.tsv.bz2'.format(sample)
                if not os.path.exists(p):
                    print(sample, file=fout)

#Unzip normalized vcf files of samples and chroms
rule gunzipNormVcf:
    input:  DATA_DISKIN + 'cgi/normalize/{sample}.{chrom}.vcf.gz'
    output: DATA + 'cgi/normalize/{sample}.{chrom}.vcf'
    shell:  'gunzip -c {input} > {output}'

#Downsample the vcf sample files to test the pipeline
rule sampleTargetSample:
    input: DATA + 'cgi/normalize/{sample}.{chrom}.vcf'
    output: DATA + 'cgi/normalizeSample/sample.{sample}.{chrom}.vcf'
    shell: 'head -n 1000 {input} > {output}'

#Get depths from vcf files
rule mkDepthVcfFull:
    """List normal and tumor depths."""
    input:  DATA + 'cgi/normalizeSample/sample.{sample}.{chrom}.vcf'
    output: DATA + 'featuresDepthAllSamplesSample/sample.{sample}.{chrom}'
    log:    LOG + 'featuresFull/sample.{sample}.{chrom}.depths'
    shell:  'python {SCRIPTS}mkDepthFromVcf.py {input} {output} &> {log}'

rule collapseDepthVcf:
    input:  DATA + 'featuresDepthAllSamplesSample/sample.{sample}.{chrom}'
    output: DATA + 'featuresDepthAllSamplesCollapseSample/sample.{sample}.{chrom}'
    shell:  'python {SCRIPTS}collapseMvarDepths.py {input} {output}'

def getSampleChroms(wc):
    # chrom has no chr
    return [ DATA + 'cgi/normalize/%s.%s.vcf' % (wc.sample, chrom)
             for chrom in mkChroms(wc.sample) ]

# for all pos, get no-calls
rule catVarsFull:
    input:  getSampleChroms
    output: DATA + 'cgi/allSampleVarBed/{sample}.bed'
    shell:  """cut -f 1,2 {input} | grep -v ^# | uniq | awk '{{print $1 "\\t" $2-1 "\\t" $2}}'| {SORT_BED} > {output}"""

rule collapseFeatsAllSamples:
    input:  DATA_DISKIN + 'cgi/noCallAllCounts',
            DATA + 'featuresDepthAllSamplesCollapseSample/sample.{sample}.{chrom}',
    output: DATA + 'features/all_mat_sample_chrom_sample/sample.{sample}.{chrom}.mat'
    shell:  'python {SCRIPTS}collapseFeaturesAllSamples.py {input} {output}'

def mkAllFinal():
    ls = []
    for sample in ALL_SAMPLES:
        for chrom in mkChroms(sample):
            ls += [DATA + 'features/all_mat_sample_chrom_sample/sample.%s.%s.mat' % (sample, chrom)]
    return ls # [x for x in ls if '.Y' in x]

rule tmp:
    input: DATA + 'features/all_mat_sample_chrom/TARGET-30-PAPKXS-10A-01D.11.mat'

rule allCgi:
    input: mkAllFinal() #expand(, sample=ALL_SAMPLES, chrom=CHROMS  )


###Also including the prediction rules for the pipeline
###from original code nb_convergence/sfPredict.py:
"""Apply the models."""

#Run decision tree prediction on snv and indel pickled trained models
rule samp_mat:
    input: DATA_DISKIN + 'features/all_mat_sample_chrom/{sample}.{chrom}.mat'
    output: DATA + 'features/all_mat_sample_chrom_sample/{sample}.{chrom}.samp.mat'
    run: 
        shell('head -n 1 {input} > {output}')
        shell('tail -n +2 {input} | shuf -n 1000 >> {output}')

rule predict:
    input:  DATA + 'features/all_mat_sample_chrom_sample/{sample}.{chrom}.samp.mat',
            MODEL + 'snv/limitFeaturesyes.limitDatayes/tree.other.pickle',
            MODEL + 'indel/limitFeaturesyes.limitDatayes/tree.other.pickle'
    output: QUICK_DIR + 'scores/{sample}.{chrom}.samp'
    run:  
        sex = 'F'
        if wildcards.sample in ALL_MALE_SAMPLES:
            sex = 'M'
        shell('python {SCRIPTS}applyTreeNew.py {wildcards.chrom} {sex} {input} {output}')

rule mkCalls:
    input:  QUICK_DIR + 'scores/{sample}.{chrom}.samp'
    output: DATA + 'tmpCalls/{sample}.{chrom}.samp'
    shell:  'python {SCRIPTS}mkSampleChromTab.py {HOM_CUT} {VAR_CUT} {wildcards.chrom} {input} {output}'

def mkAllFinal():
    ls = []
    for sample in ALL_SAMPLES:
        for chrom in mkChroms(sample):
            ls += [DATA + 'tmpCalls/%s.%s.samp' % (sample, chrom)]
    return ls

def mkAllDbsFinal():
    ls = []
    for sample in ALL_SAMPLES:
        for chrom in mkChroms(sample):
            ls += [DATA + 'hdf5/%s.%s.samp.hdf' % (sample, chrom)]
    return ls


rule allCgiScores:
    input:  mkAllFinal()
    output: DATA + 'varsWgsCgi.hdf'
    run: 
        varLs = 'varLs'
        with open(varLs, 'w') as fout:
            for afile in list(input):
                print(afile, file=fout)
        shell('python {SCRIPTS}mkHdfVars.py {varLs} {output}')        

#Rule to make hdf5 database for sample and chromosome variant calls
rule mkSampleChromDb:
    input:  DATA + 'tmpCalls/{sample}.{chrom}.samp'
    output: DATA + 'hdf5/{sample}.{chrom}.samp.hdf'
    run:  
        shell('wc -l {input} > {input}.tmp')
        with open(list(input)[0] + '.tmp') as f:
            rows = f.readline().split()[0]
        shell('python {SCRIPTS}mkHdfVarsForSampleChrom.py {rows} {input} {output}')
        shell('rm {input}.tmp')

rule tmp2:
    input: mkAllDbsFinal() #DATA_LOCAL + 'hdf5/TARGET-30-PAPKXS-10A-01D.11.hdf'

#rule mkVarDb
