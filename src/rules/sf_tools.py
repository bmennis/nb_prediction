import socket

if socket.gethostname() == 'reslnrefo01':
    prefix = '/nas/is1/'
    LUSTRE = '/home/evansj/me/tmp/'
    FQTOOLS = prefix + 'perry/tools/fqtools/bin/fqtools'
    TRIM_GALORE = '/home/evansj/me/bin/trim_galore_zip/trim_galore'
    CUT_ADAPT = 'cutadapt'
else: # respublica
    prefix = '/mnt/isilon/cbmi/variome/'
    LUSTRE = '/mnt/lustre/users/evansj/'
    FQTOOLS = '/home/evansj/me/respublicaTools/anaconda3/envs/py3/bin/fqtools'
    TRIM_GALORE = '~/me/respublicaTools/anaconda3/envs/py27/bin/trim_galore'
    CUT_ADAPT = '~/me/respublicaTools/anaconda3/envs/py27/bin/cutadapt'

HEADER_HCKR = prefix + 'perry/projects/me/vcfHeaderHckr/vcfHeadrHckr.py'
VCFANNO = '/mnt/isilon/cbmi/variome/bin/vcfanno/0.0.11/bin/vcfanno'
VCFTODB = '/mnt/isilon/cbmi/variome/bin/vcf2db_Oct13_2016/vcf2db/vcf2db.py'

GEMINI_ANNO = '/mnt/isilon/cbmi/variome/bin/gemini/data/gemini_data'
GEMINI = '/mnt/isilon/cbmi/variome/bin/gemini/tools/bin/gemini'
GPY = '/mnt/isilon/cbmi/variome/bin/gemini/tools/bin/gemini_python'

GSEA = prefix + 'perry/bin/gsea2-2.2.2.jar'
FASTQC = prefix + 'bin/FastQC/0.10.1/fastqc'
GATK = prefix + 'bin/GenomeAnalysisTK-1.6-2-gc2b74ec/3.3-0/GenomeAnalysisTK.jar'
PICARD = prefix + 'bin/picard-tools-1.121/'
HOMER = prefix + 'bin/homer/bin/'
BED_TOOLS = prefix + 'bin/BEDTools/bedtools2/bin/'
SIFT = '/mnt/isilon/cbmi/variome/perry/tools/snpEff/SnpSift.jar'
SIFT_DBNSFP = '/mnt/isilon/cbmi/variome/perry/tools/snpEff/data/dbNSFP/dbNSFP2.4.txt.gz'
EFF = '/mnt/isilon/cbmi/variome/perry/tools/snpEff/snpEff.jar'
EFF_CONFIG = '/mnt/isilon/cbmi/variome/perry/tools/snpEff/snpEff.config'
SAMTOOLS = prefix + 'bin/samtools1.3/bin/samtools'
BWA = prefix + 'bin/bwa'
KALLISTO = '/home/evansj/me/tools/kallisto_linux-v0.42.4/kallisto'
BCFTOOLS = prefix + 'devkotab/bin/bcftools/bcftools'
GVCFTOOLS = '/home/evansj/me/tools/gvcftools-0.15/bin/'
#SAMTOOLS = '/home/evansj/me/tools/tools2/samtools-0.1.19/samtools'
IGVTOOLS = '/nas/is1/bin/IGVTools/igvtools.jar'

NOVO_SORT = '/nas/is1/bin/Novoalign/3.02.07/novosort'
NOVO_ALIGN = '/nas/is1/bin/Novoalign/3.02.07/novoalign'

PINDEL = prefix + 'perry/tools/pindel/'
SORT_BED = 'sort -k 1,1 -k2,2n'
SORT_VCF = 'sort -k1,1d -k2,2n'

PY27 = '/mnt/isilon/cbmi/variome/perry/miniconda3/envs/p27/bin/python'
JAVA = prefix + 'bin/java'
RSCRIPT = '/usr/bin/Rscript'

MOUSE_FA = prefix + 'reference/mm9.fa'
MOUSE_FA_DICT = prefix + 'reference/mm9.dict'
MOUSE_FA_IDX = prefix + 'reference/mm9.fa.fai'
MM9_SNPS = prefix + 'reference/mm9/mgp.v2.snps.annot.reformat.vcf.gz'
MM9_INDELS = prefix + 'reference/mm9/mgp.v2.indels.annot.reformat.vcf.gz'

HG19_FA = prefix + 'reference/human/hg19/hg19.fa'
HG19_FA_NOCHR = prefix + 'reference/human/hg19/hg19NoChr.fa'
HG19_FA_DICT = prefix + 'reference/human/hg19/hg19.dict'
HG19_FA_IDX = prefix + 'reference/human/hg19/hg19.fa.fai'
HG19_GENOME_FILE_SORTED = prefix + 'perry/data/ucsc/hg19_simple.chrom.sizes'
DBSNP_135 = prefix + 'reference/human/dbsnp135.vcf'
DBSNP_138 = prefix + 'perry/data/human_variation_vcf/00-All_dbsnp138.vcf'
DBSNP_138_GATK_HG19 = prefix + 'perry/data/gatk/dbsnp_138.hg19.coreContig.vcf'
MILLS_HG19 = prefix + 'perry/data/gatk/Mills_and_1000G_gold_standard.indels.hg19.sites.coreContig.vcf'

REF_GENE_HG19 = prefix + 'reference/refGene_refseqhg19.txt'
REF_GENE_HG19_BED = prefix + 'reference/refGene_refseqhg19.bed'

KAVIAR = '/home/evansj/me/projects/diskin/noncoding_nbl_regions/data/kaviar/full/new/Kaviar-160204-Public/vcfs/Kaviar-160204-Public-hg19.vcf.gz'


