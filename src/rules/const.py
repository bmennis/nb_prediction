import os, sys, csv
import pandas as pd
from sfConstInit import *
from sfFuncs import *
from collections import defaultdict

p = os.getcwd()
if 'src' in p:
    PWD = p.split('src/')[0]
else:
    PWD = p + '/'

DATA_LOCAL = PWD + '/data/'
WORK = PWD + 'work/'
FILES = PWD + 'docs/'
SCRIPTS = PWD + '/src/scripts/'
DATA = '/mnt/isilon/diskin_lab/Brian/projects/nb_convergence/data/'
LOG = PWD + 'log/'
CONFIG = PWD + 'configs/'
GEMINI_ANNO = '/mnt/isilon/cbmi/variome/bin/gemini/data/gemini_data'
MODEL = PWD + 'models2/'
OTHER_SCRIPTS = ''

DATA_DISKIN = '/mnt/isilon/diskin_lab/target_pe/nb_convergence/data/'
DATA_QC = '/home/evansj/me/projects/diskin/target_qc/data/'
QUICK_DIR = '/mnt/isilon/diskin_lab/target_pe/nb_data/'
#TARGET = 'nas/nbl3/'

REF_FILE = '/mnt/isilon/maris_lab/target_nbl_ngs/shared_resources/refs/hg19_FA_UCSC/ucsc.hg19_no_chr.fa'
VCFANNO_LUA_FILE = '/mnt/isilon/cbmi/variome/perry/projects/me/vcfanno_lua/scripts/target.lua'
CGI_VCF_HEADER = PWD + 'files/cgiVcfHeader'
ANN_VCF_HEADER = PWD + 'files/vcfHeaderForAnn'
NCI_VCF_HEADER = PWD + 'files/vcfHeaderForNci'
JUNK_PED = PWD + 'files/JUNK_PED.ped'

INDEL_TYPES = ['D', 'TD', 'SI',] # INV
CHROMS = list(range(1,23)) + ['X', 'M']
CHROMS_NO_M = list(range(1,23)) + ['X']
CHROMS_WITH_Y = CHROMS + ['Y']

SAMPLES = getSamples() #['TARGET-30-PAPKXS', 'TARGET-30-PAPTAN', 'TARGET-30-PATFXV']
ALL_SAMPLES = getAllNormalSamples()
ALL_NBL_SAMPLES = getAllNormalSamplesNbl()
ALL_MALE_SAMPLES = loadMaleSamples()
#getSamples()
#['TARGET-30-PAPKXS', 'TARGET-30-PAPTAN', 'TARGET-30-PATFXV'] #getSamples() # #getSamples(
NAMES = getSampleToName()

G_PYTHON = '/mnt/isilon/cbmi/variome/bin/gemini/bin/gemini_python'
#GPY = '/mnt/isilon/cbmi/variome/bin/gemini/anaconda/bin/python'
GPY = '/mnt/isilon/cbmi/variome/bin/gemini/tools/bin/gemini_python'
GEMINI = '/mnt/isilon/cbmi/variome/bin/gemini/tools/bin/gemini'
VCFTODB =  '/mnt/isilon/cbmi/variome/bin/vcf2db/vcf2db.py'

HOM_CUT = '.8'
VAR_CUT = '.7'

