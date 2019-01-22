import sys, os
sys.path.append('/home/evansj/me/projects/me/tool_dirs/')
#from tools import *

FULL_SAMPLE_FILE = '/mnt/isilon/diskin_lab/target_pe/target_meta/working/SAMPLES'
PATHS = '/mnt/isilon/diskin_lab/target_pe/target_meta/working/evidenceDnbs.paths'
PINDEL_DIR = '/mnt/isilon/diskin_lab/zalmanv/pindel_data/BI_WXS_244_SRA_pindel/'
WG_WX_FILE = '/mnt/isilon/diskin_lab/target_pe/recover_trevor/working/wgAndExome.samples.ls'

PWD = os.getcwd().split('code')[0]
WORK = PWD + 'working/'
SCRIPTS = PWD + 'code/scripts/'
QC_SCRIPTS = '/nas/is1/perry/projects/diskin/target_qc/code/'

