######################################################

'''
xenofind_pipe.py 

By: J. Sumabat, N. Lai, Y. Huang, S. Peck 

'''
######################################################

import os 
import numpy as np 
import glob
from pathlib import Path 

######################################################

# Consensus Paths
working_dir = '/home/marchandlab/github/jay/capstone/yujia/XenoFind/xenofind_test/240303_sorting_test_2/' #Input desired working/ file output directory here 
raw_dir = '/home/marchandlab/DataAnalysis/Kaplan/raw/fast5/10.4.1/240104_BSn_90mer_xr_train/50fast5' #Input either fast5 or pod5 containing directory here
placeholder_fasta = '/home/marchandlab/github/jay/capstone/reference/xBSn_90mer_fake_randomer.fa' #Input a dummy fasta containing randomer regions and flank barcodes here

# Detection Paths
working_dir = '' #Input desired working/ file output directory here
raw_dir = '' #Input either fast5 or pod5 containing directory here 
con_fasta = '' #Input consensus fasta or canonical ground truth fasta here
######################################################


######################################################

consensus = True 
lq_detect = False #placeholder variables 

######################################################

# Generate consensus fasta 
if consensus == True: 
	cmd = 'python xenofind.py consensus -w ' + working_dir+' -f '+ raw_dir+' -r '+placeholder_fasta
	os.system(cmd)

# Low Quality Score Detection 
if lq_detect == True: 
	cmd = 'python xenofind.py low_qual -w '+working_dir+' -f '+raw_dir+' -r '+con_fasta
	os.system(cmd)

