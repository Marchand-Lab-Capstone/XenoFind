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
working_dir = '/home/marchandlab/github/jay/capstone/XenoFind/xenofind_test/240313_labelled_headers/' #Input desired working/ file output directory here 
raw_dir = '/home/marchandlab/DataAnalysis/Kaplan/raw/fast5/10.4.1/240104_BSn_90mer_xr_train/50fast5' #Input either fast5 or pod5 containing directory here
#raw_dir = '/home/marchandlab/DataAnalysis/Sumabat/230725_PZ_lib_v4_r10/20230725_1220_MN37138_APH167_a204cb54/subset_fast5' #Input either fast5 or pod5 containing directory here
placeholder_fasta = '/home/marchandlab/github/jay/capstone/reference/xBSn_90mer_fake_randomer.fa' #Input a dummy fasta containing randomer regions and flank barcodes here

# Detection Paths
#working_dir_detect = '/home/marchandlab/github/jay/capstone/XenoFind/xenofind_test/240311_lq_test/' #Input desired working/ file output directory here
#raw_dir_detect = '/home/marchandlab/DataAnalysis/Kaplan/raw/fast5/10.4.1/240104_BSn_90mer_xr_train/50fast5' #Input either fast5 or pod5 containing directory here 
#con_fasta = '/home/marchandlab/github/jay/capstone/XenoFind/xenofind_test/240309_cluster_collapse/outputs/polished_consensus.fasta' #Input consensus fasta or canonical ground truth fasta here

working_dir_detect = '/Users/hyj/Curriculum/Capstone_Project/xenofind_test/240312_lq_test_1fast5'
raw_dir_detect = '/Users/hyj/Curriculum/Capstone_Project/datasets/240104_BSn_90mer_xr_train_capstone_set/1fast5'
con_fasta = '/Users/hyj/Curriculum/Capstone_Project/consensus.fasta'

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
	cmd = 'python xenofind.py low_qual -w '+working_dir_detect+' -f '+raw_dir_detect+' -r '+con_fasta
	os.system(cmd)

