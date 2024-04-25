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
'''
DATASETS TO TESTS
/home/marchandlab/DataAnalysis/Kaplan/raw/pod5/240215_PZ_NB25_xr_Train <-- single sequence context, PZ, length of reads ~ 139 bp
/home/marchandlab/DataAnalysis/Sumabat/230725_PZ_lib_v4_r10/20230725_1220_MN37138_APH167_a204cb54 <---- all sequence context, PZ, length??? 

'''
# Consensus Paths
working_dir = '/home/marchandlab/github/jay/capstone/XenoFind/xenofind_test/240416_medaka_iteration_tests/' #Input desired working/ file output directory here  <--- fix directory function to not neeed '/' at end
raw_dir = '/home/marchandlab/DataAnalysis/Kaplan/raw/fast5/10.4.1/240104_BSn_90mer_xr_train/50fast5' #Input either fast5 or pod5 containing directory here
#raw_dir = '/home/marchandlab/DataAnalysis/Sumabat/230725_PZ_lib_v4_r10/20230725_1220_MN37138_APH167_a204cb54/subset_fast5' #Input either fast5 or pod5 containing directory here
#raw_dir = '/home/marchandlab/DataAnalysis/Kaplan/raw/pod5/240215_PZ_NB25_xr_Train/20240215_1810_MN37138_ARS988_4bbd5246/pod5_0-72'
#raw_dir = '/home/marchandlab/DataAnalysis/Sumabat/230725_PZ_lib_v4_r10/20230725_1220_MN37138_APH167_a204cb54/fast5'
placeholder_fasta = '/home/marchandlab/github/jay/capstone/reference/xBSn_90mer_fake_randomer.fa' #Input a dummy fasta containing randomer regions and flank barcodes here
#placeholder_fasta = '/home/marchandlab/github/jay/capstone/reference/xPZ_NB25_xr_Train_fake_randomer.fasta' #PZ 139 mer data set
#placeholder_fasta = '/home/marchandlab/github/jay/capstone/reference/xxref_libv2_PZ_CxDx-_fake_randomer.fasta'#PZ xenomorph library

# Detection Paths
working_dir_detect = '/home/marchandlab/github/jay/capstone/XenoFind/xenofind_test/240412_BSn_90mer_find_tests' #Input desired working/ file output directory here
raw_dir_detect = '/home/marchandlab/DataAnalysis/Kaplan/raw/fast5/10.4.1/240104_BSn_90mer_xr_train/50fast5' #Input either fast5 or pod5 containing directory here 
con_fasta = '/home/marchandlab/github/jay/capstone/XenoFind/xenofind_test/240309_cluster_collapse/outputs/polished_consensus.fasta' #Input consensus fasta or canonical ground truth fasta here


######################################################


######################################################

consensus = True
find = False #placeholder variables 

######################################################

# Generate consensus fasta 
if consensus == True: 
	cmd = 'python xenofind.py consensus -w ' + working_dir+' -f '+ raw_dir+' -r '+placeholder_fasta
	os.system(cmd)

# Low Quality Score Detection 
if find == True: 
	cmd = 'python xenofind.py find -w '+working_dir_detect+' -f '+raw_dir_detect+' -r '+con_fasta
	os.system(cmd)

