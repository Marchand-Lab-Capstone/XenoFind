########################################################################
########################################################################
"""
xr_params.py 

Title: Unpublished work

By: H. Kawabe, N. Kaplan, J. Sumabat, J. A. Marchand

Updated: 11/28/23
"""
########################################################################

############################################################
##Analysis instructions 

#Re-basecall pod5 file. Required if new reference files are being used. 
basecall_pod = False

#Perform Quality Score Analysis 
analyze_fastq = True
############################################################



############################################################
#Guppy Base caller configuration

#Path to guppy basecaller
#basecaller_path = '~/dorado-0.5.3-linux-x64/bin/dorado' #update to latest version tesing right now 
basecaller_path = '~/dorado-0.5.3-osx-arm64/bin/dorado'

#GPU enabled 
device_type = 'cuda:all' 

#Guppy q-score threshold for pass/fail 
min_qscore = 5

#Config file 
#guppy_config_file = 'dna_r9.4.1_450bps_hac.cfg'
guppy_config_file = 'dna_r10.4.1_e8.2_400bps_hac.cfg'
#guppy_config_file = 'dna_r10.4.1_e8.2_260bps_hac.cfg'
#guppy_config_file = 'dna_r10.4.1_e8.2_260bps_sup.cfg'
