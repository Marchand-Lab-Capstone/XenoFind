########################################################################
########################################################################
"""
xr_params.py 

Title: Unpublished work

By: H. Kawabe, N. Kaplan, J. Sumabat, J. A. Marchand

Updated: 11/28/23
"""
########################################################################
########################################################################

import numpy as np



#Standard basepairs written in 'purine pyrimidine' order
standard_base_pairs = ['AT','GC', 'NN']

#Convert this to set
standard_bases = np.concatenate(list(list(i) for i in standard_base_pairs))

#Alternative basepairs written in 'purine pyrimidine' order
xna_base_pairs = ['BS','PZ','JV','XK']

#Specify canonical base substitution desired for xFASTA generation here
confounding_pairs =  ['BA','SA','PG','ZC','JC','VG','XA','KG'] 

#If XNAs are given different standard base substitutions, set them up as seperate (e.g, ['P','Z'])
xna_segmentation_model_sets = ['B','S','PZ','JV','X', 'K', 'QW','ER']

#Possible XNA bases
xna_bases = np.concatenate(list(list(i) for i in xna_base_pairs))


######################XFASTA GENERATION######################

#Fasta2x - write sequences to xfasta even if they do no contain XNAs. Default = False 
write_no_xna_seq = False

#Fasta2x - Write gaps in place of XNAs in fasta reference file for null testing
write_gaps = False

############################################################
##Analysis instructions 

#Re-basecall pod5 file. Required if new reference files are being used. 
basecall_pod = False

#Perform Quality Score Analysis 
analyze_fastq = True
############################################################


############################################################
##Model Training and Basecalling Parameters

#kmer table 
#kmer_table_path = 'models/remora/4mer_9.4.1.csv'
kmer_table_path = 'models/remora/9mer_10.4.1.csv'

#ml model (ConvLSTM_w_ref.py or Conv_w_ref.py')
ml_model_path = 'models/ConvLSTM_w_ref.py'

#Modified base in Fasta sequence you wish to train model or use model to basecall
mod_base = 'Z'

#Most similar substituted canonical base you will be comparing against 
can_base = 'C'

#Extent of Kmer content (-,+) to store for model training
kmer_context ='4 4' 

#Extent of chunk context (centered around modified base) 
chunk_context = '50 50' 

#Proportion of reads to use for validation 
val_proportion = '0.2'

#Number of chunks for training (in thousands: e.g.: '200' = 200,000 chunks) 
chunk_num = '500000'




############################################################
#Guppy Base caller configuration

#Path to guppy basecaller
basecaller_path = '~/dorado-0.5.3-linux-x64/bin/dorado' #update to latest version tesing right now 
#basecaller_path = '~/dorado-0.5.3-osx-arm64/bin/dorado'

#GPU enabled 
device_type = 'cuda:all' 

#Guppy q-score threshold for pass/fail 
min_qscore = 5

#Config file 
#guppy_config_file = 'dna_r9.4.1_450bps_hac.cfg'
guppy_config_file = 'dna_r10.4.1_e8.2_400bps_hac.cfg'
#guppy_config_file = 'dna_r10.4.1_e8.2_260bps_hac.cfg'
#guppy_config_file = 'dna_r10.4.1_e8.2_260bps_sup.cfg'
