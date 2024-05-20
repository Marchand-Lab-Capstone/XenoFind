########################################################################
########################################################################
"""
xr_params.py 

Title: Unpublished work

By: H. Kawabe, N. Kaplan, J. Sumabat, J. A. Marchand

Updated: 11/28/23
"""
########################################################################
#imports 
import numpy as np


############################################################
##Analysis instructions 

#Remerge fast5 or pod5 files into single pod5 file 
regenerate_pod5 = True

#Re-basecall pod5 file. Required if new reference files are being used. 
basecall_pod = False

#Perform Quality Score Analysis 
analyze_fastq = True
######################XFASTA GENERATION######################
#Fasta2x - write sequences to xfasta even if they do no contain XNAs. Default = False 
write_no_xna_seq = False

#Fasta2x - Write gaps in place of XNAs in fasta reference file for null testing
write_gaps = False

############################################################
##Consensus Parameters

#Starting similarity ID, initial similarity ID prior to doing VSEARCH iterations
starting_similarity = 0.7

#Increment, how much to increase the similarity ID by during VSEARCH iterations
similarity_increment = 0.05

#Number of Weighted VSEARCH Rounds to Perform
vsearch_iterations = 3

##Find Parameters 
shannon_entropy_graphs = False 

#cluster size filter threshold
cluster_size_threshold = 100

############################################################
######################XFASTA GENERATION#####################
#Fasta2x - write sequences to xfasta even if they do no contain XNAs. Default = False 
write_no_xna_seq = False

#Fasta2x - Write gaps in place of XNAs in fasta reference file for null testing
write_gaps = False

##Model Training Parameters
regenerate_json = False

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

############################################################
#Dorado Base caller configuration

#Path to dorado basecaller
#basecaller_path = 'dorado'
basecaller_path = '/home/marchandlab/dorado-0.6.0-linux-x64/bin/dorado'

#Dorado Model Parameters
min_qscore = 5

#Let Dorado automatically choose a model 
auto_model = True 

#Model type dorado will search for the best model based on input (options: fast, hac, sup)
auto_model_type = 'hac' 

#If manual model selection is desired
dorado_model_path = '' 



