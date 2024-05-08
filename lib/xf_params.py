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
basecall_pod = True

#Perform Quality Score Analysis 
analyze_fastq = True
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

##Model Training Parameters
############################################################
#Dorado Base caller configuration

#Path to guppy basecaller
basecaller_path = '~/dorado-0.6.0-linux-x64/bin/dorado' #consider adding Guppy compatibility 

#Dorado Model Parameters
min_qscore = 5

#Let Dorado automatically choose a model 
auto_model = True 

#Model type dorado will search for the best model based on input (options: fast, hac, sup)
auto_model_type = 'hac' 

#If manual model selection is desired
dorado_model_path = '' 



