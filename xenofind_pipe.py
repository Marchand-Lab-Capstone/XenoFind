######################################################

'''
pipe.py 

By: J. Sumabat, N. Lai, Y. Huang, S. Peck 

'''
######################################################

import os 
import numpy as np 
import glo 
from pathlib import Path 

######################################################

# Paths 
working_dir = '' 
raw_dir = ''
ref_fasta = '' #This should be a placeholder reference fasta not necessarily ground truth, consider adding a new variable for consensus fasta.
######################################################


######################################################

consensus = True 
lq_detect = False #placeholder variables 
method_2 = False
method_3 = False

######################################################

# Generate consensus fasta 
if consensus == True: 
		cmd = 'python xenofind.py consensus -w ' + working_dir+' -f '+ raw_dir+' -r '+ref_fasta
		os.system(cmd)
# Low Quality Score Detection 
if lq_detect == True: 
	cmd = 'python xenofind.py low_qual -w '+working_dir+' -f '+raw_dir+' -r '+ref_fasta
	os.system(cmd)

# Calls for method 3 
if method_3 == True: 
	cmd = ''
	os.system(cmd)
