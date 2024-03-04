import os
import glob
import sys
import pysam
import numpy as np
from pathlib import Path
from xf_params import * 
from xf_tools import *

#Initialize
#working_dir = os.path.expanduser(sys.argv[1]) #will be final impolemented version, will be testing in large scale script later 
#xna_raw_dir = os.path.expanduser(sys.argv[2])


#Initialize (manual input of datasets) 

#using dataset nick generated to test 
working_dir = '/home/marchandlab/github/jay/capstone/xenofind/xenofind_test/240220_lq_tests'
raw_dir = '/home/marchandlab/DataAnalysis/Kaplan/raw/fast5/10.4.1/240104_BSn_90mer_xr_train/50fast5'
ref_fasta = '/home/marchandlab/DataAnalysis/Kaplan/raw/fast5/10.4.1/240104_BSn_90mer_xr_train/reference/BSn_90mer_xr_train_can.fa' #this is ground truth for now, will be substituted by consensus sequence formation pipeline once that is done (inside lab project)


# Generate Required Directories
working_dir = check_make_dir(working_dir)
processing_dir = check_make_dir(os.path.join(working_dir, 'processing/'))
pod_dir = check_make_dir(os.path.join(processing_dir, 'pod5/'))
bc_dir = check_make_dir(os.path.join(processing_dir, 'basecall/'))
bam_dir = check_make_dir(os.path.join(processing_dir, 'bam/'))
fastq_dir = check_make_dir(os.path.join(processing_dir, 'fastq/'))

'''
NOTE: 
WE CAN ADD A PART WHICH WILL AUTO INSTALL DORADO INTO THIS WORKING DIRECTORY AFTER THE FIRST TIME THE PROGRAM IS CALLED. MIGHT BE WORTH DOING
'''
#Step 1: Generate or merge pod5 files if needed
file_type = os.listdir(raw_dir)
# Check if the directory is not empty
if file_type:
    # Get the first file in the directory
    first_file = file_type[0]
    # Check if the first file is a .pod5 file
    if first_file.endswith(".pod5"):
        if os.path.isfile(os.path.join(pod_dir,os.path.basename(raw_dir))+'.pod5')==False:
            pod5_merge(get_pod5_subdir(raw_dir), os.path.join(pod_dir,os.path.basename(raw_dir))+'.pod5')
        else:
            print('XenoFind [STATUS] - POD5 files merged. Skipping merge')
    else:
        if os.path.isfile(os.path.join(pod_dir,os.path.basename(raw_dir))+'.pod5')==False:
            cod5_to_fast5(get_fast5_subdir(raw_dir), os.path.join(pod_dir,os.path.basename(raw_dir))+'.pod5')
        else: 
            print('XenoFind [STATUS] - Converted POD5 file for modified base found. Skipping POD5 coversion')
else:
    print('XenoFind [ERROR] - Modified Fast5/POD5 directory empty, please check input directory')
    sys.exit()
    
if basecall_pod == True:
    #Perform an initial basecalling with a reference file to generate a bam file with all necessary read information 
    print('XenoFind [STATUS] - Performing basecalling using Dorado')
    cmd = os.path.expanduser(basecaller_path)+ ' basecaller hac  --no-trim  ' + pod_dir + ' > '+os.path.join(bc_dir, 'bc.bam') + ' --reference ' + ref_fasta#can probably do this in a bam file as well 
    os.system(cmd)

    #Data filtering here maybe, leaving all data reads for now 
if analyze_fastq == True: 
    
    def extract_read_info(bam_file_path):
        """ 
        This function takes in the bam file generated and extracts the readID, basecalled sequence, start of reference sequence, and the quality score from the bamfile 
        """
        read_info = []
        # Open the BAM file
        with pysam.AlignmentFile(bam_file_path, "rb") as bamfile:
            for read in bamfile:
                # Check if the read is mapped
                if not read.is_unmapped:
                    qual = read.query_qualities 
                    qual = np.array(qual, dtype=int)
                    avg_qual = np.mean(qual)
                    features = [
                    read.query_name,  # Query name of the read
                    read.query_sequence,  # Sequence basecalled
                    read.reference_start,  # Position of the read relative to the reference
                    qual,  # Quality scores of the read (numerical)
                    avg_qual
                    ]
                    read_info.append(features)
                else:
                    print(f"Read {read.query_name} is unmapped and does not have a reference position.")
            return read_info
    read_info = extract_read_info(os.path.join(bc_dir, 'bc.bam'))
    
# Test if the every base pair gets its own quality score
if len(read_info[0][1]) == len(read_info[0][3]):
    print('Yes, the length of the sequences matches the length of quality score string')

else:
    print('No, the length of the sequences and quality score doesnt match!')
    
    #print(read_info) #little more than 50% alignment to ground truth for single sequence context#little more than 50% alignment to ground truth for single sequence context

    #Predict XNA position using quality string analysis 
    #Need to make it so reads are grouped by which sequence its aligned to (our current dataset is single sequence context)
    def xna_guess(read_info):
        """
        This function will perform a statistical test on a per read basis 
        """
        q_score = read_info
