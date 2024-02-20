import os
import glob
import sys
import pysam
from pathlib import Path
from xf_params import * 
from xf_tools import *

#Initialize
#working_dir = os.path.expanduser(sys.argv[1]) #will be final impolemented version, will be testing in large scale script later 
#xna_raw_dir = os.path.expanduser(sys.argv[2])


#Initialize (manual input of datasets) 
working_dir = '/home/marchandlab/github/jay/capstone/xenofind/xenofind_test/240216_lq_tests'
raw_dir = '/home/marchandlab/DataAnalysis/Sumabat/231117_PZn_libv4_minion/20231117_1644_MN37138_FAX70078_76df2751/subset_fast5'
ref_fasta = '' #this is ground truth for now, will be substituted by consensus sequence formation pipeline once that is done (inside lab project)

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
''''
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
    print('XenoFind [STATUS] - Performing basecalling using Dorado')
    cmd = os.path.expanduser(basecaller_path)+ ' basecaller hac  --no-trim  ' + pod_dir + ' > '+os.path.join(bc_dir, 'bc.bam') #can probably do this in a bam file as well 
    os.system(cmd)

if analyze_fastq == True: 
    
    def extract_read_info(bam_file_path):
        """ 
        ChatGPT first pass for relevant bam info extraction
        """
        # Open the BAM file
        with pysam.AlignmentFile(bam_file_path, "rb") as bamfile:
            for read in bamfile:
                # Check if the read is mapped
                if not read.is_unmapped:
                    # Position of the read relative to the reference
                    position = read.reference_start
                    
                    # Query name of the read
                    query_name = read.query_name
                    
                    # Quality string of the read
                    quality_string = read.qual
                    
                    # Sequence basecalled
                    sequence_basecalled = read.query_sequence
                    
                    # Quality score string (ASCII representation)
                    quality_score_string = read.query_qualities
                    
                    # Converting quality score string from array of integers to ASCII
                    quality_score_string_ascii = ''.join([chr(q + 33) for q in quality_score_string])
                    
                    # You can process these variables further as per your requirements
                    print(f"Query Name: {query_name}, Position: {position}, Sequence: {sequence_basecalled}, Quality String: {quality_string}, Quality Score String: {quality_score_string_ascii}")
                else:
                    print(f"Read {read.query_name} is unmapped and does not have a reference position.")
    extract_read_info(os.path.join(bc_dir, 'bc.bam'))
