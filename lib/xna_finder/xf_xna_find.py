from Bio import SeqIO
import os
import pysam
import sys
import setup_methods as setup
import raw_read_merger as rrm
import xf_basecall as bc

script_dir = os.path.dirname(os.path.abspath(__file__))
parent_dir = os.path.dirname(script_dir)
if parent_dir not in sys.path:
    sys.path.append(parent_dir)

import xf_params as xfp

"""
Main script to call submethods to generate features for xna detection after 
consensus sequences are built. 

WIP
"""


def preprocessing(working_directory, raw_data, reference_fasta):
    """
    preprocessing will generate necessary files needed for feature extraction 
    or generation 
    
    WIP
    """

    # Defualt filenames:
    p5_fname = "merged"
    dorado_path = xfp.basecaller_path # assumes dorado is in user's home directory, make it a variable somewhere maybe
    basecall_fname = 'basecall' # fq file

    #1. working directory
    #2. xna find 
    #3. pod5 directory 
    #4. basecall directory 
    directories_list = setup_methods(working_directory)

    #File path string for the merged pod5
    merged_pod5_path = directories_list[2] + p5_fname + '.pod5'
    # add a parameter at top that takes in forced for all the functions 
    if not (os.path.exists(merged_pod5_path)):
        # Using RRM, generate the pod5 from the data directory
        rrm.generate_merged_pod5(reads,
                                 directories_list[2],
                                 p5_fname)

    #-------Basecalling---------
    
    #Make this toggleable if data needs to be rebasecalled 
    if not (os.path.exists(basecalled_path)) or xfp.basecall_pod == True:
        print('Xenofind [STATUS] - Basecalling using Dorado')
        if xfp.auto_model == True:
            #Generate basecall command and run it 
            bccmd, bc_bam_path = bc.dorado_command(dorado_path,
                                     xfp.auto_model_type,
                                     merged_pod5_path,
                                     directories_list[4],
                                     basecall_fname)
            st = os.system(bccmd)
        else:
            #Generate basecall command and run it
            bccmd, bc_bam_path = bc.dorado_command(dorado_path,
                                     xfp.dorado_model_path,
                                     merged_pod5_path,
                                     directories_list[4],
                                     basecall_fname)
            st = os.system(bccmd)
    else: 
        print('Xenofind [STATUS] - basecalls found, skipping basecalling')

    bc_bam_path = directories_list[4] + basecall_fname
    #Primary Alignment Filtering
    filtered_bam_path = filter_primary_alignments(bc_bam_path)
    
    return filtered_bam_path

def main():
    #in_w_dir = input("Please provide working directory path: ")
    #in_r_dir = input("Please provide read directory path: ")
    #in_f_dir = input("Please provide reference fasta directory path: ")
    
    in_w_dir = sys.argv[1]
    in_r_dir = sys.argv[2]
    in_f_dir = sys.argv[3]
    
    print('XenoFind [STATUS] - Performing preprocessing')
    filtered_bam_path = preprocessing(in_w_dir, in_r_dir, in_f_dir)
    
    """
    set 'consensus_generation' above to a variable to return consensus pathway
    """


if __name__ == '__main__':
    main()