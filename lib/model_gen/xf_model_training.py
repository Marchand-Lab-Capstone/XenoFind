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


    directories_list = setup.setup_directory_system(working_directory)
    
    #File path string for the merged pod5
    merged_pod5_path = directories_list[2] + p5_fname + '.pod5'
    # add a parameter at top that takes in forced for all the functions 
    if not (os.path.exists(merged_pod5_path)):
        # Using RRM, generate the pod5 from the data directory
        rrm.generate_merged_pod5(raw_data,
                                 directories_list[2],
                                 p5_fname)

    #-------Basecalling---------
    bc_bam_path = os.path.join(directories_list[3],basecall_fname)+'.bam'
    #Make this toggleable if data needs to be rebasecalled 
    if not (os.path.exists(bc_bam_path)) or xfp.basecall_pod == True:
        print('Xenofind [STATUS] - Basecalling using Dorado')
        if xfp.auto_model == True:
            #Generate basecall command and run it 
            bccmd = bc.dorado_bc_command(dorado_path,
                                     xfp.auto_model_type,
                                     xfp.min_qscore,
                                     merged_pod5_path,
                                     directories_list[3],
                                     basecall_fname)
            st = os.system(bccmd)
        else:
            #Generate basecall command and run it
            bccmd = bc.dorado_bc_command(dorado_path,
                                     xfp.dorado_model_path,
                                     xfp.min_qscore,
                                     merged_pod5_path,
                                     directories_list[3],
                                     basecall_fname)
            st = os.system(bccmd)
    else: 
        print('Xenofind [STATUS] - basecalls found, skipping basecalling')
    
    #Alignment 
    align_cmd, aligned_bam_path = bc.alignment_command(bc_bam_path, reference_fasta, directories_list[3])
    os.system(align_cmd)
    
    return aligned_bam_path

def main():
    #in_w_dir = input("Please provide working directory path: ")
    #in_r_dir = input("Please provide read directory path: ")
    #in_f_dir = input("Please provide reference fasta directory path: ")
    
    #in_w_dir = sys.argv[1]
    #in_r_dir = sys.argv[2]
    #in_f_dir = sys.argv[3]
    in_w_dir = '/home/marchandlab/github/jay/capstone/XenoFind/xenofind_test/240508_PZ_xm_lib_model_training_development/' #Input desired working/ file output directory here
    in_r_dir = '/home/marchandlab/DataAnalysis/Sumabat/230725_PZ_lib_v4_r10/20230725_1220_MN37138_APH167_a204cb54/fast5' #Input either fast5 or pod5 containing directory here 
    in_f_fwd_dir = '/home/marchandlab/github/jay/capstone/reference/xPZ_xm_libv4_full.fa'
    in_f_rev_dir = '/home/marchandlab/github/jay/capstone/reference/xPZ_xm_libv4_full.fa;
    #Step 0: Rebasecalling/aligning data 
    fwd_aligned_bam_path = preprocessing(in_w_dir, in_r_dir, in_f_rev_dir)
    rev_aligned_bam_path = preprocessing(in_w_dir, in_r_dir, in_f_rev_dir)
    '''
    Need to decide if we simply call the shannon_entropies.py script or import functions
    '''


if __name__ == '__main__':
    main()
