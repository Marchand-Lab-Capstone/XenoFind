from Bio import SeqIO
import os
import pysam
import sys
import warnings
import setup_methods as setup
import raw_read_merger as rrm
import xf_basecall as bc
import feature_extraction as fe

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
    merged_pod5_path = os.path.join(directories_list[2], p5_fname + '.pod5')
    # add a parameter at top that takes in forced for all the functions 
    if not (os.path.exists(merged_pod5_path)):
        # Using RRM, generate the pod5 from the data directory
        rrm.generate_merged_pod5(raw_data,
                                 merged_pod5_path)

    #-------Basecalling---------
    bc_bam_path = os.path.join(directories_list[3], basecall_fname+'.bam')
    #Make this toggleable if data needs to be rebasecalled 
    if not (os.path.exists(bc_bam_path)) or xfp.basecall_pod == True:
        print('Xenofind [STATUS] - Basecalling using Dorado')
        if xfp.auto_model == True:
            #Generate basecall command and run it 
            bccmd = bc.dorado_bc_command(dorado_path,
                                     xfp.auto_model_type,
                                     xfp.min_qscore,
                                     merged_pod5_path,
                                     bc_bam_path)
            st = os.system(bccmd)
        else:
            #Generate basecall command and run it
            bccmd = bc.dorado_bc_command(dorado_path,
                                     xfp.dorado_model_path,
                                     xfp.min_qscore,
                                     merged_pod5_path,
                                     bc_bam_path)
            st = os.system(bccmd)
    else: 
        print('Xenofind [STATUS] - basecalls found, skipping basecalling')
    
    #Alignment 
    align_cmd, aligned_bam_path = bc.alignment_command(bc_bam_path, reference_fasta, directories_list[3])
    os.system(align_cmd)
    
    return merged_pod5_path, aligned_bam_path
    
def consensus_features(working_dir, merged_pod5, aligned_bam, ref_fasta):
    """
    """
    directories_list = setup.setup_directory_system(working_dir)
    
    if xfp.regenerate_json or not os.listdir(directories_list[5]):
        cmd = 'python lib/model_gen/data_concatenation.py '+merged_pod5+' '+aligned_bam+' '+ref_fasta+' '+ directories_list[5]+'/'
        os.system(cmd)
    
    warnings.filterwarnings("ignore") # stops a warning from spamming your output
    sys.path.append('..//') # path to directory holding feature_extraction
    json_dir = directories_list[5]+'/'
    json_file_names = os.listdir(json_dir)
    cons_features_list = []
    for i in range(len(json_file_names)): # can be adjusted to the number of files you want
        json_file_path = os.path.join(json_dir, json_file_names[i])
        consensus_features = fe.feature_extraction(json_file_path, verbose=False)
        cons_features_list.append(consensus_features.T)
        print('Consensus sequence', i, 'features', consensus_features.T)
    
    return cons_features_list
def main():
    #in_w_dir = input("Please provide working directory path: ")
    #in_r_dir = input("Please provide read directory path: ")
    #in_f_dir = input("Please provide reference fasta directory path: ")
    
    '''
    Need to chagne this to either: A, take in both a forward and reverse reference fasta file based on work from last week. B take in a fasta file containing both forward and reverse strands and auto split it.
    '''
    in_w_dir = sys.argv[1]
    in_r_dir = sys.argv[2]
    in_f_dir = sys.argv[3]
    '''
    if we had both forward and reverse reads decoupled, the preprocessing function would get called twice 
    '''
    #xfasta_path = 
    merged_pod5_path, aligned_bam_path = preprocessing(in_w_dir, in_r_dir, in_f_dir)
    '''
    Need to decide if we simply call the shannon_entropies.py script or import functions
    '''
    consensus_features_list = consensus_features(in_w_dir, merged_pod5_path, aligned_bam_path, in_f_dir)

if __name__ == '__main__':
    main()
