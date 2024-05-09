from Bio import SeqIO
import os
import pysam
import sys
import setup_methods as setup
import raw_read_merger as rrm
import shannon_entropies as sn
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
def filter_primary_alignments(bam_path, out_name): #can edit this function to take in an input file name 
    """
    Filters a BAM file to retain only primary alignments and outputs the result
    to a new BAM file named 'primary_alignments.bam' in the same directory as
    the input file.
    
    Parameters: 
    bam_path: path to the BAM file as a string.
    
    Returns: 
    Path to the filtered BAM file containing only primary alignments.
    """
    
    # Generating the output BAM file path
    directory = os.path.dirname(bam_path)
    output_bam = os.path.join(directory, '{}.bam'.format(out_name))
    
    # Using pysam to filter for primary alignments
    with pysam.AlignmentFile(bam_path, "rb") as infile, \
         pysam.AlignmentFile(output_bam, "wb", header=infile.header) as outfile:

        for read in infile:
            if not read.is_secondary and not read.is_supplementary and not read.is_unmapped:
                outfile.write(read)
                
    print('XenoFind [STATUS] - Primary Only BAM file generated.')

    return output_bam

def preprocessing(working_directory, raw_data, reference_fasta, direction):
    """
    preprocessing will generate necessary files needed for feature extraction 
    or generation 
    
    WIP
    """

    # Defualt filenames:
    p5_fname = "merged"
    dorado_path = xfp.basecaller_path # assumes dorado is in user's home directory, make it a variable somewhere maybe
    basecall_fname = 'basecall' # fq file

    #----------Setup----------------------#
    # Set up the working directory 0
    #               xna_find_directory, 1
    #               merged_pod5, 2
    #               basecall_directory, 3
    #               forward_reads_directory, 4
    #               forward_reads_quality_score_directory, 5
    #               forward_reads_shannon_entropy_directory, 6
    #               forward_reads_signal_level_directory, 7
    #               reverse_reads_directory, 8
    #               reverse_reads_quality_score_directory, 9
    #               reverse_reads_shannon_entropy_directory, 10
    #               reverse_reads_signal_level_directory, 11
    #               reverse_unflipped_reads_directory, 12
    #               reverse_unflipped_reads_quality_score_directory, 13
    #               reverse_unflipped_reads_shannon_entropy_directory, 14
    #               reverse_unflipped_reads_signal_level_directory, 15
    #               total_reads_directory, 16
    #               total_reads_quality_score_directory, 17
    #               total_reads_shannon_entropy_directory, 18
    #               total_reads_signal_level_directory, 19

    directories_list = setup.setup_directory_system(working_directory)
    
    if direction == 'forward':
        index = 4
    if direction == 'reverse':
        index = 8
    if direction == 'reverse_unflipped':
        index = 12
    if direction == 'total':
        index = 16

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
                                     merged_pod5_path,
                                     directories_list[3],
                                     basecall_fname)
            st = os.system(bccmd)
        else:
            #Generate basecall command and run it
            bccmd = bc.dorado_bc_command(dorado_path,
                                     xfp.dorado_model_path,
                                     merged_pod5_path,
                                     directories_list[3],
                                     basecall_fname)
            st = os.system(bccmd)
    else: 
        print('Xenofind [STATUS] - basecalls found, skipping basecalling')
    
    #Alignment 
    align_cmd, aligned_bam_path = bc.alignment_command(bc_bam_path, reference_fasta, directories_list[index])
    os.system(align_cmd)
    
    #Primary Alignment Filtering
    filtered_bam_path = filter_primary_alignments(aligned_bam_path, 'filtered')
    
    return filtered_bam_path

def main():
    #in_w_dir = input("Please provide working directory path: ")
    #in_r_dir = input("Please provide read directory path: ")
    #in_f_dir = input("Please provide reference fasta directory path: ")
    
    #in_w_dir = sys.argv[1]
    #in_r_dir = sys.argv[2]
    #in_f_dir = sys.argv[3]
    in_w_dir = '/home/marchandlab/github/jay/capstone/XenoFind/xenofind_test/240419_BSn_90mer_SE_tests/' #Input desired working/ file output directory here
    in_r_dir = '/home/marchandlab/DataAnalysis/Kaplan/raw/fast5/10.4.1/240104_BSn_90mer_xr_train/50fast5' #Input either fast5 or pod5 containing directory here 
    con_fasta_fwd = '/home/marchandlab/github/jay/capstone/XenoFind/xenofind_test/240418_four_output_tests/consensus_files/xf_consensus_output/forward_reads/labeled_consensus.fasta' #Input consensus fasta or canonical ground truth fasta here
    con_fasta_rev = '/home/marchandlab/github/jay/capstone/XenoFind/xenofind_test/240418_four_output_tests/consensus_files/xf_consensus_output/reverse_reads/labeled_consensus.fasta' #Input consensus fasta or canonical ground truth fasta here
    con_fasta_rev_unflipped = '/home/marchandlab/github/jay/capstone/XenoFind/xenofind_test/240418_four_output_tests/consensus_files/xf_consensus_output/reverse_unflipped_reads/labeled_consensus.fasta' #Input consensus fasta or canonical ground truth fasta here
    con_fasta_all = '/home/marchandlab/github/jay/capstone/XenoFind/xenofind_test/240418_four_output_tests/consensus_files/xf_consensus_output/total_reads/labeled_consensus.fasta' #Input consensus fasta or canonical ground truth fasta here
    
    #Step 0: Rebasecalling/aligning data 
    read_types = ['forward', 'reverse', 'reverse_unflipped', 'total']
    consensus_list = [con_fasta_fwd, con_fasta_rev, con_fasta_rev_unflipped, con_fasta_all]
    
    #----------Setup----------------------#
    # Set up the working directory 0
    #               xna_find_directory, 1
    #               merged_pod5, 2
    #               basecall_directory, 3
    #               forward_reads_directory, 4
    #               forward_reads_quality_score_directory, 5
    #               forward_reads_shannon_entropy_directory, 6
    #               forward_reads_signal_level_directory, 7
    #               reverse_reads_directory, 8
    #               reverse_reads_quality_score_directory, 9
    #               reverse_reads_shannon_entropy_directory, 10
    #               reverse_reads_signal_level_directory, 11
    #               reverse_unflipped_reads_directory, 12
    #               reverse_unflipped_reads_quality_score_directory, 13
    #               reverse_unflipped_reads_shannon_entropy_directory, 14
    #               reverse_unflipped_reads_signal_level_directory, 15
    #               total_reads_directory, 16
    #               total_reads_quality_score_directory, 17
    #               total_reads_shannon_entropy_directory, 18
    #               total_reads_signal_level_directory, 19
    directories_list = setup.setup_directory_system(in_w_dir)
        

    print('XenoFind [STATUS] - Performing preprocessing')
    for i in range(len(consensus_list)):
        direction = read_types[i]
        if direction == 'forward':
            qs_outputs = directories_list[5]
            shannon_outputs = directories_list[6]
            signal_outputs = directories_list[7]
        if direction == 'reverse':
            qs_outputs = directories_list[9]
            shannon_outputs = directories_list[10]
            signal_outputs = directories_list[11]
        if direction == 'reverse_unflipped':
            qs_outputs = directories_list[13]
            shannon_outputs = directories_list[14]
            signal_outputs = directories_list[15]
        if direction == 'total':
            qs_outputs = directories_list[17]
            shannon_outputs = directories_list[18]
            signal_outputs = directories_list[19]
        filtered_bam_path = preprocessing(in_w_dir, in_r_dir, consensus_list[i], read_types[i])

        #Step 1: shannon entropies
        entropy_df, entropy_path = sn.wrapper(filtered_bam_path, consensus_list[i], direction, shannon_outputs, n=10, verbose=False) #this about how to subdivide this later 
        print(entropy_df)
        
    '''
    Need to decide if we simply call the shannon_entropies.py script or import functions
    '''


if __name__ == '__main__':
    main()
