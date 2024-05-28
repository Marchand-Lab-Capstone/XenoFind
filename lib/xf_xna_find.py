from Bio import SeqIO
import os
import pysam
import sys
import warnings
import setup_methods as setup
import raw_read_merger as rrm
import xf_basecall as bc
import feature_extraction as fe
from alive_progress import alive_bar
import xf_params as xfp

"""
Main script to call submethods to generate features for xna detection after 
consensus sequences are built. 

WIP
"""
def split_fasta(input_file_path, output_directory):

    fwd_fasta_path = os.path.join(output_directory, 'fwd_cons_seq.fasta')
    rev_fasta_path = os.path.join(output_directory, 'rev_cons_seq.fasta')
    fwd_seqs = []
    rev_seqs = []
    
    print('XenoFind [STATUS] - Splitting consensus fasta by direction')
    with open(input_file_path, 'r') as file:
        lines = file.readlines()
        current_header = ''
        current_seq = ''
        
        for line in lines:
            if line.startswith('>'):
                if current_header:
                    if 'fwd' in current_header:
                        fwd_seqs.append((current_header, current_seq))
                    elif 'rev' in current_header:
                        rev_seqs.append((current_header, current_seq))
                current_header = line.strip()
                current_seq = ''
            else:
                current_seq += line.strip()
        
        # Add the last sequence to the appropriate list
        if current_header:
            if 'fwd' in current_header:
                fwd_seqs.append((current_header, current_seq))
            elif 'rev' in current_header:
                rev_seqs.append((current_header, current_seq))
    
    with open(fwd_output_file_path, 'w') as fwd_file:
        for header, seq in fwd_seqs:
            fwd_file.write(f"{header}\n{seq}\n")
    
    with open(rev_output_file_path, 'w') as rev_file:
        for header, seq in rev_seqs:
            rev_file.write(f"{header}\n{seq}\n")

    return fwd_fasta_path, rev_fasta_path

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


    directories_list = setup.setup_directory_system_find(working_directory)
    
    #File path string for the merged pod5
    merged_pod5_path = os.path.join(directories_list[3], p5_fname + '.pod5')
    # add a parameter at top that takes in forced for all the functions 
    if xfp.regenerate_pod5 ==True or not (os.path.exists(merged_pod5_path)):
        # Using RRM, generate the pod5 from the data directory
        rrm.generate_merged_pod5(raw_data,
                                 merged_pod5_path)

    #-------Basecalling---------
    bc_bam_path = os.path.join(directories_list[4], basecall_fname+'.bam')
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
            os.system(bccmd)
        else:
            #Generate basecall command and run it
            bccmd = bc.dorado_bc_command(dorado_path,
                                     xfp.dorado_model_path,
                                     xfp.min_qscore,
                                     merged_pod5_path,
                                     bc_bam_path)
            os.system(bccmd)
    else: 
        print('Xenofind [STATUS] - basecalls found, skipping basecalling')
    
    #Alignment 
    align_cmd, aligned_bam_path = bc.alignment_command(bc_bam_path, reference_fasta, directories_list[4], direction)
    os.system(align_cmd)
    
    #filter out for primary only, forward direction reads

    filtered_bam_path = os.path.join(directories_list[4], direction+'_filtered.bam')
    #filter_cmd = f"samtools view -h -F 0x10 -bo {filtered_bam_path} {aligned_bam_path}"
    filter_cmd = f"samtools view -h -F 16 -F 2048 -F 2064 -bo {filtered_bam_path} {aligned_bam_path}"
    os.system(filter_cmd)
    
    return merged_pod5_path, filtered_bam_path
    
def feature_extraction(working_dir, merged_pod5, fwd_aligned_bam, rev_aligned_bam, fwd_xfasta, rev_xfasta):
    """
    feature_extractions runs the script to extract features from raw data & 
    sequence space and merge them. Outputs each reference sequence as a json 
    file containing the reads with their merged features. 
    """
    #Set up directory system
    directories_list = setup.setup_directory_system_find(working_dir)
    json_dir = directories_list[5] +'/'
    
    #fwd reference sequences 
    cmd = 'python lib/model_gen/data_concatenation.py '+merged_pod5+' '+fwd_aligned_bam+' '+fwd_xfasta+' '+json_dir
    os.system(cmd)
    
    #reverse_reference sequences 
    cmd = 'python lib/model_gen/data_concatenation.py '+merged_pod5+' '+rev_aligned_bam+' '+rev_xfasta+' '+json_dir
    os.system(cmd)
    
    return json_dir

def consensus_features(json_dir):
    """
    Takes in a the file path containing the json files with merged raw and 
    sequence space data. Calculates consensus features for each of json file. 
    Currently, data is stored in a list of pandas dataframes. 
    
    Inputs: 
    jsor_dir - directory containing json filse 
    
    return 
    cons_features_list - list of dataframes containing consensus level features 
    """
    warnings.filterwarnings("ignore")  # stops a warning from spamming your output
    sys.path.append('..//')  # path to directory holding feature_extraction
    json_file_names = os.listdir(json_dir)
    cons_features_list = []
    
    with alive_bar(len(json_file_names), title="Processing JSON files") as bar:
        for i in range(len(json_file_names)):
            json_file_path = os.path.join(json_dir, json_file_names[i])
            consensus_features = fe.feature_extraction(json_file_path, verbose=False)
            cons_features_list.append(consensus_features.T)
            print('Consensus sequence', i, 'features', consensus_features.T)
            bar()  # Update the progress bar
    return cons_features_list

def main():
    in_w_dir = sys.argv[1]
    in_r_dir = sys.argv[2]
    in_f_dir = sys.argv[3]
    
    #Create list of directories
    directories_list = setup.setup_directory_system_find(in_w_dir)
    ref_dir = directories_list[2]
    
   #Consensus fasta splitting 
    fwd_fasta, rev_fasta = split_fasta(in_f_dir)
        
    #fwd reads
    merged_pod5_path, fwd_filtered_bam_path = preprocessing(in_w_dir, in_r_dir, fwd_fasta, 'fwd')
    
    #rev reads
    merged_pod5_path, rev_filtered_bam_path = preprocessing(in_w_dir, in_r_dir, rev_fasta, 'rev')

    #Generate json files for forward and reverse reads 
    if xfp.regenerate_json or not os.listdir(directories_list[5]):
        json_dir = feature_extraction(in_w_dir, merged_pod5_path, fwd_filtered_bam_path, rev_filtered_bam_path, fwd_xfasta, rev_xfasta)
    else: 
        json_dir = directories_list[5]

    #Extract list consensus features 
    consensus_features_list = consensus_features(json_dir)

if __name__ == '__main__':
    main()
