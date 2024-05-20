from Bio import SeqIO
from Bio.Seq import Seq
import os
import pysam
import sys
import pandas as pd

import raw_read_merger as rrm
import setup_methods as setup
import consensus_methods_split as cs

script_dir = os.path.dirname(os.path.abspath(__file__))
parent_dir = os.path.dirname(script_dir)
if parent_dir not in sys.path:
    sys.path.append(parent_dir)

import xf_params as xfp

def consensus_generation(working_dir, reads, barcoded_fasta):
    """
    bc_align takes in working directory, reads, and a barcoded fasta file,
    and performs basecalling on them using Dorado then perform alignment using 
    minimap2 set with custom parameters to maximize the amount of reads that align 
    to the initial placeholder fasta
    
    Parameters
    working_dir: desired place where preprocessing files and final outputs will be, inputted as a string 
    reads: fast5 or pod5 directory, inputted as a string
    barcoded_fasta: File pathway to reference file for minimap2 alignment, inputted as a string
    """

    # Default filenames:
    p5_fname = "merged"
    dorado_path = xfp.basecaller_path # assumes dorado is in user's home directory, make it a variable somewhere maybe
    basecall_fname = 'basecall' # fq file
    minimap2_path = 'minimap2' #called from conda environment
    aligned_barcode_fname = 'barcode_aligned' # SAM file
    trimmed_fname = 'trimmed_' # Fasta 
    sorted_fname = 'sorted_' # Fasta
    vsearch_path = 'vsearch' # should be variable
    direction = ['fwd', 'rev']

    #indexes:
    fasta_index = 5
    vsearch_index = 6
    output_index = 8

    #----------Setup----------------------#
    # Set up the working directory 0
    #             consensus_files_directory, 1
    #             merged_pod5, 2
    #             basecall_directory, 3
    #             preprocessing_directory, 4
    #             preprocessing_fasta_directory, 5
    #             vsearch_processing, 6

    directories_list = setup.setup_directory_system(working_dir)

    #File path string for the merged pod5
    merged_pod5_path = directories_list[2] + p5_fname + '.pod5'
    # add a parameter at top that takes in forced for all the functions 
    if not (os.path.exists(merged_pod5_path)):
        # Using RRM, generate the pod5 from the data directory
        rrm.generate_merged_pod5(reads,
                                 directories_list[2],
                                 p5_fname)

    #-------Basecalling and Sorting ---------


    # Filepath string for the basecalled fq 
    basecalled_path = directories_list[3] + basecall_fname + '.fq'
    
    #Make this toggleable if data needs to be rebasecalled 
    if not (os.path.exists(basecalled_path)) or xfp.basecall_pod == True:
        print('Xenofind [STATUS] - Basecalling using Dorado')
        if xfp.auto_model == True:
            #Generate basecall command and run it 
            bccmd = cs.basecall_command(dorado_path,
                                     xfp.auto_model_type,
                                     merged_pod5_path,
                                     directories_list[3],
                                     basecall_fname)
            st = os.system(bccmd)
        else:
            #Generate basecall command and run it
            bccmd = cs.basecall_command(dorado_path,
                                     xfp.dorado_model_path,
                                     merged_pod5_path,
                                     directories_list[3],
                                     basecall_fname)
            st = os.system(bccmd)
    else: 
        print('Xenofind [STATUS] - basecalls found, skipping basecalling')
        
    # ensure the barcode is absolute.
    barcoded_fasta = str(os.path.abspath(barcoded_fasta))
    
    # use minimap2 to align the basecalled to the basecalled fq
    
    #add another conditional here to stop or force minimap2 alignment <--- still needs to be done #<---- this needs to be fixed, remove secondary alignments maybe
    map2ref_cmd = cs.map_to_reference(minimap2_path,
                                  barcoded_fasta,
                                  basecalled_path,
                                  directories_list[3],
                                  aligned_barcode_fname)
    # Run the minimap2 command
    st = os.system(map2ref_cmd)
    sam_file_path = '{}{}.sam'.format(directories_list[3],aligned_barcode_fname)
    
    #Generating a bam file from the sam file above 
    bam2sam_cmd, bam_file_path = cs.sam_to_bam(sam_file_path,
                                directories_list[3],
                                aligned_barcode_fname)
    st = os.system(bam2sam_cmd)
    
    #filter the bam file with primiary aligned reads
    filtered_bam_file_path = cs.filter_primary_alignments(bam_file_path, 'primary_read_filtered')

    #--------Trimming and Sorting Steps----------#

    # trim down the samfile to a trimmed fasta using default of 95% margin
    trimmed_lists_list = cs.read_trim(filtered_bam_file_path)
    
    consesus_fasta_list = []
    filtered_fasta_files = []
    
    for index, direction in enumerate(direction):
        similarity_id = xfp.starting_similarity

        # Retrieve the correct trimmed list for the current direction
        current_trimmed_list = trimmed_lists_list[index]

        # Generate trimmed file fasta for the current direction
        trimmed_fasta_path = cs.write_to_fasta(directories_list[fasta_index],
                                               trimmed_fname + f'{direction}',
                                               current_trimmed_list)

        # Sort the trimmed fasta, write it out.
        sorted_records_list = cs.sort_fasta(trimmed_fasta_path)
        sorted_fasta_path = cs.write_to_fasta(directories_list[fasta_index],
                                              sorted_fname,
                                              sorted_records_list)

        # Vsearch Steps
        vsearch_cmd, cluster_path = cs.vsearch_command(vsearch_path, 
                                                       sorted_fasta_path, 
                                                       directories_list[vsearch_index], 
                                                       f'first_vsearch_{direction}', 
                                                       similarity_id)
        os.system(vsearch_cmd)

        # Initialize dataframe
        df = pd.DataFrame()
        df = cs.cluster_size_df(cluster_path, df)

        # Proceed with further operations as before...
        for j in range(xfp.vsearch_iterations):
        
            similarity_id += xfp.similarity_increment
            # Further commands involving realignment and reanalysis as previously described
            cluster_align_cmd, sam_path = cs.mm2_cluster_aligner(minimap2_path, cluster_path, trimmed_fasta_path, directories_list[vsearch_index], f'aligned_{direction}_cluster_{j+1}')

            # Run the minimap2 command
            st = os.system(cluster_align_cmd)
            #Convert the outputted sam file into a bam file 
            sam2bam_cmd, bam_path = cs.sam_to_bam(sam_path, directories_list[vsearch_index], f'aligned_{direction}_cluster_{j+1}')
            st = os.system(sam2bam_cmd)

            #Filter the bam file for only primary alignments 
            filtered_bam = cs.filter_primary_alignments(bam_path, f'primary_{direction}_aligned_{j+1}')

            reference_counts = cs.weight_generation(filtered_bam)
            weighted_fasta = cs.weighted_fasta_gen(cluster_path, reference_counts, f'weighted_{direction}_{j+1}')

            #Now reperform VSEARCH steps using weighted fasta 

            #need to sort weighted fasta file 
            weighted_sorted_list = cs.sort_fasta(weighted_fasta)
            weighted_sorted_fasta_path = cs.write_to_fasta(
                                                        directories_list[vsearch_index],
                                                        f'weighted_sorted_{direction}_{j+1}',
                                                         weighted_sorted_list)

            
            vsearch_cmd, cluster_path = cs.vsearch_command(vsearch_path,
                                                           weighted_sorted_fasta_path,
                                                           directories_list[vsearch_index],
                                                          f'weighted_consensus_{direction}_{j+1}',
                                                          similarity_id) #get rid of hard coded 0.9, 
            st = os.system(vsearch_cmd)

            print(f'The {j+1} Rounds of consensus generation is complete for direction {direction}')

        # Save dataframe to CSV
        df.to_csv(os.path.join(os.path.dirname(cluster_path), f'cluster_size_{direction}.csv'))

        # Handle consensus labeling and final processing
        lab_cons_path = directories_list[vsearch_index] + f'labeled_consensus_{direction}.fasta'
        
        lab_cons_path = cs.rename_consensus_headers(cluster_path, lab_cons_path)
        
        filtered_fasta_files.append(cs.filter_cluster_size(lab_cons_path, xfp.cluster_size_threshold, direction))
    
    final_fasta_path = cs.merge_fasta_files(filtered_fasta_files, os.path.join(directories_list[1], 'final_merged_consensus.fasta'))

    return final_fasta_path

def main():
    #in_w_dir = input("Please provide working directory path: ")
    #in_r_dir = input("Please provide read directory path: ")
    #in_f_dir = input("Please provide reference fasta directory path: ")
    
    in_w_dir = sys.argv[1]
    in_r_dir = sys.argv[2]
    in_f_dir = sys.argv[3]
    
    total_consensus_fasta = consensus_generation(in_w_dir, in_r_dir, in_f_dir)
     
    print('Xenofind [STATUS] - Generated consensus fasta containing forward and reverse strands at', total_consensus_fasta) 
    


if __name__ == '__main__':
    main()
