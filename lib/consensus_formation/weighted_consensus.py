from Bio import SeqIO
import os
import pysam
import sys
import raw_read_merger as rrm
import setup_methods as setup
import consensus_methods as cs

script_dir = os.path.dirname(os.path.abspath(__file__))
parent_dir = os.path.dirname(script_dir)  # Move up to 'lib' from 'lib2'
if parent_dir not in sys.path:
    sys.path.append(parent_dir)

from xf_params import *

def bc_align(working_dir, reads, barcoded_fasta):
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
    
    # Defualt filenames:
    p5_fname = "merged"
    #dorado_path = "dorado" # Should be variable, this assumes dorado is in user's PATH
    dorado_path = ' ~/dorado-0.5.3-linux-x64/bin/dorado' # assumes dorado is in user's home directory, make it a variable somewhere maybe
    basecall_fname = 'basecall' # fq file
    minimap2_path = 'minimap2' #called from conda environment
    aligned_barcode_fname = 'barcode_aligned' # SAM file
    trimmed_fname = 'trimmed' # Fasta 
    sorted_fname = 'sorted' # Fasta
    vsearch_path = 'vsearch' # should be variable
    vs_cons_fname = 'cons' # Fasta
    filtered_fname = 'represented_seq' # Fasta
    medaka_path = 'medaka_consensus' # should be variable
    polished_fname = 'polished_consensus' # Fasta

    # Get the barcode indexes using extract_n_indexes
    n_positions = cs.extract_n_indexes(barcoded_fasta)


    #----------Setup----------------------#
    # Set up the working directory 0
    #             basecall_directory, 1
    #             forward_reads_directory, 2
    #             forward_reads_fasta_directory, 3
    #             reverse_reads_directory, 4
    #             reverse_reads_fasta_directory, 5
    #             total_reads_directory, 6
    #             total_reads_fasta_directory, 7
    #             merged_pod5, 8
    #             vsearch_forward_processing, 9
    #             vsearch_reverse_processing, 10
    #             vsearch_total_processing, 11
    #             xf_consensus_output, 12
    directories_list = setup.setup_directory_system(working_dir)

    #File path string for the merged pod5
    merged_pod5_path = directories_list[8] + p5_fname + '.pod5'
    # add a parameter at top that takes in forced for all the functions 
    if not (os.path.exists(merged_pod5_path)):
        # Using RRM, generate the pod5 from the data directory
        rrm.generate_merged_pod5(reads,
                                 directories_list[8],
                                 p5_fname)

    #-------Basecalling and Sorting ---------


    # Filepath string for the basecalled fq 
    basecalled_path = directories_list[1] + basecall_fname + '.fq'
    
    #Make this toggleable if data needs to be rebasecalled 
    if not (os.path.exists(basecalled_path)) or basecall_pod == True:
    # Generate the dorado basecall command 
        print('Xenofind [STATUS] - Basecalling using Dorado')
        bccmd = cs.basecall_command(dorado_path,
                                 merged_pod5_path,
                                 directories_list[1],
                                 basecall_fname)
        # Run the basecall command
        st = os.system(bccmd)
    else: 
        print('Xenofind [STATUS] - basecalls found, skipping basecalling')
        
    # ensure the barcode is absolute.
    barcoded_fasta = str(os.path.abspath(barcoded_fasta))
    
    # use minimap2 to align the basecalled to the basecalled fq
    
    #add another conditional here to stop or force minimap2 alignment <--- still needs to be done 
    map2ref_cmd = cs.map_to_reference(minimap2_path,
                                  barcoded_fasta,
                                  basecalled_path,
                                  directories_list[1],
                                  aligned_barcode_fname)
    # Run the minimap2 command
    st = os.system(map2ref_cmd)
    sam_file_path = '{}{}.sam'.format(directories_list[1],aligned_barcode_fname)
    
    #Generating a bam file from the sam  file above 
    bam2sam_cmd, bam_file_path = cs.sam_to_bam(sam_file_path,
                                directories_list[1],
                                aligned_barcode_fname)
    st = os.system(bam2sam_cmd)

    return bam_file_path

def decouple(bam_file_path):
    """
    Separates reads from an aligned BAM file into forward and reverse strand reads
    using samtools. 
    
    Parameters: 
    bam_file_path: path to the BAM file, inputted as a string.
    
    Returns: 
    Paths to forward and reverse strand BAM files as strings. Also generates these files in the same directory.
    """
    # Correct variable names for clarity
    forward_out_path = os.path.join(os.path.dirname(bam_file_path), 'forward.bam')
    reverse_out_path = os.path.join(os.path.dirname(bam_file_path), 'reverse.bam')
    
    # Open the input BAM file for reading
    with pysam.AlignmentFile(bam_file_path, "rb") as infile:  # Note "rb" for reading BAM

        # Open two files for writing: one for forward strand reads, another for reverse strand reads
        # Note "wb" for writing BAM files
        with pysam.AlignmentFile(forward_out_path, "wb", header=infile.header) as outfile_forward, \
             pysam.AlignmentFile(reverse_out_path, "wb", header=infile.header) as outfile_reverse:

            # Iterate through reads in the input file
            for read in infile:
                # Check if the read is mapped to the reverse strand
                if read.is_reverse:
                    # Write the read to the reverse strand reads file
                    outfile_reverse.write(read)
                else:
                    # Write the read to the forward strand reads file (implicitly forward if not reverse)
                    outfile_forward.write(read)
    
    return forward_out_path, reverse_out_path
    
def consensus_generation(working_dir, bam_file_path, direction):

    #Default File names 
    trimmed_fname = 'trimmed' # Fasta 
    sorted_fname = 'sorted' # Fasta
    vsearch_path = 'vsearch' # should be variable
    vs_cons_fname = 'cons' # Fasta
    filtered_fname = 'represented_seq' # Fasta
    medaka_path = 'medaka_consensus' # should be variable
    polished_fname = 'polished_consensus' # Fasta
    minimap2_path = 'minimap2'
    similarity_id = '0.80'
    
    #----------Setup----------------------#
    # Set up the working directory 0
    #             basecall_directory, 1
    #             forward_reads_directory, 2
    #             forward_reads_fasta_directory, 3
    #             reverse_reads_directory, 4
    #             reverse_reads_fasta_directory, 5
    #             total_reads_directory, 6
    #             total_reads_fasta_directory, 7
    #             merged_pod5, 8
    #             vsearch_forward_processing, 9
    #             vsearch_reverse_processing, 10
    #             vsearch_total_processing, 11
    #             xf_consensus_output, 12
    directories_list = setup.setup_directory_system(working_dir)
    
    if direction == 'forward': #forward directory isnt calling anything right now, need to fix 
        fasta_index = 3 
        vsearch_index = 9
        #insert a new required index here as necessary
    elif direction == 'reverse':
        fasta_index = 5
        vsearch_index = 10
        #insert a new required index here as necessary
    else: 
        fasta_index = 7
        vsearch_index = 11
        #insert a new required index here as necessary

    #--------Trimming and Sorting Steps----------#

    # trim down the samfile to a trimmed fasta using default of 95% margin
    read_trims_list = cs.read_trim(bam_file_path)
    trimmed_fasta_path = cs.write_to_fasta(directories_list[fasta_index],
                                        trimmed_fname,
                                        read_trims_list)

    # Sort the trimmed fasta, write it out.
    sorted_records_list = cs.sort_fasta(trimmed_fasta_path)
    sorted_fasta_path = cs.write_to_fasta(directories_list[fasta_index],
                                       sorted_fname,
                                       sorted_records_list)
    #NOTE: ADD DATA  VIS IN BETWEEN EACH ROUND OF VSEARCH

    #--------Vsearch Steps-------------#
    # use rename_consensus_headers to assign the adequate headers to the consensus sequence
 
    # generate the vsearch command
    vsearch_cmd, cluster_path = cs.vsearch_command(vsearch_path, 
                                  sorted_fasta_path, 
                                  directories_list[vsearch_index], 
                                  'first_vsearch', 
                                  similarity_id)
    # run the vsearch command
    st_vsearch = os.system(vsearch_cmd)
    
    for i in range(3):
    '''
    TO DO: remove hard coding here, make sure to keep each iteration of the consensus fastas
    Use NCBI Multiple Sequence Aligner to visualy check how each cluster is doing with each iteration 
    Add the positions where 'N' bases start and end and rename clusters after this for loop 
    Push this to github and have sebastian update to the main branch 
    '''
        similarity_id = similarity_id + 0.05
        # realign the reads to the cluster generated by vsearch
       # generate the minimap2 command
        cluster_align_cmd, sam_path = cs.mm2_cluster_aligner(minimap2_path, cluster_path, trimmed_fasta_path, directories_list[vsearch_index], f'aligned_cluster_{i+1}')

        # Run the minimap2 command
        st = os.system(cluster_align_cmd)
        #Convert the outputted sam file into a bam file 
        sam2bam_cmd, bam_path = cs.sam_to_bam(sam_path, directories_list[vsearch_index], f'aligned_cluster_{i+1}')
        st = os.system(sam2bam_cmd)

        #Filter the bam file for only primary alignments 
        first_filtered_bam = cs.filter_primary_alignments(bam_path)

        reference_counts = cs.weight_generation(first_filtered_bam)
        weighted_fasta = cs.weighted_fasta_gen(cluster_path, reference_counts)

        #Now reperform VSEARCH steps using weighted fasta 

        #need to sort weighted fasta file 
        weighted_sorted_list = cs.sort_fasta(weighted_fasta)
        weighted_sorted_fasta_path = cs.write_to_fasta(
                                                    directories_list[vsearch_index],
                                                    'weighted_sorted',
                                                     weighted_sorted_list)

        vsearch_cmd, cluster_path = cs.vsearch_command(vsearch_path,
                                                       weighted_sorted_fasta_path,
                                                       directories_list[vsearch_index],
                                                      f'weighted_consensus_{i+1}',
                                                      similarity_id) #get rid of hard coded 0.9, 
        st = os.system(vsearch_cmd)

'''
    j, k = n_positions[0]


    # return the path to the polished fasta.
    return cs.rename_consensus_headers(medak_cons_path, j, k , lab_cons_path)
'''

def main():
    #in_w_dir = input("Please provide working directory path: ")
    #in_r_dir = input("Please provide read directory path: ")
    #in_f_dir = input("Please provide reference fasta directory path: ")
    
    in_w_dir = sys.argv[1]
    in_r_dir = sys.argv[2]
    in_f_dir = sys.argv[3]
    
    alignment_path = bc_align(in_w_dir, in_r_dir, in_f_dir)
    forward_alignments, reverse_alignments = decouple(alignment_path)
    print('Xenofind [STATUS] - Forward and Reverse strands decoupled')
    
    print('Xenofind [STATUS] - Preparing to create consensus fastas') 
    consensus_generation(in_w_dir, alignment_path, 'all') #full dataset 
    consensus_generation(in_w_dir, forward_alignments, 'forward')
    consensus_generation(in_w_dir, reverse_alignments, 'reverse') 
    
    """
    set 'consensus_generation' above to a variable to return consensus pathway
    """
    #MAJOR NOTE CURRENTLY CODE CANNOT GENERATE FILES WITH DIFFERENT FILE NAMES SINCE IT USES DEFAULT FILE PATH NAMES IN VSEARCH CONSIDERING CREATING INDIVIDAL DIRECTORIES FOR ALL 3
    
    #All reads 
    #functioncall 1 
    
    #forward 
    # function call 2 
    
if __name__ == '__main__':
    main()
