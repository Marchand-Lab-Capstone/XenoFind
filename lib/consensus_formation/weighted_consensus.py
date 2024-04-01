from Bio import SeqIO
import os
import pysam
import sys
import raw_read_merger as rrm
import setup_methods as setup
import consensus_methods as cs
from .. import xf_params # <-- -this is one directory up, figure this out yujia 


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
    #             all_reads_directory, 6
    #             all_reads_fasta_directory, 7
    #             merged_pod5, 8
    #             rough_consensus_output, 9
    #             xf_consensus_output 10
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
    
    #add another conditional here to stop or force minimap2 alignment 
    map2refcmd = cs.map_to_reference(minimap2_path,
                                  barcoded_fasta,
                                  basecalled_path,
                                  directories_list[1],
                                  aligned_barcode_fname)
    # Run the minimap2 command
    st = os.system(map2refcmd)
    sam_file_path = '{}{}.sam'.format(directories_list[1],aligned_barcode_fname)
    return sam_file_path

def decouple(sam_file_path):
    """
    decouple takes in a fasta file generated from minimap2 and separates 
    it by forward and reverse strand reads using samtools. 
    
    Parameters: 
    sam_file_path: path to the sam file, inputted as a string,
    
    Returns: 
    forward and reverse strand  sam files as a string. ALso generates these files
    in a directory. 
    """
    forward_out_path = os.path.join(os.path.dirname(sam_file_path), 'forward.sam')
    reverse_out_path = os.path.join(os.path.dirname(sam_file_path), 'reverse.sam')
    #This section actually generates the forward and reverse only reads 
    # Open the input SAM file for reading
    with pysam.AlignmentFile(sam_file_path, "r") as infile:

        # Open two files for writing: one for forward strand reads, another for reverse strand reads
        with pysam.AlignmentFile(forward_out_path, "w", header=infile.header) as outfile_forward, \
             pysam.AlignmentFile(reverse_out_path, "w", header=infile.header) as outfile_reverse:

            # Iterate through reads in the input file
            for read in infile:
                # Check if the read is mapped to the reverse strand
                if read.is_reverse and not read.is_unmapped and not read.is_secondary and not read.is_supplementary:
                    # Write the read to the reverse strand reads file
                    outfile_reverse.write(read)
                elif read.is_forward and not read.is_unmapped and not read.is_secondary and not read.is_supplementary:
                    # Otherwise, write the read to the forward strand reads file
                    outfile_forward.write(read)
    return forward_out_path, reverse_out_path 
    

def consensus_generation(working_dir, sam_file_path):

    #Default File names 
    trimmed_fname = 'trimmed' # Fasta 
    sorted_fname = 'sorted' # Fasta
    vsearch_path = 'vsearch' # should be variable
    vs_cons_fname = 'cons' # Fasta
    filtered_fname = 'represented_seq' # Fasta
    medaka_path = 'medaka_consensus' # should be variable
    polished_fname = 'polished_consensus' # Fasta

    #----------Setup----------------------# <----- doing resetup since this list isnt global, maybe ask sebastian if making it a global variable is a bad idea 
    # Set up the working directory 0
    #             basecall_directory, 1
    #             fasta_directory, 2
    #             merged_pod5, 3
    #             rough_consensus_output, 4
    #             xf_consensus_output 5
    directories_list = setup.setup_directory_system(working_dir)

    #--------Trimming and Sorting Steps----------#

    # trim down the samfile to a trimmed fasta using default of 95% margin
    read_trims_list = cs.read_trim(sam_file_path)
    trimmed_fasta_path = cs.write_to_fasta(directories_list[2],
                                        trimmed_fname,
                                        read_trims_list)

    # Sort the trimmed fasta, write it out.
    sorted_records_list = cs.sort_fasta(trimmed_fasta_path)
    sorted_fasta_path = cs.write_to_fasta(directories_list[2],
                                       sorted_fname,
                                       sorted_records_list)
    #NOTE: ADD DATA  VIS IN BETWEEN EACH ROUND OF VSEARCH
'''
    #--------Vsearch Steps-------------#
    # Generate and use vsearch on the fasta, 3 rounds from 85 to 95%.

    for i in range(3):

        # Get the degree to be estimated by loop iteration
        degree = 85 + i*5

        # Get a string version of the degree
        strdeg = str(degree)

        # Get the previous path of the data, default is sorted_fasta_path.
        prevpath = sorted_fasta_path

        # If the iteration is not zero, 
        if (i != 0):

            # The previous path is actually dependant on previous degree.
            prevpath = directories_list[2] + 'filtered{}'.format(degree-5) + '.fasta'

        # generate the vsearch command with the current degree and previous path
        vscmd = cs.vsearch_command(vsearch_path,
                                prevpath,
                                directories_list[2],
                                vs_cons_fname+strdeg,
                                degree/100)

        # run the command.
        st = os.system(vscmd) #asking why os.system is being assigned to a variable


        # from this vsearch, sort it. 
        subsearch_sort = cs.sort_fasta(directories_list[2] + vs_cons_fname + strdeg + '.fasta')

        # Output the sorted one to fasta.
        subsearch_sorted_fasta_path = cs.write_to_fasta(directories_list[2],
                                                     'subsort'+strdeg,
                                                      subsearch_sort)
        if (i!=2):
            # Output the subsearch cluster 
            subsearch_cluster = cs.filter_cluster_size(directories_list[2] + 'subsort{}.fasta'.format(strdeg))
    
            # Write it to a fasta.
            clustered_fasta_path = cs.write_to_fasta(directories_list[2],
                                                 'filtered' +strdeg,
                                                  subsearch_cluster)
    
    # filepath string for the final vsearch-sort-filter fasta
    vsearch_cons_path = directories_list[2] + 'subsort95'+ '.fasta'

    #--------------Run Medaka ---------------------#
    # Generate the medaka command to perform consensus
    mdkacmd = cs.medaka_consensus_command(medaka_path,
                                          trimmed_fasta_path,
                                          vsearch_cons_path,
                                          directories_list[4])
    st = os.system(mdkacmd)
    
    
    # filepath string for the medaka consensus:
    medak_cons_path = directories_list[4] + 'consensus.fasta'
    lab_cons_path = directories_list[4] + 'labeled_consensus.fasta'

    # use rename_consensus_headers to assign the adequate headers to the consensus sequence
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
    consensus_generation(in_w_dir, alignment_path) #full dataset 
    #MAJOR NOTE CURRENTLY CODE CANNOT GENERATE FILES WITH DIFFERENT FILE NAMES SINCE IT USES DEFAULT FILE PATH NAMES IN VSEARCH CONSIDERING CREATING INDIVIDAL DIRECTORIES FOR ALL 3
    
    
    
if __name__ == '__main__':
    main()
