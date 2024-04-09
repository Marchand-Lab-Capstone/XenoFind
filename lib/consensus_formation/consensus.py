from Bio import SeqIO
import os
import pysam
import raw_read_merger as rrm
import setup_methods as setup
import consensus_methods as cs

def first_consensus(working_dir, reads, barcode_fasta):
    """
    first_consensus takes in working directory, reads, and a barcode fasta,
    and runs the steps needed to generate a first-pass consensus with no
    xenobases.
    """
    
    # Defualt filenames:
    p5_fname = "merged"
    #dorado_path = "dorado" # Should be variable, this assumes dorado is in user's PATH
    dorado_path = ' ~/dorado-0.5.3-linux-x64/bin/dorado' # assumes dorado is in user's home directory, make it a variable somewhere maybe
    basecall_fname = 'basecall' # fq file
    minimap2_path = 'minimap2' #called from conda environment
    aligned_bc_fname = 'bc_aligned' # SAM file
    trimmed_fname = 'trimmed' # Fasta 
    sorted_fname = 'sorted' # Fasta
    vsearch_path = 'vsearch' # should be variable
    vs_cons_fname = 'cons' # Fasta
    filtered_fname = 'represented_seq' # Fasta
    medaka_path = 'medaka_consensus' # should be variable
    polished_fname = 'polished_consensus' # Fasta

    # Get the barcode indexes using extract_n_indexes
    n_positions = cs.extract_n_indexes(barcode_fasta)


    #----------Setup----------------------#
    # Set up the working directory 0
    #             basecall_directory, 1
    #             fasta_directory, 2
    #             merged_pod5, 3
    #             rough_consensus_output, 4
    #             xf_consensus_output 5
    directories_list = setup.setup_directory_system(working_dir)

    #File path string for the merged pod5
    merged_pod5_path = directories_list[3] + p5_fname + '.pod5'
    
    if not (os.path.exists(merged_pod5_path)):
        # Using RRM, generate the pod5 from the data directory
        rrm.generate_merged_pod5(reads,
                                 directories_list[3],
                                 p5_fname)

    #-------Basecalling and Sorting ---------


    # Filepath string for the basecalled fq 
    basecalled_path = directories_list[1] + basecall_fname + '.fq'
    
    if not (os.path.exists(basecalled_path)):
    # Generate the dorado basecall command 
        bccmd = cs.basecall_command(dorado_path,
                                 merged_pod5_path,
                                 directories_list[1],
                                 basecall_fname)
        # Run the basecall command
        st = os.system(bccmd)

    # ensure the barcode is absolute.
    barcode_fasta = str(os.path.abspath(barcode_fasta))
    
    # use minimap2 to align the basecalled to the basecalled fq
    map2refcmd = cs.map_to_reference(minimap2_path,
                                  barcode_fasta,
                                  basecalled_path,
                                  directories_list[1],
                                  aligned_bc_fname)
    # Run the minimap2 command
    st = os.system(map2refcmd)

    #--------Trimming and Sorting Steps----------#
    # Filepath string for the sam file.
    samfile_path = directories_list[1] + aligned_bc_fname + '.sam'

    # trim down the samfile to a trimmed fasta using default of 95% margin
    read_trims_list = cs.read_trim(samfile_path)
    trimmed_fasta_path = cs.write_to_fasta(directories_list[2],
                                        trimmed_fname,
                                        read_trims_list)

    # Sort the trimmed fasta, write it out.
    sorted_records_list = cs.sort_fasta(trimmed_fasta_path)
    sorted_fasta_path = cs.write_to_fasta(directories_list[2],
                                       sorted_fname,
                                       sorted_records_list)

    #--------Vsearch Steps-------------#
    # Generate and use vsearch on the fasta, 3 rounds from 85 to 95%.
    '''
    Ask Sebastian where primary aligned filtering went
    '''
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


def main():
    in_w_dir = input("Please provide working directory path: ")
    in_r_dir = input("Please provide read directory path: ")
    in_f_dir = input("Please provide reference fasta directory path: ")
    
    consensus_path = first_consensus(in_w_dir, in_r_dir, in_f_dir)
    
    print("Consensus fasta located at: {}").format(consensus_path)
    
    
if __name__ == '__main__':
    main()
