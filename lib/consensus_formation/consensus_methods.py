"""
consensus_methods.py 
J. Sumabat, N. Lai, S. Peck, Y. Huang, 3/14/2024

consensus_methods.py contains the primary methods involved in forming a first-pass
consensus on a DNA sequence with potential XNA bases. 

basecall_command() - Generates a command to run basecalling on pod5 reads using dorado
map_to_reference() - Generates a command to map a basecalled sequence to a reference fasta
vsearch_command() - Generates a command to run cluster vsearch on a fasta file
medaka_consensus_command() - Generates a command to run medaka's consensus on a fasta file.
read_trim() - Generates a list of reads that are trimmed within 95% of the reference.
sort_fasta() - Generates a list of reads, sorted by length.
filter_cluster_size() - Generates a list of reads filtered to a cluster size -- NONFUNCTIONAL
write_to_fasta() - Is able to take a list of reads generated by the previous three methods, and write them to a fasta file. 
first_consensus() - Runs all of the steps involved with generating a consensus.

Danger Zone: these were added specifically for interplay with xf_low_qual, they could be considered last minute additions, and out of scope.
extract_n_indexes() - Pulls the first and last index of the unknown bases from a fasta reference file.
rename_consensus_headers() - renames the headers of a medaka consensus file to have the first and last
                             indexes of the unknown bases.
                             

"""


from Bio import SeqIO
import os
import pysam
import raw_read_merger as rrm
import setup_methods as setup


def basecall_command(basecaller_path, pod5_path, out_path, out_name):
    """
    Basecall_command generates a command to run the basecaller dorado.
    
    Parameters:
    basecaller_path: path to the basecaller, as str
    pod5_path: path to the pod5 file to be basecalled, as str
    out_path: path to the output directory, as str
    out_name: name the basecalled fq file will be given, as str
    
    Returns:
    a command string.
    """
    # Currently only supports Dorado
    cmd = "{} basecaller hac --no-trim --emit-fastq {} > {}{}.fq".format(basecaller_path, pod5_path, out_path, out_name)
    print('[Basecalling]: Command Generated: "{}"'.format(cmd))
    return cmd


def map_to_reference(mapper_path, reference_path, basecall_path, out_path, out_name):
    """
    map_to_reference generates a command to use minimap2 to align the basecall with
    a reference sequence.
    
    Parameters:
    mapper_path: path to minimap2 as str
    reference_path: path to reference fasta (typically the barcodes) as str
    basecall_path: path to basecalled .fq file, as str
    out_path: path to the output directory, as str
    out_name: name the output sam file will be given, as str
    """
    # Currently only supports minimap2
    cmd = "{} -ax map-ont --score-N 0 --MD --min-dp-score 10 {} {} > {}{}.sam".format(mapper_path, reference_path, basecall_path, out_path, out_name)

    print('[Mapping]: Command Generated: "{}"'.format(cmd))
    return cmd


def read_trim(sam_path):
    """
    read_trim takes in a samfile and returns a list of the reads
    in that sam file, with the query name and alignment sequence,
    so long as they are mapped and >=95% of the reference length.
    
    Parameters:
    sam_path: path to the sam file in question as a string.
    
    Returns:
    list of queries and alignments as a string.
    """
    # Sam path is the basecalled to reference
    output_list = []
    
    # open the alignment sam using pysam, open the fasta path as a fasta file
    with pysam.AlignmentFile(sam_path, 'r') as samfile:
 
        # Set up variables for #unmapped reads & #outside length
        num_unmapped = 0
        num_outside = 0
        
        # for each read in the samfile,
        for read in samfile.fetch():

            # Check that the read is mapped
            if (not read.is_unmapped):

                # Get the reference sequence length and the alignment length
                reference_length = samfile.get_reference_length(read.reference_name)
                aligned_length = read.reference_length

                # if the aligned length is greater or equal to 95% of the reference length,
                if aligned_length >= reference_length * .95:
                    
                    # append it to the output.
                    output_list.append(f">{read.query_name}\n{read.query_alignment_sequence}")
                else:
                    num_outside += 1

            else:
                num_unmapped += 1
    print("[Trimming]: {} unmapped reads, {} short reads removed.".format(num_unmapped, num_outside))
    
    # return the output list
    return output_list


# How to validate the read_trim worked?

def sort_fasta(fasta_path):
    """
    sort_fasta takes in a fasta filepath and sorts it by length.
    
    Parameters:
    fasta_path: a path to a trimmed fasta file.

    
    Returns:
    a list of the records, sorted by length.
    """
    # get SeqRecord iterator as list by parsing the fasta file and passing fasta format.
    records = list(SeqIO.parse(fasta_path, "fasta"))
    
    # Sort the records by length using in-built sorting by sequence length.
    sorted_records = sorted(records, key=lambda x: len(x.seq))
    
    # create an empty list to hold the records to be output.
    output_records = []

    for record in sorted_records:
        output_records.append(">{}\n{}".format(record.id, record.seq))

    print("[Sorting]: {} records resorted.".format(len(output_records)))
    
    # return the sorted_records list
    return output_records


def vsearch_command(vsearch_path, fasta_path, out_path, out_name, sim_id):
    """
    vsearch_command takes in a fasta path, a similarity value, an output_path, and an output filename,
    and generates a command to perform clustering/ rough consensus formation on the fasta file.
    
    Parameters:
    vsearch_path: path to vsearch as a string.
    fasta_path: path to the fasta file in question, as a string.
    out_path: path to output the final consensus fasta.
    out_name: filename for the fasta file. 
    sim_id: float value representing how similar the centroids can be (between 0 and 1)
    
    Returns: 
    command to perform vsearch with given directories as a string.
    """
    
    # Generate the command.
    cmd = "{} --cluster_fast {} --id {} --clusterout_sort --consout {}{}.fasta".format(vsearch_path, fasta_path, sim_id, out_path, out_name)

    print('[Vsearching]: Command Generated: "{}"'.format(cmd))
    return cmd
            
    
def filter_cluster_size(fasta_path, threshold=1):
    """
    filter_cluster_size takes in a fasta filepath and a threshold,
    and removes all clusters that are not within that size threshold.
    
    Parameters:
    fasta_path: path to clustered/consensus fasta file, as string.
    threshold: int value representing minimum size. Default is 1.
    
    Returns:
    a list of strings containing the record id and sequence.
    """
    # create a list to store filtered records
    filtered_records = []
    
    # parse through the passed fasta's records
    for record in SeqIO.parse(fasta_path, "fasta"):
    
        # split the description of the record by semicolon
        parsed_record = record.description.split(';')
        
        # Check that the parsed record contains more than one part, get the second segment,
        # and check it starts with the string 'seqs'
        if len(parsed_record) > 1 and parsed_record[-1].startswith('seqs='):
            
            # Get the size of the sequence by the last value
            size = int(parsed_record[-1].split('=')[-1])
            
            # if the size is larger than the threshold,
            if size > threshold:
                
                # add it to the filtered records.
                filtered_records.append(">{}\n{}".format(record.id, record.seq))
    
    return filtered_records
    
    
def write_to_fasta(out_path, out_name, list_data):
    """
    write_to_fasta takes in an output path and file name, as well as a list
    of data to be written, and writes that data out as a fasta file at that
    path.
    
    Parameters:
    out_path: filepath to output directory, as str.
    out_name: name for the fasta file, as str.
    list_data: a list of the data to be written, formatted with each value
               being ">{id}\n{sequence}"
    
    Returns:
    the final path to the output fasta.
    """
    # create the output filename.
    out_file = out_path + out_name + ".fasta"

    if os.path.exists(out_file):
        os.remove(out_file)
    
    # open the output file in write mode.
    with open(out_file, 'w') as output_file:

        for datum in list_data:
            # write each value in the fasta to a new line. 
            output_file.write("{}\n".format(datum))
        print("[Writing Fasta]: {} lines written to {}".format(len(list_data), out_file))
        
    # return the output filepath.
    return out_file
        

def medaka_consensus_command(medaka_path, trim_fasta, filtered_fasta, out_path):
    """
    medaka_consensus_command takes in a filepath for medaka, a trimmed fasta file,
    a filtered cluster/consensus fasta, an output filepath, and then polishes/forms a consensus fasta.
    
    Parameters:
    medaka_path: path, as str, to medaka_consensus
    trim_fasta: path, as str, to trimmed fasta file
    filtered_fasta: path, as str, to filtered fasta file
    out_path: path, as str, to output directory
    """

    # Generate the command.
    cmd = "{} -i {} -d {} -o {} -m r1041_e82_400bps_hac_v4.2.0 -f -b 300".format(medaka_path, trim_fasta, filtered_fasta, out_path)
    print('[Consensus Forming]: Command Generated: "{}"'.format(cmd))
    return cmd


def extract_n_indexes(n_fasta_file):
    """
    extract_n_indexes Extracts the indexes j and k for each sequence in a barcode fasta
    file with unknown regions.
    j: the index of the first N (0 indexed)
    k: the index of the last N (0 indexed)

    Designed specifically for interplay with xf_lowqual

    Parameters:
    n_fasta_file: path to the barcode/reference fasta file

    Returns:
    the indexes, as a list containing touples???? Why:TODO
    
    """
    n_indexes = []
    for record in SeqIO.parse(n_fasta_file, "fasta"):
        sequence_str = str(record.seq)
        first_n_index = sequence_str.find('N')
        last_n_index = sequence_str.rfind('N')
        if first_n_index != -1:
            j = first_n_index
            k = last_n_index
        else:
            j = k = None
        n_indexes.append((j, k))
    return n_indexes


def rename_consensus_headers(consensus_fasta_file, j, k, output_file):
    """
    rename_consensus_headers will
    Rename the headers in the consensus FASTA file based on the provided first and last N indexes.
    All headers will be renamed using the same j and k values.

    Designed specifically for interplay with xf_lowqual

    Parameters:
    consensus_fasta_file: the consensus file output by medaka as consensus.fasta
    j: index of first unknown base
    k: index of last unknown base
    output_file: the filepath to the desired output file

    TODO:NEEDS TO BE BROUGHT UP TO SPEED WITH OTHER WRITING METHDS

    Returns:
    Path to the output file

    
    """
    with open(output_file, 'w') as outfile:
        for i, record in enumerate(SeqIO.parse(consensus_fasta_file, "fasta"), start=1):
            record.description = f"consensus {i}- BC 1: {j}, BC 2: {k}"
            record.description = f" - BC 1: {j}, BC 2: {k}"
            record.id = f"consensus_{i}"
            SeqIO.write(record, outfile, "fasta")

    return str(os.path.abspath(output_file))


def first_consensus(working_dir, reads, barcode_fasta):
    """
    first_consensus takes in working directory, reads, and a barcode fasta,
    and runs the steps needed to generate a first-pass consensus with no
    xenobases.
    """
    
    # Defualt filenames:
    p5_fname = "merged"
    dorado_path = "dorado" # Should be variable, this assumes dorado is in user's PATH
    basecall_fname = 'basecall' # fq file
    minimap2_path = 'minimap2' # should be variable, this assumes dorado is in user's PATH
    aligned_bc_fname = 'bc_aligned' # SAM file
    trimmed_fname = 'trimmed' # Fasta 
    sorted_fname = 'sorted' # Fasta
    vsearch_path = 'vsearch' # should be variable
    vs_cons_fname = 'cons' # Fasta
    filtered_fname = 'represented_seq' # Fasta
    medaka_path = 'medaka_consensus' # should be variable
    polished_fname = 'polished_consensus' # Fasta

    # Get the barcode indexes using extract_n_indexes
    n_positions = extract_n_indexes(barcode_fasta)


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
        bccmd = basecall_command(dorado_path,
                                 merged_pod5_path,
                                 directories_list[1],
                                 basecall_fname)
        # Run the basecall command
        st = os.system(bccmd)

    # ensure the barcode is absolute.
    barcode_fasta = str(os.path.abspath(barcode_fasta))
    
    # use minimap2 to align the basecalled to the basecalled fq
    map2refcmd = map_to_reference(minimap2_path,
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
    read_trims_list = read_trim(samfile_path)
    trimmed_fasta_path = write_to_fasta(directories_list[2],
                                        trimmed_fname,
                                        read_trims_list)

    # Sort the trimmed fasta, write it out.
    sorted_records_list = sort_fasta(trimmed_fasta_path)
    sorted_fasta_path = write_to_fasta(directories_list[2],
                                       sorted_fname,
                                       sorted_records_list)

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
        vscmd = vsearch_command(vsearch_path,
                                prevpath,
                                directories_list[2],
                                vs_cons_fname+strdeg,
                                degree/100)

        # run the command.
        st = os.system(vscmd)


        # from this vsearch, sort it. 
        subsearch_sort = sort_fasta(directories_list[2] + vs_cons_fname + strdeg + '.fasta')

        # Output the sorted one to fasta.
        subsearch_sorted_fasta_path = write_to_fasta(directories_list[2],
                                                     'subsort'+strdeg,
                                                      subsearch_sort)
        if (i!=2):
            # Output the subsearch cluster 
            subsearch_cluster = filter_cluster_size(directories_list[2] + 'subsort{}.fasta'.format(strdeg))
    
            # Write it to a fasta.
            clustered_fasta_path = write_to_fasta(directories_list[2],
                                                 'filtered' +strdeg,
                                                  subsearch_cluster)
    
    # filepath string for the final vsearch-sort-filter fasta
    vsearch_cons_path = directories_list[2] + 'subsort95'+ '.fasta'

    #--------------Run Medaka ---------------------#
    # Generate the medaka command to perform consensus
    mdkacmd = medaka_consensus_command(medaka_path,
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
    renamed_consensus = rename_consensus_headers(medak_cons_path, j, k , lab_cons_path)



    #-------Consensus Basecalling ---------
        # This section is to generate alignment with the consensus file so that we have a basecalled, and aligned, consensus file.

    # Filepath string for the basecalled fq 
    basecalled_path_2 = directories_list[4] + 'basecall_consensus.fq'
    
    if not (os.path.exists(basecalled_path_2)):
    # Generate the dorado basecall command  #do we use dorado aligner? Use dorado aligner rather than minimap2. 
        bccmd_2 = basecall_command(dorado_path,
                                 merged_pod5_path,
                                 directories_list[4],
                                 'basecall_consensus')
        # Run the basecall command
        st = os.system(bccmd_2)

    # ensure the barcode is absolute.
    barcode_fasta = str(os.path.abspath(barcode_fasta))
    
    # use minimap2 to align the basecalled to the basecalled fq
    map2refcmd_2 = map_to_reference(minimap2_path,
                                  barcode_fasta,
                                  basecalled_path_2,
                                  directories_list[4],
                                  'bc_consensus_aligned')
    # Run the minimap2 command
    st = os.system(map2refcmd_2)
    
    return map2refcmd_2


def main():
    in_w_dir = input("Please provide working directory path: ")
    in_r_dir = input("Please provide read directory path: ")
    in_f_dir = input("Please provide reference fasta directory path: ")
    
    consensus_path = first_consensus(in_w_dir, in_r_dir, in_f_dir)
    
    print("Consensus fasta located at: {}").format(consensus_path)
    
    
if __name__ == '__main__':
    main()
    
