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
import subprocess
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
    cmd = "{} basecaller hac --no-trim --emit-moves --emit-fastq {} > {}{}.fq".format(basecaller_path, pod5_path, out_path, out_name)
    print('[Basecalling]: Command Generated: "{}"'.format(cmd))
    return cmd


def map_to_reference(mapper_path, reference_path, basecall_path, out_path, out_name):
    """
    map_to_reference generates a command to use minimap2 to align the basecall with
    a reference sequence.
    
    Parameters:
    mapper_path: path to minimap2 as str
    reference_path: path to reference fasta as str
    basecall_path: path to basecalled .fq file, as str
    out_path: path to the output directory, as str
    out_name: name the output sam file will be given, as str
    """
    # Currently only supports minimap2
    cmd = "{} -ax map-ont --score-N 0 --MD --min-dp-score 10 {} {} > {}{}.sam".format(mapper_path, reference_path, basecall_path, out_path, out_name) #probably doesn't need minimap2 as a separate path

    print('[Mapping]: Command Generated: "{}"'.format(cmd))
    return cmd

def sam_to_bam(input_sam, out_path, out_name):
    """
    sam_to_bam is a function that takes in a sam file generated from minimap2 
    and generates a samtools command that can be run using os.system to 
    generated a bam file containing the sam information 
    
    Parameters: 
    input_sam: file path way to generated sam file as a string 
    out_path: file path where bam file should be generated 
    out_name: desired name of the bam file 
    """
    cmd = 'samtools view -bho {}{}.bam {}'.format(out_path, out_name, input_sam)
    bam_path = '{}{}.bam'.format(out_path, out_name)
    return cmd, bam_path
    

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
        
#def strand_decouple(primary_sam_path, forward_out_path, reverse_out_path): #hae outpaths be GENERATED
def strand_decouple(primary_bampath):
    """
    strand_decouple takes in a bam file generated & converted using minimap2 and samtools 
    and separates it by forward and reverse strand reads using samtools. 
    
    Parameters: 
    sam_path: path to the sam file as a string,
    
    Returns: 
    forward and reverse strand  sam files as a string. ALso generates these files
    in a directory. 
    
    NOTE: NEED TO EDIT THIS FUNCTION INTO TWO FUNCTIONS, GENERATE PRIMARY SAM FILE PATH STRING AND RUN IT IN THE FIRST PASS FUNCTION. SECOND FUNCTION TO GENERATE THE FORWARD AND REVERSE STRINGS
    """
    forward_out_path = os.path.join(os.path.dirname(primary_sam_path), 'forward.sam')
    reverse_out_path = os.path.join(os.path.dirname(primary_sam_path), 'reverse.sam')
    #This section actually generates the forward and reverse only reads 
    # Open the input SAM file for reading
    with pysam.AlignmentFile(output_sam, "r") as infile:

        # Open two files for writing: one for forward strand reads, another for reverse strand reads
        with pysam.AlignmentFile(forward_out_path, "w", header=infile.header) as outfile_forward, \
             pysam.AlignmentFile(reverse_out_path, "w", header=infile.header) as outfile_reverse:

            # Iterate through reads in the input file
            for read in infile:
                # Check if the read is mapped to the reverse strand
                if read.is_reverse:
                    # Write the read to the reverse strand reads file
                    outfile_reverse.write(read)
                else:
                    # Otherwise, write the read to the forward strand reads file
                    outfile_forward.write(read)
        return forward_out_path, reverse_out_path #maybe dont need to return these but will leave this here for now 

def read_trim(bam_path):
    """
    read_trim takes in a BAM file, sorts and indexes it, then returns a list of the reads
    in that BAM file, with the query name and alignment sequence,
    so long as they are mapped and >=95% of the reference length.
    
    Parameters:
    bam_path: path to the BAM file in question as a string.
    
    Returns:
    list of queries and alignments as a string.
    """
    output_list = []

    # Generate paths for sorted and indexed bam
    sorted_bam_path = bam_path + ".sorted.bam"
    
    # Sorting the bam file
    sort_cmd = f'samtools sort {bam_path} -o {sorted_bam_path}'
    subprocess.run(sort_cmd, shell=True, check=True)
    
    # Indexing the sorted bam file
    index_cmd = f'samtools index {sorted_bam_path}'
    subprocess.run(index_cmd, shell=True, check=True)
    
    # Open the sorted and indexed bam using pysam
    with pysam.AlignmentFile(sorted_bam_path, 'rb') as bamfile:
        num_unmapped = 0
        num_outside = 0
        
        # for each read in the bamfile,
        for read in bamfile.fetch():
            # Check that the read is mapped
            if not read.is_unmapped:
                # Get the reference sequence length and the alignment length
                reference_length = bamfile.get_reference_length(read.reference_name)
                aligned_length = read.reference_length

                # if the aligned length is greater or equal to 95% of the reference length,
                if aligned_length >= reference_length * 0.95:
                    # Append it to the output.
                    output_list.append(f">{read.query_name}\n{read.query_alignment_sequence}")
                else:
                    num_outside += 1
            else:
                num_unmapped += 1

    print(f"[Trimming]: {num_unmapped} unmapped reads, {num_outside} short reads removed.")
    
    # Remove the sorted and indexed BAM file after processing
    os.remove(sorted_bam_path)
    # Also remove the BAM index file
    os.remove(sorted_bam_path + ".bai")
    
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
    cluster_path = "{}{}.fasta".format(out_path, out_name)
    print('[Vsearching]: Command Generated: "{}"'.format(cmd))
    return cmd, cluster_path
            
    
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

def mm2_cluster_aligner(mapper_path, reference_path, trimmed_reads, out_path, out_name):
    """
    mm2_cluster_aligner takes in a cluster fasta from VSEARCH as well as the 
    original trimmed data set and performs realignment on the data using 
    minimap2
    
    Parameters:
    mapper_path: path to minimap2 as str
    reference_path: path to cluster/consensus reference fasta
    trimmed_reads: path to basecalled and trimmed reads in fasta format
    out_path: path to the output directory, as str
    out_name: name the output sam file will be given, as str
    
    Returns:
    cmd: command with apppropriate minimap2 parameters and inputs
    """
    #NOTES FOR SELF: should probably have it generate the minimap2 string and not run it??? 
    #other note: since sam file is outputted, can use one of the methods above to filter it, then use the length of the sam file (assuming no headers are outputted) to get the weight of the particular cluster 
    cmd = "{} -ax map-ont --MD {} {} > {}{}.sam".format(mapper_path, reference_path, trimmed_reads, out_path, out_name) #probably doesn't need minimap2 as a separate path
    sam_path = "{}{}.sam".format(out_path, out_name)
    
    print('[Mapping]: Command Generated: "{}"'.format(cmd))
    return cmd, sam_path

def weight_generation(aligned_bam_path):
    """
    Takes in an aligned BAM (primary reads only) and returns a dictionary containing 
    the references (clusters) from VSEARCH and the associated count of primary 
    aligned reads.
    
    Parameters: 
    aligned_bam_path: File path to inputted BAM file generated from mm2_cluster_aligner.
    
    Returns: 
    reference_counts: Dictionary containing the different reference sequence names 
    and the count of reads that aligned to each.
    """
    # Open the BAM file for reading ("rb" mode)
    with pysam.AlignmentFile(aligned_bam_path, "rb") as bamfile:
        # Initialize a dictionary to hold the count of reads per reference sequence
        reference_counts = {}
        
        # Iterate over each read in the BAM file
        for read in bamfile:
            # Check if the read is mapped and is a primary alignment (not secondary/supplementary)
            if not read.is_unmapped and not read.is_secondary and not read.is_supplementary:
                ref_name = bamfile.get_reference_name(read.reference_id)  # Get the reference sequence name
                
                if ref_name in reference_counts:
                    # If the reference sequence is already in the dictionary, increment the count
                    reference_counts[ref_name] += 1
                else:
                    # If it's the first time we see this reference sequence, initialize its count to 1
                    reference_counts[ref_name] = 1
                    
    return reference_counts

def weighted_fasta_gen(cluster_fasta, reference_counts, out_name):
    """
    weighted_fasta_gen will take a cluster fasta outputted from VSEARCH
     as well as a dictionary containing the reference 
    sequences in that sam file (clusters outputted from VSEARCH) and the number 
    of reads aligned to that particular reference sequence. It will generate 
    a fasta file with the same reference sequence except multiplied by the number 
    of reads that aligned to perfor weight cluster formation downstream. 
    
    Parameters: 
    cluster_fasta: used to extract the actual sequence of the reference 
    reference_counts: dictionary containing the names of the reference sequences (should be the same names in cluster_fasta) as well as the number of reads that aligned to that cluster after minimap2
    
    Returns:
    weighted_fasta_file_path: fasta file containing the same reference sequences as cluster_fasta except multiplied by the number of counts from reference_counts 
    """
    weighted_fasta_file_path = os.path.join(os.path.dirname(cluster_fasta), '{}.fasta'.format(out_name))
    
    with open(weighted_fasta_file_path, "w") as output_file:
        #write into the weighted fasta file
        for record in SeqIO.parse(cluster_fasta, "fasta"):
            reference_name = record.id
            sequence = str(record.seq)
            # find the reference name in the cluster fasta
            
            if reference_name in reference_counts:
                count = reference_counts[reference_name]
                i = 0
                for i in range(count):
                 # multiply the reference sequence by the number of reads that aligned to it
                    output_file.write(f">{reference_name}_weighted_{i+1}\n{sequence}\n")
                
    return weighted_fasta_file_path

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
        

def medaka_consensus_command(medaka_path, trim_fasta, cluster_fasta, out_path):
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
    cmd = "{} -i {} -d {} -o {} -m r1041_e82_400bps_hac_v4.2.0 -f -b 300".format(medaka_path, trim_fasta, filtered_fasta, out_path) #note, allow specific model selection
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

