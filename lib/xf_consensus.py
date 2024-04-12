'''
consensus.py is the start of the Xenovo pipeline from the Xenobiology lab at the University of Washington. This script requires fast5 files from Nanopore Sequencing as inputs and will output a reference fasta called 'consensus.fasta'. This outputted fasta file will be used downstream as the ATGC reference file for the location of potential regions that may contain a modified base.
'''
import os
import random 
import subprocess
from Bio import SeqIO
import pysam
import sys
from xf_tools  import *
from xf_params import *

basecall = False
trim = False
sort = False
clustering_VSEARCH = False
#cluster_filter = False
medaka_consensus = True

working_dir = os.path.expanduser(sys.argv[1])
raw_dir = os.path.expanduser(sys.argv[2])
ref_fasta = os.path.expanduser(sys.argv[3])

# Making directories 
working_dir = check_make_dir(working_dir)
ref_dir = check_make_dir(os.path.join(working_dir, 'ref'))
bc_dir = check_make_dir(os.path.join(working_dir, 'basecall'))
fastq_dir = check_make_dir(os.path.join(working_dir, 'fastq/'))
bam_dir = check_make_dir(os.path.join(working_dir, 'bam/'))
pod_dir = check_make_dir(os.path.join(working_dir, 'pod5/'))
processing_dir = check_make_dir(os.path.join(working_dir, 'processing/'))
output_dir = check_make_dir(os.path.join(working_dir, 'outputs/'))

#Parameters (move to xr_tools later) 
min_length = 70
max_length = 300 #MAXximum length reads should be, filters out super long reads that are messing with cluster formation
sample_cluser_size = 5000
min_cluster_seq = 100 #minimum number of reads that need to be in a cluster to be allowed for consensus formation

#Step 0: Generate or merge pod5 files if needed
file_type = os.listdir(raw_dir)
# Check if the directory is not empty
if file_type:
    # Get the first file in the directory
    first_file = file_type[0]
    # Check if the first file is a .pod5 file
    if first_file.endswith(".pod5"):
        if os.path.isfile(os.path.join(pod_dir,os.path.basename(raw_dir))+'.pod5')==False:
            pod5_merge(get_pod5_subdir(raw_dir), os.path.join(pod_dir,os.path.basename(raw_dir))+'.pod5')
        else:
            print('XenoFind [STATUS] - POD5 files merged. Skipping merge')
    else:
        if os.path.isfile(os.path.join(pod_dir,os.path.basename(raw_dir))+'.pod5')==False:
            cod5_to_fast5(get_fast5_subdir(raw_dir), os.path.join(pod_dir,os.path.basename(raw_dir))+'.pod5')
        else: 
            print('XenoFind [STATUS] - Converted POD5 file for modified base found. Skipping POD5 coversion')
else:
    print('XenoFind [ERROR] - Fast5/POD5 directory empty, please check input directory')
    sys.exit()

#Step 1: Basecalling
if basecall == True:
    print('XenoFind [STATUS] - Performing basecalling using Dorado')
    #cmd = os.path.expanduser(basecaller_path)+ ' basecaller hac  --no-trim  ' + pod_dir + ' > '+os.path.join(bc_dir, 'bc.bam') + ' --reference ' + ref_fasta #can probably do this in a bam file as well 
    cmd = os.path.expanduser(basecaller_path)+ ' basecaller hac  --no-trim  --emit-fastq ' + pod_dir + ' > '+os.path.join(bc_dir, 'bc.fq')
    os.system(cmd)
    
    cmd = 'minimap2 -ax map-ont --score-N 0 --MD --min-dp-score 10 '+ref_fasta+' '+os.path.join(bc_dir, 'bc.fq')+ ' > ' +os.path.join(processing_dir,'bc_aligned.sam')
    os.system(cmd)
    
#Step 2: Read Trimming and fasta generation
if trim == True:     
    #CIGAR parser function 
    def read_trim(sam_file_path, fasta_output):
        """
        read_trim is a function that takes in the sam file generated after minimap2 as a string path input & a desired 
        fasta output path and outputs a fasta file to be used downstream. This function iterates through the sam file using pysam 
        scanning for aligned reads that are at least 95% aligned to the dummy reference file, removes soft clipped bases (bases that show up in the sequence but aren't aligned)
        and writes them to the inputted fasta file path. These reads are effectively 'trimmed.' Any reads that are unaligned or less than 95% aligned are ignored. 
        The outputted fasta file contains sequences that are effeictively the same length and is used downstream to create sequence clusters using VSEARCH. 
        """
        # Open the SAM/BAM file
        with pysam.AlignmentFile(sam_file_path, "r") as samfile, open(fasta_output, 'w') as fasta_file:
            # Iterate over each read in the SAM/BAM file
            for read in samfile.fetch():
                if read.is_unmapped:
                    # Skip unmapped reads
                    continue

                # Get the total length of the reference sequence for this read
                reference_length = samfile.get_reference_length(read.reference_name)

                # Ensure the read has an alignment (is not unmapped)
                if not read.is_unmapped:
                    # Compare the aligned portion of the read to the total reference length
                    aligned_length = read.reference_length  # Length of the alignment on the reference
                    if aligned_length < reference_length * 0.95: # If the read is shorter than 95% of the total reference length, toss it
                        print(f"Read {read.query_name} significantly shorter than the reference: ({aligned_length} < {reference_length}). Discarding")
                        continue

                    else:
                        # If the read is not shorter, write to fasta here
                        fasta_file.write(f">{read.query_name}\n{read.query_alignment_sequence}\n")

    read_trim(os.path.join(processing_dir, 'bc_aligned.sam'), os.path.join(processing_dir, 'trimmed.fasta'))
    
#Step 3: Sorting by length 
if sort == True:
    #using biopython
    def sort_fasta(input_fasta, output_fasta):
        records = list(SeqIO.parse(input_fasta, "fasta"))

        sorted_records = sorted(records, key=lambda x: len(x.seq)) #chatGPT first passs, change variable names

        with open(output_fasta, "w") as output_file:
            for record in sorted_records:
                output_file.write(f">{record.id}\n{record.seq}\n")
                
    sort_fasta(os.path.join(processing_dir, 'trimmed.fasta'), os.path.join(processing_dir, 'sorted.fasta'))
    
#Step 4: VSEARCH Clustering 
if clustering_VSEARCH == True: 
    print('XenoFind [STATUS] - Clustering reads with VSEARCH')
    
    # first round of VSEARCH
    cmd = 'vsearch --cluster_fast ' + os.path.join(processing_dir, 'sorted.fasta') + ' --id 0.85 --clusterout_sort --consout ' + os.path.join(processing_dir, 'vs_cons_85.fasta')
    os.system(cmd)
    
    # second round of VSEARCH
    cmd = 'vsearch --cluster_fast ' + os.path.join(processing_dir, 'vs_cons_85.fasta') + ' --id 0.90 --clusterout_sort --consout ' + os.path.join(processing_dir, 'vs_cons_90.fasta')
    os.system(cmd)
    
    # third round of VSEARCH
    cmd = 'vsearch --cluster_fast ' + os.path.join(processing_dir, 'vs_cons_90.fasta') + ' --id 0.95 --clusterout_sort --consout ' + os.path.join(processing_dir, 'vs_cons_95.fasta')
    os.system(cmd)
    
    #cons.fasta made a consensus that was 98.34% similar to GT reference file
    
#Step 5: setting a threshold for amount of reads needed in a cluster for cluster to be considering for consensus sequence formation 
'''
if cluster_filter == True: 
    print('Xenovo [STATUS] - Filtering clusters and choosing representative sequences')

    #need to fix this function, take in a threshold for how many reads in a cluster and filter out, have minimum amount be variable so when diversity is introduced, threshhold should be decreased 
    def cluster_size_filter(input_fasta, output_fasta, threshold):
        filtered_records = []
        for record in SeqIO.parse(input_fasta, "fasta"):
            parts = record.description.split(';')
            if len(parts) > 1 and parts[1].startswith('seqs='):
                size = int(parts[1].split('=')[1])
                if size > threshold:
                    filtered_records.append(f">{record.id}\n{record.seq}")
        with open(output_fasta, "w") as output_file:
            output_file.write("\n".join(filtered_records))
        return output_fasta

    cluster_size_filter(os.path.join(processing_dir, 'cons.fasta'), os.path.join(processing_dir, 'represented_seq.fasta'), min_cluster_seq)
    print('Xenovo [STATUS] - Represented Clusters outputted in', os.path.join(processing_dir, 'represented_seq.fasta'))
'''

#Step 6: Medaka Consensus Sequence Formation 
if medaka_consensus == True: 
    print('XenoFind [STATUS] - Performing consensus sequence formation with Medaka')
    cmd = 'medaka_consensus -i ' + os.path.join(processing_dir, 'trimmed.fasta') + ' -d ' + os.path.join(processing_dir, 'vs_cons_95.fasta') + ' -o ' + output_dir + ' -m r1041_e82_400bps_hac_v4.2.0 -f -b 300' # used with cluster filter
    os.system(cmd)
    
    def extract_n_indexes(n_fasta_file):
        """
        Extract the indexes j and k for each sequence in the fasta file with N regions.
        j: the index just before the first N (0 indexed)
        k: the index of the first N (0 indexed)
        """
        n_indexes = []
        for record in SeqIO.parse(n_fasta_file, "fasta"):
            sequence_str = str(record.seq)
            first_n_index = sequence_str.find('N')
            if first_n_index != -1:
                # Index just before the first N and index of the first N
                j = first_n_index - 1
                k = first_n_index
            else:
                j = k = None
            n_indexes.append((j, k))
        return n_indexes

    def rename_consensus_headers(consensus_fasta_file, n_indexes, output_file):
        """
        Rename the headers in the consensus FASTA file based on the provided N indexes.
        """
        with open(output_file, 'w') as outfile:
            for i, record in enumerate(SeqIO.parse(consensus_fasta_file, "fasta"), start=1):
                if i-1 < len(n_indexes):
                    j, k = n_indexes[i-1]
                    record.description = f"consensus {i}- BC 1: {j}, BC 2: {k}"
                    record.id = f"consensus_{i}"
                    SeqIO.write(record, outfile, "fasta")
                else:
                    print(f"No N index information for consensus {i}. Skipping.")

    n_positions = extract_n_indexes(ref_fasta)
    rename_consensus_headers(os.path.join(output_dir, 'consensus.fasta'), n_positions, os.path.join(output_dir,'labelled_consensus.fasta'))