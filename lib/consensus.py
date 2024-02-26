'''
consensus.py is the start of the Xenovo pipeline from the Xenobiology lab at the University of Washington. This script requires fast5 files from Nanopore Sequencing as inputs and will output a reference fasta called 'consensus.fasta'. This outputted fasta file will be used downstream as the ATGC reference file for the location of potential regions that may contain a modified base.
'''
import os
import random 
import subprocess
from Bio import SeqIO
import pysam
from xf_tools  import *
from xf_params import *

basecall = True
trim = True
merge_fastq = False
merge_bam = False
VSEARCH_filter = False
sort = False
fastq_to_fasta = False
clustering_VSEARCH = False
cluster_filter = False
medaka_consensus = False

#working_dir = '/home/marchandlab/DataAnalysis/Kawabe/231204_Xenovo_FfAME'
#fast5_dir = '/home/marchandlab/DataAnalysis/Sumabat/230608_XPCR_Taq/20230608_1726_MN37138_AOD257_b27f754c/fast5'
#fast5_dir = '/home/marchandlab/DataAnalysis/Kawabe/231204_Xenovo_FfAME/fast5'
#fast5_dir = '/home/marchandlab/DataAnalysis/Sumabat/xenovo_datasets/mixed_fast5'
#make a reference directory with ground truths for single dataset
#ref_fasta = '/home/marchandlab/github/jay/xenovo-js/ref/230308_PZ_Xemora_train.fa'
#ref_dir = '/home/marchandlab/DataAnalysis/Kawabe/231204_Xenovo_FfAME/reference/231204_FfAME_duplexPCR_adapters.fa'

working_dir = '/home/marchandlab/github/jay/capstone/XenoFind/xenofind_test/240226_trimming_work'
raw_dir = '/home/marchandlab/DataAnalysis/Kaplan/raw/fast5/10.4.1/240104_BSn_90mer_xr_train/50fast5' #dataset to start performing trimming on 
ref_fasta = '/home/marchandlab/DataAnalysis/Kaplan/raw/fast5/10.4.1/240104_BSn_90mer_xr_train/reference/BSn_90mer_xr_train_can.fa'
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
max_length = 300#MAXximum length reads should be, filters out super long reads that are messing with cluster formation
sample_cluser_size = 5000
similarity_id = 0.90 #between 0 - 1, want this higher 
min_cluster_seq = 2 #minimum number of reads that need to be in a cluster to be allowed for consensus formation




#Step 0: FASTA to xFASTA conversion
if os.path.isfile(os.path.expanduser(ref_fasta)): 
    print('pass')
    cmd = 'python xm_fasta2x_rc.py '+os.path.expanduser(ref_fasta)+' '+os.path.join(ref_dir,'x'+os.path.basename(ref_fasta))
    os.system(cmd)
else: 
    print('Xemora  [ERROR] - Reference fasta xna file not file. Please check file exist or file path.')
    sys.exit()

'''
referenc in theory should be [constant region] - variable region (NNNN) - [constant region] 
want to trim constant regions out maybe? 
'''
#Step 1: Generate or merge pod5 files if needed
file_type = os.listdir(raw_dir)
# Check if the directory is not empty
if file_type:
    # Get the first file in the directory
    first_file = file_type[0]
    # Check if the first file is a .pod5 file
    if first_file.endswith(".pod5"):
        if os.path.isfile(os.path.join(mod_pod_dir,os.path.basename(raw_dir))+'.pod5')==False:
            pod5_merge(get_pod5_subdir(raw_dir), os.path.join(pod_dir,os.path.basename(raw_dir))+'.pod5')
        else:
            print('Xenovo [STATUS] - POD5 files merged. Skipping merge')
    else:
        if os.path.isfile(os.path.join(pod_dir,os.path.basename(raw_dir))+'.pod5')==False:
            cod5_to_fast5(get_fast5_subdir(raw_dir), os.path.join(pod_dir,os.path.basename(raw_dir))+'.pod5')
        else: 
            print('Xemora [STATUS] - Converted POD5 file for modified base found. Skipping POD5 coversion')
else:
    print('Xenovo [ERROR] - Fast5/POD5 directory empty, please check input directory')
    sys.exit()

#Step 1: Basecalling
if basecall == True:
    #think about ways to get 100% alignmet compared to reference file
    #cmd=os.path.expanduser(basecaller_path)+' -i '+pod_dir+' -s '+fastq_dir+' -c '+guppy_config_file+' -x auto --bam_out --index --moves_out  -a '+os.path.join(ref_dir,'x'+os.path.basename(ref_fasta)) #Guppy basecalling with alignment
    #cmd=os.path.expanduser(basecaller_path)+' -i '+pod_dir+' -s '+ fastq_dir+' -c '+guppy_config_file+' -x auto --bam_out --index'
    print('Xenovo [STATUS] - Performing basecalling using Dorado')
    cmd = os.path.expanduser(basecaller_path)+ ' basecaller hac  --no-trim  ' + pod_dir + ' > '+os.path.join(bc_dir, 'bc.bam') + ' --reference ' + ref_fasta #can probably do this in a bam file as well 
    os.system(cmd)
    
#need a merge bam instead of a merge fastq if i want to use alignment check, neither of these needed if using dorado
if merge_fastq == True:
    print('Xenovo [STATUS] - Merging fastq files')
    cmd = 'find ' + fastq_dir + " -type f -name '*.fastq' -exec cat {} + > " + processing_dir + 'merged.fastq'
    os.system(cmd) 
if merge_bam == True: 
    cmd = 'samtools merge ' + os.path.join(bam_dir, os.path.basename(bam_dir))+'.bam'+' ' +os.path.join(mod_fastq_dir,'pass/*.bam -f')
    os.system(cmd)
#Read Trimming 
if trim == True: 
    print('test')
    
    #Reads are trimmed to have all sequences be approximately the same length, will make error analysis be constant as  the ends would not be significatly adding to error 
    #Trim fastq files maybe? 
    '''
    def read_trim(input_bam, output_fasta):
        bamfile = pysam.AlignmentFile(input_bam, "rb")
        for read in bamfile: 
            read_length = 
            if len(
    '''
    '''
    Pseudo code: 
    for reads 
    check length of read compared to reference at START 
    1. decide if keep reads if length of query is shorter than reference (might be a question for jorge) 
    2. if the read is Loner than the reference sequence go through the followin
    
    Check start of read compared to alignment, reference file will have which base it starts aligning to 
        if there are bases before that point, add those basees to bases kept using indexing till reaching samtool index 1 
        if there arent enough bases 
            ????
        
    at end, not sure how to decide trimming cut off, maybe just go to end of sequence and see if the read length reaches the end 
        if shorter than desired length keep
        if longer trim
        
    consider case where total length of read is < length of reference sequence 
        toss read? keep? threshold? 
    
    append sequences to a new fasta file with query name and sequence
    
    
    '''
#Step 2: Merge fastq 


# Length and Phred score filtering using VSEARCH, thisstep might not be needed anymore?
if VSEARCH_filter == True: 
    cmd = 'vsearch --fastx_filter ' + os.path.join(processing_dir, 'merged.fastq') + ' --fastq_qmax 90 --fastq_maxlen ' + str(max_length) + ' --fastq_minlen ' + str(min_length) + ' --fastqout ' + os.path.join(processing_dir, 'filtered.fastq') 
    os.system(cmd)

if fastq_to_fasta == True:
    print('Xenovo [STATUS] - Convering from fastq to fasta')
    cmd = 'seqtk seq -a ' + os.path.join(processing_dir, 'filtered.fastq') + ' > ' + os.path.join(processing_dir, 'filtered.fasta')
    os.system(cmd)

#Sorting by length 
if sort == True:
    #using biopython
    records = list(SeqIO.parse("your_fasta_file.fasta", "fasta"))
    sorted_records = sorted(records, key=lambda x: len(x.seq))
    SeqIO.write(sorted_records, "sorted_fasta_file.fasta", "fasta")
    
#Step 3: VSEARCH Clustering 
if clustering_VSEARCH == True: 
    print('Xenovo [STATUS] - Clustering reads with VSEARCH')
    # Input FASTA file and the number of sequences to select

    def random_fasta_sample(input_fasta, output_fasta):
        # Read the input FASTA file and store sequences in a list
        with open(input_fasta, "r") as infile:
            fasta_lines = infile.readlines()

        #Create a list containing paired headers and sequences 
        head_seq_pair = [fasta_lines[i] + fasta_lines[i+1] for i in range(0, len(fasta_lines), 2)]
        
        # Randomly select X sequences without replacement
        selected_pairs = random.sample(head_seq_pair, sample_cluser_size)

        # Write the selected sequences to an intermediate FASTA file
        with open(output_fasta, "w") as outfile:
            outfile.writelines(selected_pairs)
    
    #call function for first round of clustering
    random_fasta_sample(os.path.join(processing_dir, 'filtered.fasta'), os.path.join(processing_dir, 'random_sample.fasta'))

    cmd = 'vsearch --cluster_fast ' + os.path.join(processing_dir, 'random_sample.fasta') + ' --id ' + str(similarity_id) + ' --centroids ' + os.path.join(processing_dir, 'clusters.fasta') + ' --uc ' + os.path.join(processing_dir, 'cluster_info.uc') + ' --sizeout --clusterout_sort --consout ' + os.path.join(processing_dir, 'cons.fasta')
    os.system(cmd)

#Step 4: setting a threshold for amount of reads needed in a cluster for cluster to be considering for consensus sequence formation 
if cluster_filter == True: 
    print('Xenovo [STATUS] - Filtering clusters and choosing representative sequences')

    def cluster_size_filter(input_fasta, output_fasta, threshold):
        modified_records = []
        cluster_number = 1
        for record in SeqIO.parse(input_fasta, "fasta"):
            parts = record.description.split(';')
            if len(parts) > 1 and parts[1].startswith('size='):
                size = int(parts[1].split('=')[1])
                if size > threshold:
                    record.id = f"Cluster{cluster_number};size={size}"
                    record.description = ""
                    modified_records.append(record)
                    cluster_number += 1
        SeqIO.write(modified_records, output_fasta, "fasta")

    cluster_size_filter(os.path.join(processing_dir, 'clusters.fasta'), os.path.join(processing_dir, 'represented_seq.fasta'), min_cluster_seq)
    print('Xenovo [STATUS] - Represented Clusters outputted in', os.path.join(processing_dir, 'represented_seq.fasta'))
#Step 4: Medaka Consensus Sequence Formation 
if medaka_consensus == True: 
    print('Xenovo [STATUS] - Performing consensus sequence formation with Medaka')
    cmd = 'medaka_consensus -i ' + os.path.join(processing_dir, 'filtered.fasta') + ' -d ' + os.path.join(processing_dir, 'represented_seq.fasta') + ' -o ' + output_dir + ' -m r1041_e82_400bps_hac_v4.2.0 -f -b 300'
    os.system(cmd)


#down stream, maybe add adapter sequences before inputting fasta into xenomorph? 
