# Components

## Big Idea/Script Components 

### Command Line 
This would be the main UI for users. A user would call ```python xenofind.py [method] [fast5/pod5 directory] [reference fasta]``` with required inputs. Progress and errors would show up here.

### xenofind_pipe.py 
This will be the other 'UI' for users. Users would input where they want the working directory to be, as well as the directories for inputted fast5/pod5 and the reference fasta for basecalling. In addition, users would set which detection method they would want to use on this script. This is essenetially pretty much the same as running the command line as mentioned above except it those inputs would be stored and easy to see. 

### xenofind.py
This is the primary script. This will call the desired detection method in the ```lib``` directory regard if it is called using the command line or xenofind_pipe.py. 

### xf_params.py
This will contain the parameters used in any of the detection methods such as the path for the basecaller. It also contains toggleable variables to only run certain parts of the script (e.g., initial basecall). This is imported into most of the scripts in this program. 

### setup_methods.py
This script contains the functions to set up the directories for the detection method that the users would want to use. This is imported in the main scripts in this program.

### weighted_consensus_split.py
This script takes in the fast5 files generated from Nanopore Sequencing, splits the forward and reverse reads, and outputs a fasta file containing the consensus sequences of the forward and reverse reads respectively as the reference sequences for downstream possible location detection of the XNAs.

### consensus_methods_split.py
This script contains all relevant functions for consensus formation to occur. 

### xf_model_training.py
This is the main script to call submethods to generate features from training dataset (forward/reverse respectively) for training the model to detect XNA position after consensus sequences are built.

### xf_xna_find.py
This script splits the reads based on the direction (forward/reverse) and calls all the submethods to generate features from the data that the users provided. The features are going to be fed into the trained model and XNA detection will be performed.

## Internal Components/ Functions 

### read_trim(xf_consensus.py)
This function takes in the sam file generated from minimap2, trim the reads according to the length of the reference sequence and outputs the reads with approximately the same length in a fasta file. This is done by truncating reads the bases that are longer than the 'dummy' reference file where we have a length approximation of the dataset. The output of read_trim is a fasta file formatted as the following 

    >readID
    trimmed sequence


### filter_cluster_size(consensus.py)
Duringt the rounds of VSEARCH clustering, this function establishes a threshold and selectively filters the clusters based on the number of reads they contain. Only clusters surpassing this specified threshold are considered eligible for the formation of a consensus sequence. The output of this function is a fasta file containing the chosen clusters with the amount of reads grouped in the header. 

### extract_read_info function (xf_low_qual.py)
This function takes in the bam file generated during Dorado basecalling and iterates through it using pysam to extract the readID, basecalled sequence, start of reference sequence, and the per base quality score. It also calculates the per read quality score by taking the average of the per base quality scores for a given ead.

### xfasta_gen (xf_low_qual.py)
xfasta_gen will take in the regions analyzed from xna_guess and will add those significant low-quality score regions to the consensus fasta header in the estabalished 'xFasta' format from the Marchand Lab. The fasta file would be formatted as shown below.

    >consensus_sequence+XPOS['X':potential_positions]
    sequence

### get_fast5_subdir (xf_tools)
get_fast5_subdir is a function that takes in a fast5 containing directory as its input and extracts that path to be given to other functions such as cod5_to_fast5. The pathway extracted is the 'raw_dir' in the string datatype from either terminal calling xenofind.py or from xenofind_pipe.py

### get_pod5_subdir (xf_tools)
get_pod5_subdir is a function that takes in a pod5 containing directory as its input and extracts that path to be given to other functions such as pod5_merge. The pathway extracted is the xna_raw_dir, dna_raw_dir, or bc_raw pathways in the string datatype from either terminal calling xenofind.py or from xenofind_pipe.py

### check_make_dir (xf_tools)
check_make_dir is a function that takes in a directory in the form of a string to actually generate that directory by using the built in os.makedirs() function. The string that is generated is either the inputted working directory string from the user or other predetermined directories that are put into the working directory. There are added print statements to help users know what directory is being generated. 

### cod5_to_fast5 (xf_tools)
cod5_to_fast5 is a function that calls the pod5 package from ONT to convert and merge a fast5 dataset to a single pod5 file within the desired pod5 directory. This is generally the generated pod5 directory in the working directory chosen by the user. 

### pod5_merge (xf_tools)
pod5_merge is a function that calls the pod5 function from ONT to merge a pod5 dataset into a single pod5 file within the desired pod5 directory. This is generally the generated pod5 directory in the working directory chosen by the user.

### filter_primary_alignments (xf_tools)
filter_primary_alignments is a function that takes in a bam file (generally a merged bam file after basecalling) and uses samtools and pysam to iterate through the bam file. During iteration, reads that are not aligned, secondary aligned, or supplementary are ignored. Read that are primary alignments are written to a new bam file, creating a primary alignments only bam file. Samtools is used to 'index' the bam file, which allows pysam to iterate through it.

### first_consensus 
first_consensus takes in a working directory, pod5/fast5 dataset, and a barcoded rough fasta file and runs the commands needed to generate a first-pass fasta file with no xenobases.

### extract_n_indexes
extract_n_indexes is function that extracts the start and end of the randomer 'N' region. This is used downstream to limit k-mer analysis to a specific area of reads. 

### medaka_consensus_command 
medaka_consensus_comand takes in a file path for Medaka by ONT, a trimmed fasta file, and a consensus clusters from VSEARCH and polishes the data. 

### write_to_fasta 
write_to_fasta takes in an output path and a file name, as well as a list of data to be written and writes that data out as a fasta file at the inputted output path. 

### vsearch_command 
vsearch_command takes in a fasta path, a similarity threshold value, and output_path, and an output file name and generates a command to perform clustering and rough consensus on the fasta file. 

### basecall_commands 
basecall_command generates a command to run the dorado basecaller. 