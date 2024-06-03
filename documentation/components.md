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

### raw_read_merger.py
This script contains methods used for taking a directory of either fast5 or pod5 files, and then merging them into a singular pod5 file.

### weighted_consensus_split.py
This script takes in the fast5 files generated from Nanopore Sequencing, splits the forward and reverse reads, and outputs a fasta file containing the consensus sequences of the forward and reverse reads respectively as the reference sequences for downstream possible location detection of the XNAs.

### consensus_methods_split.py
This script contains all relevant functions for consensus formation to occur. 

### xf_model_training.py
This is the main script to call submethods to generate features from training dataset (forward/reverse respectively) for training the model to detect XNA position after consensus sequences are built.

### xf_xna_find.py
This script splits the reads based on the direction (forward/reverse) and calls all the submethods to generate features from the data that the users provided. The features are going to be fed into the trained model and XNA detection will be performed.

### xf_basecall.py
This is a script imported in ```xf_model_training.py``` and ```xf_xna_find.py``` to perform basecalling using Dorado and align the basecalled reads to a reference fasta by minimap2.

### feature_extraction.py
This script contains methods useful for extracting possible machine learning features from a consensus of reads and their aggregate information as produced by ```aggregate_reads.py```. This is designed to operate on one json file at a time, and is hitherto incomplete with respect to actually extracting all features. Presently, this includes base-level and quality-score-level feature extraction, with signal extraction in escrow.

### aggregate_reads.py
This script is designed around importing a bamfile of reads, a consensus fasta, and a pod5 file, 
and then combining the information contained into separate json files organized by consensus sequence and constituent reads. This was penned with the hope of avoiding/mitigating memory leakage from a variety of sources. 

### xr_fasta2x_rc.py

### xf_rc_fasta_gen.py 
This script takes in the input fasta file and outputs the reverse complement of the fasta file.

## Internal Components/ Functions 

### check_make_dir (setup_methods.py)
check_make_dir is a function that takes in a directory in the form of a string to actually generate that directory by using the built in os.makedirs() function. The string that is generated is either the inputted working directory string from the user or other predetermined directories that are put into the working directory. There are added print statements to help users know what directory is being generated.

### setup_directory_system_[con/model/find] (setup_methods.py)
setup_directory_system takes in a working directory filepath for the chosen method, and ensures that the required subdirectories are extant.

### contains_xna_bases (setup_methods.py)
contains_xna_bases checks if a sequence has xna bases.

### merge_reads_command (raw_read_merger.py)
merge_reads_command merges the subcommands to generate a full command of converting all the pod5 or fast5 in the reads directory to a merged pod5 file.

### validate_read_directory (raw_read_merger.py)
validate_read_directory checks the existence of the reads directory.

### validate_target_directory (raw_read_merger.py)
validate_target_directory checks the existence of the target directory.

### generate_merge_pod5 (raw_read_merger.py)
generate_merge_pod5 takes in the reads directory and target directory, runs the command generated by the ```merge_read_command``` and merges the pod5/fast5 files.

### consensus_generation (weighted_consensus_split.py)
consensus_generation takes in working directory, reads, and a barcoded fasta file and performs basecalling on them using Dorado. It then performs alignment using minimap2 set with custom parameters to maximize the amount of reads that align to the initial placeholder fasta, filters with the primary aligned reads, and trims and sorts the read using default of 95% margin. The trimmed and sorted reads are then fed into vsearch clustering and the outcome consensus of each cluster is weighted by the number of reads that gets aligned to. Realign the reads to the weighted consensus and reperform the vsearch steps until the final consensus is generated.

### basecall_commands (consensus_methods_split.py)
basecall_command generates a command to run the dorado basecaller.

### map_to_reference (consensus_methods_split.py)
map_to_reference generates a command to use minimap2 to align the basecall with a reference sequence.

### filter_primary_alignments (consensus_methods_split.py)
Filters a BAM file to retain only primary alignments and outputs the result to a new BAM file named 'primary_alignments.bam' in the same directory as the input file.

### read_trim (consensus_methods_split.py)
This function takes in the sam file generated from minimap2, trim the reads according to the length of the reference sequence and outputs the reads with approximately the same length in a fasta file. This is done by truncating reads the bases that are longer than the 'dummy' reference file where we have a length approximation of the dataset. The output of read_trim is a fasta file formatted as the following 

    >readID
    trimmed sequence

### write_to_fasta (consensus_methods_split.py)
write_to_fasta takes in an output path and a file name, as well as a list of data to be written and writes that data out as a fasta file at the inputted output path. 

### sort_fasta (consensus_methods_split.py)
sort_fasta takes in a fasta filepath and sorts it by length.

### vsearch_command (consensus_methods_split.py)
vsearch_command takes in a fasta path, a similarity threshold value, and output_path, and an output file name and generates a command to perform clustering and rough consensus on the fasta file.

### mm2_cluster_aligner (consensus_methods_split.py)
mm2_cluster_aligner takes in a cluster fasta from VSEARCH as well as the original trimmed data set and performs realignment on the data using minimap2.

### weight_generation (consensus_methods_split.py)
weight_generation takes in an aligned BAM (primary reads only) and returns a dictionary containing the references (clusters) from VSEARCH and the associated count of primary aligned reads.

### weighted_fasta_gen (consensus_methods_split.py)
weighted_fasta_gen will take a cluster fasta outputted from VSEARCH as well as a dictionary containing the reference sequences in that sam file (clusters outputted from VSEARCH) and the number of reads aligned to that particular reference sequence. It will generate a fasta file with the same reference sequence except multiplied by the number of reads that aligned to perfor weight cluster formation downstream. 

### rename_consensus_headers (consensus_methods_split.py)
rename_consensus_headers will rename the headers in the consensus FASTA file based on the provided first and last N indexes. All headers will be renamed using the same j and k values.

### cluster_size_df (consensus_methods_split.py)
cluster_size_df takes in the cluster fasta generated by VSEARCH, extracts the size of each cluster during rounds of VSEARCH and stores the size of each cluster in every iteration in a dataframe for a later csv file generation.

### filter_cluster_size (consensus_methods_split.py)
Duringt the rounds of VSEARCH clustering, this function establishes a threshold and selectively filters the clusters based on the number of reads they contain. Only clusters surpassing this specified threshold are considered eligible for the formation of a consensus sequence. The output of this function is a fasta file containing the chosen clusters with the amount of reads grouped in the header. 

### merge_fasta_files (consensus_methods_split.py)
Merges multiple FASTA files into one.

### preprocessing (xf_model_training.py and xf_xna_find.py)
preprocessing takes in a working directory, the raw data, the reference fasta file and the direction (e.g. forward/reverse) and will generate necessary files (e.g. merged pod5 files and bam files) needed for feature extraction or generation.

### raw_basecall_features (xf_model_training.py and xf_xna_find.py)
raw_basecall_features runs the script to extract features from raw data & sequence space and merge them. Outputs each reference sequence as a json file containing the reads with their merged features. 

### consensus_features (xf_model_training.py and xf_xna_find.py)
consensus_features takes in a the file path containing the json files with merged raw and sequence space data. Calculates consensus features for each of json file. Currently, data is stored in a list of pandas dataframes. 

### split_fasta (xf_xna_find.py)
split_fasta takes in the input consensus fasta files from the users and splits it by direction for the following preprocessing and feature extraction.

### dorado_bc_command (xf_basecall.py)
dorado_bc_command generates a command to run the basecaller dorado.

### alignment_command (xf_basecall.py)
alignment_command generates a command to run the minimap2 alignment.

### reverse_complement_fasta (xf_rc_fasta_gen.py)
reverse_complement_fasta reads the input fasta file and creat the reverse complement.