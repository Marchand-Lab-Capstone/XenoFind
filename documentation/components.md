# Components

## Big Idea/Script Components 

### Command Line 
This would be the main UI for users. A user would call ```python xenofind.py [method] [fast5/pod5 directory] [reference fasta]``` with required inputs. Progress and errors would show up here.

### xenofind_pipe.py 
```xenofind_pipe.py``` will be the other 'UI' for users. Users would input where they want the working directory to be, as well as the directories for inputted FAST5/POD5 and the reference FASTA for basecalling. In addition, users would set which detection method they would want to use on this script. This is essenetially pretty much the same as running the command line as mentioned above except it those inputs would be stored and easy to see. 

### xenofind.py
```xenofind.py``` is the primary script. This will call the desired detection method in the ```lib``` directory regard if it is called using the command line or xenofind_pipe.py. 

### xf_params.py
```xf_params.py``` will contain the parameters used in any of the detection methods such as the path for the basecaller. It also contains toggleable variables to only run certain parts of the script (e.g., initial basecall). This is imported into most of the scripts in this program. 

### setup_methods.py
```setup_methods.py``` contains the functions to set up the directories for the detection method that the users would want to use. This is imported in the main scripts in this program.

### raw_read_merger.py
```raw_read_merger.py``` contains methods used for taking a directory of either FAST5 or POD5 files, and then merging them into a singular POD5 file.

### weighted_consensus_split.py
```weighted_consensus_split.py``` takes in the fast5 files generated from Nanopore Sequencing, splits the forward and reverse reads, and outputs a FASTA file containing the consensus sequences of the forward and reverse reads respectively as the reference sequences for downstream possible location detection of the XNAs.

### consensus_methods_split.py
```consensus_methods_split.py``` contains all relevant functions for consensus formation to occur. 

### xf_model_training.py
```xf_model_training.py``` is the main script to call submethods to generate features from training dataset (forward/reverse respectively) for training the model to detect XNA position after consensus sequences are built.

### xf_xna_find.py
```xf_xna_find.py``` splits the reads based on the direction (forward/reverse) and calls all the submethods to generate features from the data that the users provided. The features are going to be fed into the trained model and XNA detection will be performed.

### xf_basecall.py
```xf_basecall.py``` is a script imported in ```xf_model_training.py``` and ```xf_xna_find.py``` to perform basecalling using Dorado and align the basecalled reads to a reference FASTA file by minimap2.

### feature_extraction.py
```feature_extraction.py``` contains methods useful for extracting possible machine learning features from a consensus of reads and their aggregate information as produced by ```aggregate_reads.py```. This is designed to operate on one JSON file at a time, and is hitherto incomplete with respect to actually extracting all features. Presently, this includes base-level and quality-score-level feature extraction, with signal extraction in escrow.

### aggregate_reads.py
```aggregate_reads.py``` is designed around importing a bamfile of reads, a consensus FASTA file, and a POD5 file, and then combining the information contained into separate JSON files organized by consensus sequence and constituent reads. This was penned with the hope of avoiding/mitigating memory leakage from a variety of sources. 

### xr_fasta2x_rc.py
```xr_fasta2x_rc.py``` takes in a inputted FASTA file containing an xna (e.g., P),  and generates a file in the xFASTA format as mentioned in the README.md. This file also generates a reverse complemented xFASTA file for the opposite XNA (e.g., Z for P). These xFASTA files are used for model training to establish XNA labels. 

### xf_rc_fasta_gen.py 
```xf_rc_fasta_gen.py``` takes in the input FASTA file and outputs the reverse complement of the FASTA file.

## Internal Components/ Functions 

### check_make_dir (setup_methods.py)
```check_make_dir``` is a function that takes in a directory in the form of a string to actually generate that directory by using the built in os.makedirs() function. The string that is generated is either the inputted working directory string from the user or other predetermined directories that are put into the working directory. There are added print statements to help users know what directory is being generated.

### setup_directory_system_[con/model/find] (setup_methods.py)
```setup_directory_system``` takes in a working directory filepath for the chosen method, and ensures that the required subdirectories are extant.

### contains_xna_bases (setup_methods.py)
```contains_xna_bases``` checks if a sequence has xna bases.

### merge_reads_command (raw_read_merger.py)
```merge_reads_command``` merges the subcommands to generate a full command of converting all the POD5 or FASTA5 in the reads directory to a merged pod5 file.

### validate_read_directory (raw_read_merger.py)
```validate_read_directory``` checks the existence of the reads directory.

### validate_target_directory (raw_read_merger.py)
```validate_target_directory``` checks the existence of the target directory.

### generate_merge_pod5 (raw_read_merger.py)
```generate_merge_pod5``` takes in the reads directory and target directory, runs the command generated by the ```merge_read_command``` and merges the pod5/fast5 files.

### consensus_generation (weighted_consensus_split.py)
```consensus_generation``` takes in working directory, reads, and a barcoded FASTA file and performs basecalling on them using Dorado. It then performs alignment using minimap2 set with custom parameters to maximize the amount of reads that align to the initial placeholder FASTA, filters with the primary aligned reads, and trims and sorts the read using default of 95% margin. The trimmed and sorted reads are then fed into vsearch clustering and the outcome consensus of each cluster is weighted by the number of reads that gets aligned to. Realign the reads to the weighted consensus and reperform the vsearch steps until the final consensus is generated.

### basecall_commands (consensus_methods_split.py)
```basecall_command``` generates a command to run the dorado basecaller.

### map_to_reference (consensus_methods_split.py)
```map_to_reference``` generates a command to use minimap2 to align the basecall with a reference sequence.

### filter_primary_alignments (consensus_methods_split.py)
```filter_primary_alignments``` filters a BAM file to retain only primary alignments and outputs the result to a new BAM file named 'primary_alignments.bam' in the same directory as the input file.

### read_trim (consensus_methods_split.py)
```read_trim``` takes in the SAM file generated from minimap2, trim the reads according to the length of the reference sequence and outputs the reads with approximately the same length in a FASTA file. This is done by truncating reads the bases that are longer than the 'dummy' reference file where we have a length approximation of the dataset. The output of read_trim is a FASTA file formatted as the following 

    >readID
    trimmed sequence

### write_to_fasta (consensus_methods_split.py)
```write_to_fasta``` takes in an output path and a file name, as well as a list of data to be written and writes that data out as a FASTA file at the inputted output path. 

### sort_fasta (consensus_methods_split.py)
```sort_fasta`` takes in a FASTA filepath and sorts it by length.

### vsearch_command (consensus_methods_split.py)
```vsearch_command``` takes in a FASTA path, a similarity threshold value, and output_path, and an output file name and generates a command to perform clustering and rough consensus on the FASTA file.

### mm2_cluster_aligner (consensus_methods_split.py)
```mm2_cluster_aligner``` takes in a cluster FASTA from VSEARCH as well as the original trimmed data set and performs realignment on the data using minimap2.

### weight_generation (consensus_methods_split.py)
```weight_generation``` takes in an aligned BAM (primary reads only) and returns a dictionary containing the references (clusters) from VSEARCH and the associated count of primary aligned reads.

### weighted_fasta_gen (consensus_methods_split.py)
```weighted_fasta_gen``` will take a cluster FASTA outputted from VSEARCH as well as a dictionary containing the reference sequences in that SAM file (clusters outputted from VSEARCH) and the number of reads aligned to that particular reference sequence. It will generate a FASTA file with the same reference sequence except multiplied by the number of reads that aligned to perfor weight cluster formation downstream. 

### rename_consensus_headers (consensus_methods_split.py)
```rename_consensus_headers``` will rename the headers in the consensus FASTA file based on the provided first and last N indexes. All headers will be renamed using the same j and k values.

### cluster_size_df (consensus_methods_split.py)
```cluster_size_df``` takes in the cluster FASTA generated by VSEARCH, extracts the size of each cluster during rounds of VSEARCH and stores the size of each cluster in every iteration in a dataframe for a later csv file generation.

### filter_cluster_size (consensus_methods_split.py)
```filter_cluster_size``` establishes a threshold and selectively filters the clusters based on the number of reads they contain during the rounds of VSEARCH clustering. Only clusters surpassing this specified threshold are considered eligible for the formation of a consensus sequence. The output of this function is a FASTA file containing the chosen clusters with the amount of reads grouped in the header. 

### merge_fasta_files (consensus_methods_split.py)
```merge_fasta_files``` merges multiple FASTA files into one.

### preprocessing (xf_model_training.py and xf_xna_find.py)
```preprocessing``` takes in a working directory, the raw data, the reference FASTA file and the direction (e.g. forward/reverse) and will generate necessary files (e.g. merged POD5 files and BAM files) needed for feature extraction or generation.

### raw_basecall_features (xf_model_training.py and xf_xna_find.py)
```raw_basecall_features``` runs the script to extract features from raw data & sequence space and merge them. Outputs each reference sequence as a JSON file containing the reads with their merged features. 

### consensus_features (xf_model_training.py and xf_xna_find.py)
```consensus_features``` takes in a the file path containing the JSON files with merged raw and sequence space data. Calculates consensus features for each of JSON file. Currently, data is stored in a list of pandas dataframes. 

### split_fasta (xf_xna_find.py)
```split_fasta``` takes in the input consensus FASTA files from the users and splits it by direction for the following preprocessing and feature extraction.

### dorado_bc_command (xf_basecall.py)
```dorado_bc_command``` generates a command to run the basecaller dorado.

### alignment_command (xf_basecall.py)
```alignment_command``` generates a command to run the minimap2 alignment.

### load_info_from_json (feature_extraction.py)
```load_info_from_json``` loads the information from a consensus JSON file that was exported from ```aggregate_reads.py```and provides them as a handful of variables.

### extract_coordinate_data (feature_extraction.py)
```extract_coordinate_data``` takes in the alignment data (a pandas series containing a list of three datasets per entry) and then extracts the data to individual pandas series.

### align_observation_ops (feature_extraction.py)
```align_observation_ops``` takes in a CIGAR and list of mapped signals, and then aligns the signal by CIGAR string operations.

### shift_to_alignment (feature_extraction.py)
```shift_to_alignment``` shifts the alignment of signals, base operations, and qualities to their proper corresponding alignment for a single read. 

### parse_read_alignment_info (feature_extraction.py)
```parse_read_alignment_info``` prepares read data, and then aligns it to a reference sequence and provides the resulting information.

### extract_alignment_from_dict (feature_extraction.py)
```extract_alignment_from_dict``` takes in a list of dictionaries of merged information of alignments, and the reference sequence for the reads in the list of dictionaries, and then generates the reference-based alignment data.

### calc_base_probs (feature_extraction.py)
```calc_base_prob``` takes in a dataframe of base sequence coordinates and calculates the numerical probability of each base at each position by the number of occurences.

### shannon_entropies (feature_extraction.py)
```shannon_entropies``` takes in a dataframe with base options as columns and the base positions as rows, with the values being probability of each base option at each position, and reuturns the shannon entropies of each position in the sequence.

### get_first_base (feature_extraction.py)
```get_first_base``` extracts the first base in a list of bases.

### check_base_match (feature_extraction.py)
```check_base_match``` checks if a base matches the reference base.

### mistmatch_prob (feature_extraction.py)
```mismatch_prob``` checks the aligned reads versus a reference sequence, and returns the probability of a mismatch at each base position.

### is_in_del (feature_extraction.py)
```is_in_del``` takes in a list of bases, and checks if the base is an insertion or deletion, and returns a list of booleans corresponding to if the base is either an insertion or deletion.

### filter_by_indel (feature_extraction.py)
```filter_by_indel``` takes in a list of booleans and values (of the same length) and returns the values where there was not an insertion or deletion.

### remove_in_dels (feature_extraction.py)
```remove_in_dels``` takes in base coordinate dataframe, and a base coordinate dataframe of data to be filtered and then filters that data by if it had insertions/deletions or not.

### getmean (feature_extraction.py)
```getmeans``` gets the mean value from a list of values, while checking if it's empty or not.

### mean_std_qual_by_base (feature_extraction.py)
```mean_std_qual_by_base``` gets the mean and standard deviation quality value by the base position.

### get_raw_median (feature_extraction.py)
```get_raw_median``` takes in a quality coordinate dataframe, and then provides the median quality of each base in the reference.
       
### convert_signals_to_objects (feature_extraction.py)
```convert_signals_to_objects``` converts each signal value in the sequence coordinates to a signal object using the convert_signal_list method.

### convert_signal_list (feature_extraction.py)
```convert_signal_list``` takes in a list of signals and a sampling frequency, and returns a list of signal objects.

### convert_signal_obs_to_groups (feature_extraction.py)
```convert_signal_obs_to_groups```maps the conversion of signal object lists to signal groups.

### convert_signal_obj_list (feature_extraction.py)
```convert_signal_obj_list``` converts a list of signal object list to signal groups.

### twas_the_monster_mash (feature_extraction.py)
```twas_the_monster_mash``` converts signal observations to usable features and squish it all together.

### feature_extraction (feature_extraction.py)
```feature_extraction``` takes in a path to a JSON file, and loads all of the features in that dataset
to a dataframe, with a column called 'XNA_PRESENT' representing the index where the xna is present.

### batch_read_alignment_data (feature_extraction.py)
```batch_read_alignment_data``` takes in a list of read alignment data, and a batch size, and generates a list of subdivided batches of the size of the batch. if there arent enough to fill the last, it populates it with what it can - occasionally resulting in smaller batches.

### extract_batch_features (feature_extraction.py)
```extract_batch_features``` takes in a read batch and the reference data, and then generates the features for that batch as if it were an entire consensus.

### batched_feature_extraction (feature_extraction.py)
```batched_feature_extraction``` takes in a JSON filepath and a batch size, and returns a list of the assembled features as batches.

### class Signal (feature_extraction.py)
class ```Signal``` contains and automatically runs a handful of valuable pieces of information about a signal list, such as mean, std, etc.

### class GroupedSignals (feature_extraction.py)
class ```GroupedSignals``` contains useful methods for processing groups of signal objects.

### export_to_json (aggregate_reads.py)
```export_to_json``` takes in a directory path and saves the ConsensusReadsData object to a JSON file at that path.

### feature_extraction (aggregate_reads.py)
```feature_extraction``` runs the feature extraction methods on the data stored within the current object.

### save_batches_to_parquet (aggregate_reads.py)
```save_batches_to_parquet``` takes in an output superdirectory, and saves the batched consensus features as parquets to that directory.

### percenttobar (aggregate_reads.py)
```percenttobar``` converts a fractional value to a loading progress bar string.

### load_bam (aggregate_reads.py)
```load_bam``` takes in a path to a basecalled and aligned consensus BAM file, and converts the data into a list of dicts format.

### load_fasta (aggregate_reads.py)
```load_fasta``` takes in a consensus fasta path, and returns a list of dicts containing the read IDs and sequences of all valid samples.

### map_signal_to moves (aggregate_reads.py)
```map_signal_to_moves``` takes in a moves table, base offeset, and a signal, and then maps the signal to each corresponding move/observation in the table.

### iterable_reads (aggregate_reads.py)
```iterable_reads``` uses a path to a POD5 file and a list of valid read IDs, and returns an iterable generator of Read objects constructed from the POD5 file.

### aggregate_reads_to_reference (aggregate_reads.py)
```aggregate_reads_to_reference``` generates a list of JointBamReads based on a POD5 path and BAM data.

### iterable_merged_data (aggregate_reads.py)
```create an generator``` of ConsensusReadsData objects from alignment data and the fasta data, as well as a list of JointBamReads objects.

### class Read (aggregate_reads.py) 
class ```Read``` stores the read ID, signal, and frequency of a given read.

### class JointBamReads (aggregate_reads.py)
class ```JointBamReads ``` stores the reference name and read iterator for a given reference sequence.

### class ConsensusReadsData (aggregate_reads.py)
class ```ConsensusReadsData``` stores the reference name, sequence,frequency, and list of reads data.

### fetch_xna_pos (xr_fasta2x_rc.py)
```fetch_xna_pos``` extracts the position of XNA from the header in the input FASTA file.

### xna_base_rc (xr_fasta2x_rc.py)
```xna_base_rc``` takes in an XNA base and outputs its reverse complement base.

### reverse_complement_fasta (xf_rc_fasta_gen.py)
```reverse_complement_fasta``` reads the input FASTA file and creates the reverse complement.