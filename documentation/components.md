# Components

## Big Idea/Script Components 

### Command Line 
This would be the main UI for users. A user would call ```python xenofind.py [method] [fast5/pod5 directory] [reference fasta]``` with required inputs. Progress and errors would show up here.

### xenofind_pipe.py 
This will be the other 'UI' for users. Users would input where they want the working directory to be, as well as the directories for inputted fast5/pod5 and the reference fasta for basecalling. In addition, you would set which detection method you would want to use on this script. This is essenetially pretty much the same as running the command line as mentioned above except it those inputs would be stored and easy to see. 

### xenofind.py
This is the primary script. This will call the desired detection method in the ```lib``` directory regard if it is called using the command line or xenofind_pipe.py. 

### xf_params
This will contain the parameters used in any of the detection methods such as the path for the basecaller. It also contains toggleable variables to only run certain parts of the script (e.g., initial basecall). This is imported into most of the scripts in this program. 

### xf_tools 
This script contains all the reusable, general functions that can be called as needed (e.g., fast5 to pod5 conversion). This script is imported to most of the scripts in this program. 

### consensus.py
This script takes in the fast5 files generated from Nanopore Sequencing and outputs a fasta file as the reference sequence for downstream possible location detection of the XNAs.

### xf_low_qual.py
This script outlines the first approach for identifying potential positions of XNA bases with limited accuracy and precision. In this method, the raw datasets will be merged and basecalled with the fasta reference sequence from ```consensus.py``` into a ```.bam``` file, and the features of the reads in the bam file will be generated. Then it conducts a statistical test to identify low-quality scores on a per-read basis. Following this assessment, the reads are recentered and grouped to the reference sequence. Finally, positions with low-quality scores indicative of XNA presence are extracted.

## Internal Components/ Functions 

### read_trim(consensus.py)
This function takes in the sam file generated from minimap2, trim the reads according to the length of the reference sequence and outputs the reads with approximately the same length in a fasta file.

### random_fasta_sample(consensus.py)
This function takes in the fasta file with the reads sorted by length, generates a list containing the paired headers and sequences and select X sequences randomly without replacement.

### cluster_size_filter(consensus.py)
Following the VSEARCH clustering, this function establishes a threshold and selectively filters the clusters based on the number of reads they contain. Only clusters surpassing this specified threshold are considered eligible for the formation of a consensus sequence.

### extract_read_info function (xf_low_qual.py)
This function takes in the merged bam file and extracts the readID, basecalled sequence, start of reference sequence, the quality socre per base and the average quality score of each read.

### xna_guess function (xf_low_qual.py)
This function performs a statistical test on a per base read to indentify low-quality scores, recenters and groups the reads to the reference sequence, and extracts the positions indicative of the XNA presence.
