# Let's take this shit step by step
'''
get the working directory
get the raw directory
get the reference fasta file

generate the working directory (how?)
generate processing directory
generate pod5 directory
generate basecalling directory
generate bam file storage
generate fastq storage

take in the filetype, if fast5, convert to combined pod5

if basecalling is needed, basecall the pod5 file using a [unknown basecaller] and convert to bam file

analyze the fastq file and extract the read information
    - take pysam and read the bam file
    - if the read is mapped,
        -get the qualities of each query
        -get an array of the qualities
        -get the mean of the qualities
        -generate the features: query name, sequence, start position relative to reference, quality scores of read, and the average quality
    - if the read is not mapped, let user know there is no mapping in the bam file and it does not have a reference position
    
    - return the list of read info, which is a list of the lists of features on a per-read basis
    

check to see if each base pair has its own quality score

        
'''

import sys
import os
import numpy as np
import consensus_refactor as crf

EXTRA_DIRECTORIES = ['pod5','basecall','bam','fastq']
MERGE_FILENAME = 'merged_reads'

def get_working_directory(path=None):
    response_val=None
    
    # if a path was passed
    if (path not None):
        
        # Check the filestructure
        file_checkup = check_filestructure(path)
        
        # if the directory passed is valid
        if (file_checkup['is_valid']):
            response_val = path
        else:
            os.make_dirs(path)
            
        # if an invalid directory, exit
        else:
            response_val = sys.exit()
            
    # if nothing pased, prompt
    else:
        path = input("Provide working directory (ends with a /): ")
        file_checkup = check_filestructure(path)
        # if the directory passed is valid
        if (file_checkup['is_valid']):
            
            # check if each subfolder is in there, and if it isn't make it.
            for folder in EXTRA_DIRECTORIES:
                if (folder not in file_checkup['subfiles']):
                    st = os.makedirs(path+folder)
            # Set the response value to the passed filepath
            response_val = path
            
        # if an invalid directory, exit
        else:
            response_val = sys.exit()
        
    return response_val
        
        
def check_filestructure(path, t='d', ext='fa'):
    is_valid = None
    is_empty = None
    subfiles = []

    if t = 'd':
        is_valid = os.path.isdir(path)
        subfiles = os.path.listdir(path)
        is_empty = (len(subfiles) == 0)
        
    elif t = 'f':
        is_valid = (path.split('.')[-1] == ext)
        
    else:
        is_valid = None
    
    out_dict = {'path':path,
                'is_valid':is_valid,
                'is_empty':is_empty,
                'subfiles':subfiles}
    
    return out_dict
        
        
    
    

def request_directories():
    '''
    request_directories is designed to prompt the user for the working directory,
    the directory with the raw read data, and the directory with the reference fasta file.
    If one or all of these do not exist, the program using this method will exit.
    Also contains segments to check if the files passed are valid for their purposes and 
    prompts for overwriting if nessecary. 
    
    Parameters:
    None
    
    Returns:
    a touple containing the working directory, read directory, and the reference file paths.
    - or -
    sys.exit() - only if invalid input or user choice
    
    TODO: Modify this so that it does proper error handling, as well as better handling of directory confirmation
          Should pass sys.exit as return instead of calling it within method
          Should have default paths as none that only requests new ones if the defaults are not handled
    
    '''
    
    #--Working Directory Segment----------------------------------------------#
    # prompt for working directory
    w_dir = input("Provide working directory (ends with a /): ")
    
    # Check if directory exists
    w_is_valid = os.path.isdir(w_dir)
    
    # Create list for all subdirectories to be generated
    fnames = EXTRA_DIRECTORIES
    
    #if valid directory:
    if w_is_valid:
        # list of directories in passed dir:
        list_working_dirs = list(os.listdir(w_dir))
        
        # Check if directory is empty
        is_empty = len(list_working_dirs)==0
        
        # if it is empty
        if is_empty:  
            # generate the working subdirectories with the folder filler method
            ex_stat = folder_filler(w_dir, fnames)
        
        # if it isnt, check if the subfolders are there
        elif (set(EXTRA_DIRECTORIES).issubset(list_working_dirs)):
            print("All directories present.")
        
        # otherwise, prompt to see if the user wants to add folders.
        else:
            overwrite = input("WARNING: {} contains files. The following folders will be added: {}. Continue? [y/n]: ".format(w_dir, fnames)) == 'y'
            # if the user wants to overwrite, 
            if (overwrite):
                #print("Overwriting...")
                ex_stat = folder_filler(w_dir, fnames)
                
            # otherwise, exit.
            else:
                print("Exiting.")
                sys.exit()
                
    # if the directory is invalid:
    else:
        # prompt the user to see if they want to make the directory.
        build = input("WARNING: {} does not exist. Should it be created? [y/n]: ".format(w_dir)) == 'y'
        
        # if the user agrees, build the directory.
        if(build):
            #print("Creating directory at {}...".format(w_dir))
            st = os.makedirs(w_dir)
            folder_filler(w_dir, fnames)
            
        # if the user doesn't, exit.
        else:
            print("Exiting.")
            sys.exit()
    
    #--Reads directory segment------------------------------------------------#
    # Prompt user for read directory
    read_dir = input("Provide directory that contains raw read data: ")
    
    # check if the read directory is valid:
    readdir_is_valid = os.path.isdir(read_dir)

    # if the read directory is valid
    if readdir_is_valid:
        
        # Check if the directory is empty.
        is_empty = (len(os.listdir(read_dir)) == 0)
        
        # if it is not valid, exit. 
        if is_empty:
            print("ERROR: {} is empty. Exiting.")
            sys.exit()
            
    #--Reference fasta segment------------------------------------------------#
    # prompt user for the reference fasta path
    ref_dir = input("Provide path to reference fasta file: ")
    
    # get the file extension
    extension = ref_dir.split('.')[-1]
    
    # if the the extension is not a fasta, exit.
    if ('fa' not in extension):
        print("Error: {} is not a valid fasta file. Exiting.".format(ref_dir))
        sys.exit()
    
    # Return the working directory, read directory, and reference directory paths
    return w_dir, read_dir, ref_dir


def folder_filler(core_dir, folder_names):
    '''
    folder_filler takes in a core working directory and a list of folder names,
    and generates the subfolders in the core directory.
    
    Parameters:
    core_dir: working directory path as a string
    folder_names: a list containing strings of subfolders to be made
    
    Returns:
    A boolean, if True, nothing went wrong.
    '''
    # Create a list of statuses
    statuses = []
    
    # loop through each folder in the folder names list
    for folder in folder_names:
        
        # Append the status of the os.makedirs command for the folder.
        statuses.append((os.makedirs(core_dir+folder)))
        
    # Create a boolean that checks if all the folders were able to be made. DOENT WORK
    #exit_stat = not((np.asarray(statuses).mean()) < 1)
    
    # Return the exit status.
    return statuses


# These are not currently working due to pysam being a lil nerd
"""
def extract_read_info(bam_file_path):
    '''
    This function takes in the bam file generated and extracts the readID, basecalled sequence, start of reference sequence, and the quality score from the bamfile 
    
    Parameters:
    bam_file_path: a path to the .bam file
    
    Returns:
    read_info: a list containing lists of features (query name, sequence, start position relative to reference, quality scores of read, and the average quality)
    '''
    read_info = []
    # Open the BAM file
    with pysam.AlignmentFile(bam_file_path, "rb") as bamfile:
        for read in bamfile:
            # Check if the read is mapped
            if not read.is_unmapped:
                qual = read.query_qualities 
                qual = np.array(qual, dtype=int)
                avg_qual = np.mean(qual)
                features = [
                read.query_name,  # Query name of the read
                read.query_sequence,  # Sequence basecalled
                read.reference_start,  # Position of the read relative to the reference
                qual,  # Quality scores of the read (numerical)
                avg_qual
                ]
                read_info.append(features)
            else:
                print(f"Read {read.query_name} is unmapped and does not have a reference position.")
    return read_info

    
def confirm_bases(basecalled, quality_scored):
    '''
    confirm_bases takes in a list of bases in a sequence and a list of quality scores, and checks if they match
    Also, prints if they match.
    Parameters:
    basecalled: basecall sequence
    quality_scored: corresponding quality scores to a basecall
    
    Returns:
    True if lengths match, False if not.
    '''
    matching_len = len(basecalled) == len(quality_scored)
    if matching_len:
        print('Yes, the length of the sequences matches the length of quality score string')

    else:
        print('No, the length of the sequences and quality score doesnt match!')
    
    return(matching_len)
""" 


# There are three core files needed:
#1: The directory where the operations and output will be performed and held
WORKING_DIR = None

#2: The directory that contains the raw read data
RAW_DIR = None

#3: The directory that contains the reference fasta file, which can either be consensus or otherwise
REFERENCE_DIR = None

# Get the directories using request_directories()
WORKING_DIR, RAW_DIR, REFERENCE_DIR = request_directories()

# use consensus-refactor.py's get_read_data method to get the read data in the directories
merged_files = crf.get_read_data(RAW_DIR, WORKING_DIR+EXTRA_DIRECTORIES[0], MERGE_FILENAME)

# Ignoring basecalling for now, assume xf has a path to the bam file
bam = WORKING_DIR+EXTRA_DIRECTORIES[2]+'/bc.bam'

'''
# Extract the read information from the bam file
bam_read = extract_read_info(bam)

# Check if the base pairs have the proper quality scores using only the first read in the bam.
# Note: This really should be checking all of them, not just the first.
bases_confirmed = confirm_bases(bam_read[0][1], bam_read[0][3])
'''
# Guess which positions have XNA
