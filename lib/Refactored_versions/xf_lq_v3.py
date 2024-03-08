import sys
import os
import numpy as np
import xf_tools as xft
import xf_fastpod_handler as xpod

WORKING_DIR = None
EXTRA_DIRECTORIES = ['pod5','basecall','bam','fastq']
RAW_DIR = None
MERGED_FILENAME = 'merged' #.pod5
REF_DIR = None
BASECALL_DIR = None

def get_directory(path, typ):
    # valid types: 'working', 'raw read', 'reference', 'basecall' 
    directory = None
    
    if (type(path) != type(None)):
        directory = xft.check_make_dir(path) # THIS ASSUMES THE USER PASSES A VALID DIRECTORY STR
        print(validate_dirtype(directory, typ)) # if false, handle based on type.
        
    else:
        path = input("Provide {} directory (ends with /): ".format(typ))
        directory = xft.check_make_dir(path)
        print(validate_dirtype(directory, typ))
        
    return directory

def generate_working_folders(path):
    for directory in EXTRA_DIRECTORIES:
        true_dir = os.path.join(path, directory)
        xft.check_make_dir(true_dir)
    return 'working folders generated'
        
def validate_dirtype(path, typ):
    output = None
    if typ == 'working': output = generate_working_folders(path)
    elif typ == 'raw read': output = (xpod.get_read_data(path, WORKING_DIR+EXTRA_DIRECTORIES[0], "/" + MERGED_FILENAME))
    elif typ == 'reference': output = ((path.split('.')[-1] == 'fa')) # if this returns false, generate a fasta
    elif typ == 'basecall': output = ((path.split('.')[-1] == 'bam')) # if this returns false, basecall the fasta
    else: output = ('unrecognized typecall "{}"'.format(typ)) # handle error
    return output


def extract_read_info(bam_path):
    # I ASSUME THIS WORKS, I CANT TEST IT ON WINDWS CUZ PYSAM DONT WORK HERE
    reads_dict = {}
    # open the alignmentfile.
    with pysam.AlignmentFile(bam_path, "rb") as bamfile:
        
        # for each read in the bam file
        for read in bamfile:
            # Check if the read is mapped:
            if not read.is_unmapped:
                
                # Get name, reference start position, and query sequence
                name = read.query_name
                refpos = read.reference_start
                sequence = read.query_sequence
                
                # turn the read qualities into ints in a numpy array
                quals = np.array(read.query_qualities, dtype=int)
                
                # get the mean quality
                mean_qual = np.mean(qual)
                reads_dict[name] = {'reference_start':refpos,
                                    'sequence':sequence,
                                    'quals':quals,
                                    'mean_qual':mean_qual}
            else:
                print(f"Read {read.query_name} is unmapped and does not have a reference position.")
        # if we do validation that the sequence lengths match, we do it before returning the read.
    return reads_dict


def run_script():        
    # SETUP 
    WORKING_DIR = get_directory(WORKING_DIR, 'working') # this needs to work for everything else to work
    RAW_DIR = get_directory(RAW_DIR, 'raw read')
    REF_DIR = get_directory(REF_DIR, 'reference')
    BASECALL_DIR = get_directory(BASECALL_DIR, 'basecall')

    # extract read information
    basecall_dict = extract_read_info(BASECALL_DIR)
    return basecall_dict

