"""
consensus-refactor.py - Steps to make an AGTC referece file for later comparison to xenobases.
"""

# Import Statements:
import os
import random
import subprocess
from Bio import SeqIO
import xf_tools as xt
import xf_params as xp
import platform
SYS = platform.system()


##--filetype handler methods--##
def fast5_handler(filepath):
    # Get the directory of the fast5 file
    file_dir = filepath.split(".fast5")[0]

    # Create the output filepath for the pod5 file
    output_filepath = file_dir + ".pod5"
    
    # Create empty status variable:
    status = None
    
    # If the pod5 conversion doesn't exist
    if (not os.path.isfile(output_filepath)):
        # Convert the fast5 at the filepath to a pod5 file
        cmd = "pod5 convert fast5 {} -o {}".format(filepath, output_filepath)
        st = os.system(cmd)
        
        # Create a status string to check in case of verbosity
        status = ("Error! failed in "*st) +"Converted {} to {}".format(filepath, output_filepath)
        print(status)
        # Return the output filepath
        return output_filepath
    
    # If pod5 conversion exists, 
    else:
        # Return the pod5 filepath
        status = "{} already exists.".format(output_filepath)
        print(status)
        return output_filepath



def pod5_handler(filepath):
    print("{} found.".format(filepath))
    return filepath

##--End Filetype handler methods--##


## STEP 1: TAKE FAST5 OR POD5 FILES AND CONDENSE THEM INTO A SINGLE POD5 FILE ##

def merge_command(in_dir, out_filename):
    # method to handle a merge command for pod5 with a specific input directory and output filename
    cmd = ""
    
    # Check OS to do conditional formatting
    if (SYS == 'Windows'):
        cmd = "pod5 merge --recursive --force-overwrite {} -o {}".format(in_dir, out_filename)
    else:
        cmd = "pod5 merge --force-overwrite {}*.pod5 -o {}".format(in_dir, out_filename)
        
    # Run the command, return the status
    st = os.system(cmd)
    return st


def get_read_data(in_dir, out_dir, filename):
    # Get all the read files in the specified directory,
    # and convert them ot a merged pod5 file in the output dir
    
    # List all the files in the directory
    list_files = os.listdir(in_dir)
    
    # Create an empty list to store all the files to be merged
    files_to_merge = []
    
    # Check if the output path exists:
    if (not os.path.exists(out_dir)):
        # if not, make it exist
        os.mkdir(out_dir)
        print("Folder created at {}".format(out_dir))
    
    # Loop through each file in the list
    for file in list_files:
        
        # get the filepath
        filepath = in_dir + file
        
        # Get the extension of the file
        extension = file.split(".")[1]
        
        # handle each extension type
        switch = {"fast5": fast5_handler,
                  "pod5": pod5_handler
                  #Extension name: extension conversion handler,
                 }
        
        # Get the path of the corresponding pod5 file based on the case
        pod5_version_path = ""
        try:
            pod5_version_path = switch[extension](filepath)
        except:
            pod5_version_path = "err"
        
        # Merge the 
        files_to_merge.append(pod5_version_path)
    
    # remove duplicates from the merge list:
    files_to_merge = list(set(files_to_merge))
    
    # Get the output filename:
    output_filename = "{}{}.pod5".format(out_dir, filename)
    
    # if the output file does not yet exist:
    if (not os.path.exists(output_filename)):
        
        # merge the files using pod5:
        st = merge_command(in_dir, output_filename)
        print("Un"*st + "successful merge into {}.".format(output_filename))
    
    # If the output file does eixist:
    else:

        # Set up a boolean to check for valid user input
        validuserinput = False
        
        # While there isnt valid input,
        while (validuserinput == False):
            # Prompt user for if the old merge file should be overwritten
            userin = input("File already exists at {}. Overwrite? [y/else]".format(output_filename))
            
            # if y, merge. if not y, ignore the merge.
            if userin == 'y':
                validuserinput = True
                st = merge_command(in_dir, output_filename)
                print("Un"*st + "successful merge into {}.".format(output_filename))
            else:
                validuserinput = True
                print("Ignoring merge.")
    
    return files_to_merge
    
## END STEP 1 ##
# Example use: get_read_data("../Data/sample_fast5/subset_fast5/", "../Data/sample_fast5/merged_folder/", "merged")

## STEP 2: BASECALL THE MERGED POD5 FILE ##

def basecall(merged_path):
    return None