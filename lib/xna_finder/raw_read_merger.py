"""
raw_read_merger.py
J. Sumabat, N. Lai, S. Peck, Y. Huang, 3/12/2024

raw_read_merger.py contains methods used for taking a directory of either fast5 or pod5 files,
and then merging them into a singular pod5 file. 

merge_reads_command() - generate the merge command dependant on the OS and filetype.
validate_read_directory() - check that the passed read directory actually exists.
validate_target_directory() - check if a passed directory exists. 
generate_merged_pod5() - use a directory of reads to generate the merged pod5, agnostic of filetype
"""


import os
import platform
SYS = platform.system()

def merge_reads_command(reads_dir, filetype, target_dir_path):
    # This already assumes the read directory and target directory are valid
    cmd = ""
    subcommand = ""
    os_command = ""
    output_filename = target_dir_path
    
    if filetype == 'fast5':
        subcommand = 'convert fast5'
    elif filetype == 'pod5':
        subcommand = 'merge'
    
    os_command = '{}/*.{}'.format(reads_dir, filetype) #need to account for when users enter their file path as A/B/C/ <-- last slash included
    cmd = "pod5 {} --force-overwrite {} -o {}".format(subcommand, os_command, output_filename)
    print(cmd)

    return cmd


def validate_read_directory(reads_dir):
    directory_exists = os.path.isdir(reads_dir)
    homogenous_files = True
    filetype = ""
    if directory_exists:
        directory_files = os.listdir(reads_dir)
        ext_list = []
        for file in directory_files:
            ext = file.split('.')[-1]
            ext_list.append(ext)
        uniques = list(set(ext_list))
        if (len(uniques) != 1):
            homogenous_files = False
            print('Passed reads directory not homogenous. Filetypes found: {}'.format(uniques))
        else:
            filetype = uniques[0]
    return filetype


def validate_target_directory(target_dir):
    directory_exists = os.path.isdir(target_dir)
    return directory_exists
        
    
def generate_merged_pod5(reads_dir, target_dir_path):
    print('XenoFind [STATUS] - Generating merged pod5 command')
    filetype = validate_read_directory(reads_dir)
    if filetype:
        cmd = merge_reads_command(reads_dir, filetype, target_dir_path)
        print('XenoFind [STATUS] - running merge command...')
        st = os.system(cmd)
    return filetype, target_dir_path
    
