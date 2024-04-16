"""
setup_methods.py
J. Sumabat, N. Lai, S. Peck, Y. Huang, 3/12/2024

setup_methods.py is a file that contains relevant tools and methods to set up consensus.py's functionality, but can be used for other purposes. 

check_make_dir() - check if a directory string exists, and if it doesn't, make it.
setup_directory_system() - setup the directory system used for consensus_methods.py
"""

import os


def check_make_dir(directory):
    """
    ADAPTED FROM XF_TOOLS.PY, By: By: J. Sumabat, N. Lai, Y. Huang, S. Peck, H. Kawabe, N. Kaplan, J. A. Marchand
    
    
    check_make_dir is a function that takes in a directory in the form of a string to actually generate that directory 
    by using the built in os.makedirs() function. The string that is generated is either the inputted working directory 
    string from the user or other predetermined directories that are put into the working directory. There are added print statements to 
    help users know what directory is being generated. 
    
    Parameters:
    directory: a directory path, as a string
    returns: 
    """
    directory = os.path.expanduser(directory)
    if not os.path.exists(directory):
        os.makedirs(directory)
    return str(os.path.abspath(directory))+"/"


def setup_directory_system(working_directory):
    """
    setup_directory_system takes in a working directory filepath, and
    ensures that the required subdirectories are extant.
    
    Parameters: 
    working_directory: a str representing a directory to do all the processes in
    
    Returns: a 11-value list with the paths for the following:
             basecall_directory,
             forward_reads_directory,
             forward_reads_fasta_directory,
             reverse_reads_directory,
             reverse_reads_fasta_directory,
             all_reads_directory,
             all_reads_fasta_directory,
             merged_pod5,
             rough_consensus_output,
             xf_consensus_output
    """
    
    # Use Check_make_dir to generate or validate the existing directories.
    wdir = check_make_dir(working_directory)
    xna_find_dir = check_make_dir(working_directory + 'XNA_find/')
    p5dir = check_make_dir(xna_find_dir + "merged_pod5")
    bcdir = check_make_dir(xna_find_dir + "basecall_directory/")
    return [wdir,xna_find_dir,p5dir,bcdir]