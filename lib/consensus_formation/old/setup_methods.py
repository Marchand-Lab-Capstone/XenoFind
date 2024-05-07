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
    
    # Use check_make_dir to generate or validate the existing directories.
    wdir = check_make_dir(working_directory)
    con_file_dir = check_make_dir(working_directory + "consensus_files")
    p5dir = check_make_dir(con_file_dir + "merged_pod5")
    bcdir = check_make_dir(con_file_dir + "basecall_directory")
    forward_dir = check_make_dir(con_file_dir + "forward_reads")
    fadir_fwd = check_make_dir(forward_dir + "fasta_directory")
    reverse_dir = check_make_dir(con_file_dir + "reverse_reads")
    fadir_rev = check_make_dir(reverse_dir + "fasta_directory")
    reverse_unflip_dir = check_make_dir(con_file_dir + 'reverse_unflipped_reads')
    fadir_unrev = check_make_dir(reverse_unflip_dir + 'fasta_directory')
    total_dir = check_make_dir(con_file_dir + "total_reads")
    fadir_tot = check_make_dir(total_dir + "fasta_directory")
    vsearch_dir_fwd = check_make_dir(forward_dir + "vsearch_processing")
    vsearch_dir_rev = check_make_dir(reverse_dir + "vsearch_processing")
    vsearch_dir_rev_unflip = check_make_dir(reverse_unflip_dir + 'vsearch_processing')
    vsearch_dir_tot = check_make_dir(total_dir + "vsearch_processing")
    xfconsdir = check_make_dir(con_file_dir + "xf_consensus_output")
    xfconsdir_fwd = check_make_dir(xfconsdir + "forward_reads")
    xfconsdir_rev = check_make_dir(xfconsdir + "reverse_reads")
    xfconsdir_rev_unflip = check_make_dir(xfconsdir + 'reverse_unflipped_reads')
    xfconsdir_all = check_make_dir(xfconsdir + "total_reads")
    return [wdir, con_file_dir, p5dir, bcdir, forward_dir, fadir_fwd, reverse_dir, fadir_rev, reverse_unflip_dir, fadir_unrev, total_dir, fadir_tot, vsearch_dir_fwd, vsearch_dir_rev, vsearch_dir_rev_unflip, vsearch_dir_tot, xfconsdir, xfconsdir_fwd, xfconsdir_rev, xfconsdir_rev_unflip, xfconsdir_all]
