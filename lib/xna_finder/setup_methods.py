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
    
    Returns: a nth-value list with the paths for the following:

    """
    
    # Use Check_make_dir to generate or validate the existing directories.
    wdir = check_make_dir(working_directory) #0
    xna_find_dir = check_make_dir(working_directory + 'XNA_find/') #1
    p5dir = check_make_dir(xna_find_dir + "merged_pod5") #2
    bcdir = check_make_dir(xna_find_dir + "basecall_directory/") #3
    fwd_dir = check_make_dir(xna_find_dir + 'forward_read') #4
    fwd_qs = check_make_dir(fwd_dir + 'quality_score') #5
    fwd_shannon = check_make_dir(fwd_dir + 'shannon_entropy') #6
    fwd_signal = check_make_dir(fwd_dir + 'signal_level') #7
    rev_dir = check_make_dir(xna_find_dir + 'reverse_read') #8
    rev_qs = check_make_dir(rev_dir + 'quality_score') #9
    rev_shannon = check_make_dir(rev_dir + 'shannon_entropy') #10
    rev_signal = check_make_dir(rev_dir + 'signal_level') #11
    rev_unflip_dir = check_make_dir(xna_find_dir + 'reverse_unflipped_read') #12
    rev_unflip_qs = check_make_dir(rev_unflip_dir + 'quality_score') #13
    rev_unflip_shannon = check_make_dir(rev_unflip_dir + 'shannon_entropy') #14
    rev_unflip_signal = check_make_dir(rev_unflip_dir + 'signal_level') #15
    total_dir = check_make_dir(xna_find_dir + 'total_read') #16
    total_qs = check_make_dir(total_dir + 'quality_score') #17
    total_shannon = check_make_dir(total_dir + 'shannon_entropy') #18
    total_signal = check_make_dir(total_dir + 'signal_level') #19
    return [wdir,xna_find_dir,p5dir,bcdir, fwd_dir, fwd_qs, fwd_shannon, fwd_signal, rev_dir, rev_qs, rev_shannon, rev_signal, rev_unflip_dir, rev_unflip_qs, rev_unflip_shannon, rev_unflip_signal, total_dir, total_qs, total_shannon, total_signal]
