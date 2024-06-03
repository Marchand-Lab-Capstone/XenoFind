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
        print('XenoFind [STATUS] - Required directory not found. Creating directory: '+directory)
    return directory


def setup_directory_system_con(working_directory):
    """
    setup_directory_system takes in a working directory filepath, and
    ensures that the required subdirectories are extant.
    
    Parameters: 
    working_directory: a str representing a directory to do all the processes in
    
    Returns: a 9-value list with the paths for the following:
            working_directory 0,
            consensus_files_directory 1,
            merged_pod5 2,
            basecall_directory 3,
            preprocessing_directory 4,
            preprocessing_fasta_directory 5,
            preprocessing_vsearch_directory 6,
    """
    
    # Use check_make_dir to generate or validate the existing directories.
    wdir = check_make_dir(working_directory)
    con_file_dir = check_make_dir(os.path.join(working_directory,"consensus_files/"))
    p5dir = check_make_dir(os.path.join(con_file_dir,"merged_pod5/"))
    bcdir = check_make_dir(os.path.join(con_file_dir,"basecall_directory/"))
    preprocess_dir = check_make_dir(os.path.join(con_file_dir,"preprocessing/"))
    fadir_pre = check_make_dir(os.path.join(preprocess_dir,"fasta_directory/"))
    vsearch_dir_pre = check_make_dir(os.path.join(preprocess_dir,"vsearch_processing/"))

    return [wdir, con_file_dir, p5dir, bcdir, preprocess_dir, fadir_pre, vsearch_dir_pre]

def setup_directory_system_model(working_directory):
    """
    setup_directory_system takes in a working directory filepath, and
    ensures that the required subdirectories are extant.
    
    Parameters: 
    working_directory: a str representing a directory to do all the processes in
    
    Returns: a nth-value list with the paths for the following:

    """
    
    # Use Check_make_dir to generate or validate the existing directories.
    wdir = check_make_dir(working_directory) #0
    model_train_dir = check_make_dir(os.path.join(working_directory, 'model_training')) #1
    ref_dir = check_make_dir(os.path.join(working_directory, 'reference')) #2
    p5dir = check_make_dir(os.path.join(working_directory, 'merged_pod5')) #3
    bcdir = check_make_dir(os.path.join(model_train_dir, 'basecall_directory')) #4
    json_dir = check_make_dir(os.path.join(model_train_dir, 'json_files')) #5
    parquet_dir = check_make_dir(os.path.join(model_train_dir, 'parquet_superdirectory')) #6
    model_dir = check_make_dir(os.path.join(model_train_dir, 'model')) #7
    return [wdir,model_train_dir, ref_dir, p5dir,bcdir, json_dir, parquet_dir, model_dir]

def setup_directory_system_find(working_directory):
    """
    setup_directory_system takes in a working directory filepath, and
    ensures that the required subdirectories are extant.
    
    Parameters: 
    working_directory: a str representing a directory to do all the processes in
    
    Returns: a nth-value list with the paths for the following:

    """
    # Use Check_make_dir to generate or validate the existing directories.
    wdir = check_make_dir(working_directory) #0
    find_dir = check_make_dir(os.path.join(working_directory, 'xna_find')) #1
    ref_dir = check_make_dir(os.path.join(working_directory, 'reference')) #2
    p5dir = check_make_dir(os.path.join(working_directory, 'merged_pod5')) #3
    bcdir = check_make_dir(os.path.join(find_dir, 'basecall_directory')) #4
    json_dir = check_make_dir(os.path.join(find_dir, 'json_files')) #5
    parquet_dir = check_make_dir(os.path.join(find_dir, 'parquet_superdirectory')) #6
    return [wdir,find_dir, ref_dir, p5dir,bcdir, json_dir, parquet_dir]
    
def contains_xna_bases(fasta_file, xna_base_pairs):
    with open(fasta_file, 'r') as f:
        for line in f:
            if not line.startswith('>'):
                for base_pair in xna_base_pairs:
                    if any(base in line for base in base_pair):
                        return True
    return False

