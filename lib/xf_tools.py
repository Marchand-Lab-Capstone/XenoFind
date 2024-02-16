#####################################################################################
''' 
tools.py

Useful functions to be called on  by other scripts.

By: By: J. Sumabat, N. Lai, Y. Huang, S. Peck, H. Kawabe, N. Kaplan, J. A. Marchand

'''
#####################################################################################
import os
from pathlib import Path

#Scans directory and subdirectory to get proper fast5 file path. Not explicitly handled with pod5 commands
def get_fast5_subdir(fast5_dir): 
    path = os.path.normpath(os.path.expanduser(fast5_dir))
    if os.path.exists(path):
        fast5_files = list(Path(path).rglob("*.fast5" ))
        if len(fast5_files)>0:
            fast5_subdir = os.path.dirname(fast5_files[0])
            print('Xemora [STATUS] - Found '+str(len(fast5_files))+' fast5 files in '+fast5_dir)
            return fast5_subdir
        else: 
            print('Xemora [ERROR] - Could not find Fast5 files in specified directory. Check .fast5 exist.')
            return False
    else: 
        print('Xemora [ERROR] - Could not find Fast5 directory. Check path')
        return False

#Scans directory and subdirectory to get proper pod5 file path. 
def get_pod5_subdir(pod5_dir): 
    path = os.path.normpath(os.path.expanduser(pod5_dir))
    if os.path.exists(path):
        pod5_files = list(Path(path).rglob("*.pod5" ))
        if len(pod5_files)>0:
            pod5_subdir = os.path.dirname(pod5_files[0])
            print('Xemora [STATUS] - Found '+str(len(pod5_files))+' POD5 files in '+pod5_dir)
            return pod5_subdir
        else: 
            print('Xemora [ERROR] - Could not find POD5 files in specified directory. Check .pod5 exist.')
            return False
    else: 
        print('Xemora [ERROR] - Could not find POD5 directory. Check path')
        return False

#Check if working directory exists, if not create it. 
def check_make_dir(directory):
    directory = os.path.expanduser(directory)
    if not os.path.isdir(directory):
        os.makedirs(directory)
        print('Xemora [STATUS] - Required directory not found. Creating directory: '+directory)
    return directory

#Fast5 to pod5 conversion if fast5 files inputted
def cod5_to_fast5(fast5_input, pod5_output):
    cmd = 'pod5 convert fast5 '+fast5_input+'/*.fast5 -o '+pod5_output
    os.system(cmd)
