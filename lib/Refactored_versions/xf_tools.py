#####################################################################################
''' 
xf_tools.py

Useful functions to be called on  by other scripts.

By: By: J. Sumabat, N. Lai, Y. Huang, S. Peck, H. Kawabe, N. Kaplan, J. A. Marchand

'''
#####################################################################################
import os
from pathlib import Path

#Scans directory and subdirectory to get proper fast5 file path. Not explicitly handled with pod5 commands
def get_fast5_subdir(fast5_dir): 
    """
    get_fast5_subdir is a function that takes in a fast5 containing directory as its input and 
    extracts that path to be given to other functions such as cod5_to_fast5. The pathway extracted 
    is the 'raw_dir' in the string datatype from either terminal
    calling xenofind.py or from xenofind_pipe.py
    """
    path = os.path.normpath(os.path.expanduser(fast5_dir))
    if os.path.exists(path):
        fast5_files = list(Path(path).rglob("*.fast5" ))
        if len(fast5_files)>0:
            fast5_subdir = os.path.dirname(fast5_files[0])
            print('XenoFind [STATUS] - Found '+str(len(fast5_files))+' fast5 files in '+fast5_dir)
            return fast5_subdir
        else: 
            print('XenoFind [ERROR] - Could not find Fast5 files in specified directory. Check .fast5 exist.')
            return False
    else: 
        print('XenoFind [ERROR] - Could not find Fast5 directory. Check path')
        return False

#Scans directory and subdirectory to get proper pod5 file path. 
def get_pod5_subdir(pod5_dir): 
    """
    get_pod5_subdir is a function that takes in a pod5 containing directory as its input and 
    extracts that path to be given to other functions such as pod5_merge. The pathway extracted 
    is the xna_raw_dir, dna_raw_dir, or bc_raw pathways in the string datatype from either terminal
    calling xenofind.py or from xenofind_pipe.py
    """
    path = os.path.normpath(os.path.expanduser(pod5_dir))
    if os.path.exists(path):
        pod5_files = list(Path(path).rglob("*.pod5" ))
        if len(pod5_files)>0:
            pod5_subdir = os.path.dirname(pod5_files[0])
            print('Xemora [STATUS] - Found '+str(len(pod5_files))+' POD5 files in '+pod5_dir)
            return pod5_subdir
        else: 
            print('XenoFind [ERROR] - Could not find POD5 files in specified directory. Check .pod5 exist.')
            return False
    else: 
        print('XenoFind [ERROR] - Could not find POD5 directory. Check path')
        return False

#Check if working directory exists, if not create it. 
def check_make_dir(directory):
    """
    check_make_dir is a function that takes in a directory in the form of a string to actually generate that directory 
    by using the built in os.makedirs() function. The string that is generated is either the inputted working directory 
    string from the user or other predetermined directories that are put into the working directory. There are added print statements to 
    help users know what directory is being generated. 
    """
    directory = os.path.expanduser(directory)
    if not os.path.exists(directory):
        os.makedirs(directory)
        print('Xemora [STATUS] - Required directory not found. Creating directory: '+directory)
    return directory

#Fast5 to pod5 conversion if fast5 files inputted
def cod5_to_fast5(fast5_input, pod5_output):
    """
    cod5_to_fast5 is a function that calls the pod5 package from ONT to convert and merge a fast5 dataset to a single pod5
    file within the desired pod5 directory. This is generally the generated pod5 directory in the working directory chosen 
    by the user. 
    """
    cmd = 'pod5 convert fast5 '+fast5_input+'/*.fast5 -o '+pod5_output
    os.system(cmd)

#Merge pod5 files
def pod5_merge(pod5_input, merged_pod5):
    """
    pod5_merge is a function that calls the pod5 function from ONT to merge a pod5 dataset into a single pod5 file within the desired
    pod5 directory. This is generally the generated pod5 directory in the working directory chosen by the user.
    """
    cmd = 'pod5 merge '+ pod5_input+'/*.pod5 -o ' + merged_pod5
    os.system(cmd)
    
#BAM file primary read filter 
def filter_primary_alignments(input_bam, output_bam):
    """
    filter_primary_alignments is a function that takes in a bam file (generally a merged
    bam file after basecalling) and uses samtools and pysam to iterate through the bam file. 
    During iteration, reads that are not aligned, secondary aligned, or supplementary are 
    ignored. Read that are primary alignments are written to a new bam file, creating a 
    primary alignments only bam file. Samtools is used to 'index' the bam file, which allows
    pysam to iterate through it.
    """
    cmd = 'samtools index ' + input_bam #creates index file for BAM
    os.system(cmd)
    with pysam.AlignmentFile(input_bam, "rb") as infile, \
         pysam.AlignmentFile(output_bam, "wb", header=infile.header) as outfile:

        for read in infile:
            if not read.is_secondary and not read.is_supplementary and not read.is_unmapped:
                outfile.write(read)

    print(f"Filtered BAM file saved to {output_bam}")


def fetch_xna_pos(xm_header):
    pos=xm_header[xm_header.find('XPOS[')+5:-1].split('-')
    xpos = [x.split(':') for x in pos]
    return xpos


def xna_base_rc(xna_base, xna_bp): 
	for x in xna_bp: 
		if xna_base in x:
			if len(x)>1:
			    xx=list(x)
			    xx.remove(xna_base)
			    xrc=xx[0]
			else:
			   xrc=False
	return xrc


def check_xfasta_format(xfasta_file,standard_bases): 
    xfasta_header=False
    xna_in_sequence = False 
    with open(xfasta_file, "r") as fh:
        for line in fh:
            #Get header
            if line[0]=='>':
                header = line
                if 'XPOS[' in header: 
                    xfasta_header=True
                
            #Get sequence
            if line[0]!='>':
            
                #Make all upper case
                uline = line.upper()
                
                #Look for non-standard base
                diff = list(set(uline.replace('\n',''))-(set(standard_bases)))
                
                #This sequence contains an XNA
                if len(diff)>0:
                    xna_in_sequence = True 
    if xfasta_header==True and xna_in_sequence == False: 
        return True 
    else: 

        return False 
