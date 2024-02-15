import os
import glob
import sys
from xd_params import * 

#Initialize
#working_dir = os.path.expanduser(sys.argv[1]) #will be final impolemented version, will be testing in large scale script later 
#xna_raw_dir = os.path.expanduser(sys.argv[2])


#Initialize (manual input of datasets) 
working_dir = ''
raw_dir = ''

# Generate Required Directories
working_dir = check_make_dir(working_dir)
bam_dir = check_make_dir(os.path.join(working_dir, 'bam/'))
fastq_dir = check_make_dir(os.path.join(working_dir, 'fastq/'))


if basecall_pod == True: 
    cmd = os.path.exanduser(basecaller_path)+ ' hac ' + raw_dir + ' > '+os.path.join(bam_dir, 'bc.bam')
    os.system(cmd)

if bam_to_fastq == True: 
    cmd = 'samtools bam2fq ' + os.path.join(bam_dir,'bc.bam') + ' > ' + os.path.join(fastq_dir, 'bc.fastq')
    os.system(cmd) 


