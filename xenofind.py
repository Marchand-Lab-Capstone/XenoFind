######################################################
'''

xenofind.py

By: J. Sumabat, N. Lai, Y. Huang, S. Peck 

'''
######################################################

######################################################

import argparse, textwrap
import os 
import sys



import argparse, textwrap
import os 
import sys 
from lib.xf_tools import *
from lib.xf_params import *


parser = argparse.ArgumentParser(
        usage='"python xenofind.py  [-h] {consensus, low_qual}',
        formatter_class=argparse.RawTextHelpFormatter,
	description= textwrap.dedent('''



______________________________________________________________________________
______________________________________________________________________________
_/\/\____/\/\_________________________________________________________________
___/\/\/\/\______/\/\/\____/\/\/\__/\/\______/\/\/\____/\/\__/\/\__/\/\/\_____
_____/\/\______/\/\/\/\/\__/\/\/\/\/\/\/\__/\/\__/\/\__/\/\/\/\________/\/\___
___/\/\/\/\____/\/\________/\/\__/\__/\/\__/\/\__/\/\__/\/\________/\/\/\/\___
_/\/\____/\/\____/\/\/\/\__/\/\______/\/\____/\/\/\____/\/\________/\/\/\/\/\_
______________________________________________________________________________
______________________________________________________________________________


********** Xemora : An XNA sequencing neural network trainer *********

Xemora is a tool pipeline used for nanopore sequencing of alternative basepairs (XNAs) that latches onto Remora 2.0 (ONT). This toolkit incorporates ONT-workflows to preprocess fast5, pod5, fastq, bam, and bed files for training a remora model. XNA sequence handling is done using the xFASTA format. Therefore, FASTA reference inputs should contain XNA bases in sequence lines. Models for basecalling can be trained on specific set of sequences (rather than entire sequence space). For optimal implementation, Xemora should be trained on at least two datasets, with and without the XNA substitutions. Xemora models can be exported and used as remora models for guppy, dorado, or bonito basecalling. Alternatively, Xemora can handle reference-based XNA identification directly. The general pipeline consists of two steps: 1) Training xemora on a reference set of reads with and without XNAs 2) basecalling using the trained model. 


Xemora command groups (additional help available within each command group):
	train		[Train] a xemora model on a set of input fast5 reads, localized to reference fasta containing XNAs. 
	basecall	[Basecall] a fast5 reads around XNA using a previously trained xemora model. 
         '''))


subparsers = parser.add_subparsers(dest='subparsers')



#Consensus Formation
parser_train = subparsers.add_parser('consensus', help='[-w working_dir] [-f raw_dir] [-r placeholder_fasta]')
parser_train.add_argument('-w',metavar = '[working_dir]', type=str,required = True, help='Path to output directory for storing output intermediates, temp files, and models.')
parser_train.add_argument('-f',metavar ='[raw_dir]', type=str,required = True, help='Path to input directories containing multi-fast5 or pod5 folder')
parser_train.add_argument('-r',metavar = '[placeholder_fasta]', type=str, required = True, help='Path to FASTA (.fa, .fasta) file of sequence or sequences with barcodes for alignment and randomer placeholders for consensus formation. ')


#Low Quality XNA Detection
parser_basecall = subparsers.add_parser('low_qual', help='[-w working_dir] [-f raw_dir] [-r consensus_fasta] ')
parser_basecall.add_argument('-w',metavar = '[working_dir]', type=str,required = True, help='Path to output directory for storing output intermediates and basecall results.')
parser_basecall.add_argument('-f',metavar ='[raw_dir]', type=str,required = True, help='Path to input directories containing multi-fast5 folders of XNA-containing sequence.')
parser_basecall.add_argument('-r',metavar = '[con_fasta]', type=str, required = True, help='Path to FASTA (.fa, .fasta) file of sequence or sequences with XNAs (e.g. BSPZKXJV). Should be same sequence context as xemora model training.')



args = parser.parse_args()
args.subparsers
exit_flag = False

#Print help if no arguments provided
if args.subparsers==None: 
    parser.print_help()
    sys.exit(0)

if args.subparsers == 'consensus': 
    if os.path.exists(os.path.expanduser(args.f))==False:
        print('XenoFind [ERROR] - Consensus requires either a Pod5 or Fast5 directory. This file path is invalid. Check to ensure pod5/fast path is correct.')
        exit_flag = True 

    if os.path.exists(os.path.expanduser(args.r))==False or os.path.exists(os.path.expanduser(args.r[1]))==False:
        print('XenoFind [ERROR] - Consensus requires a dummy fasta to be inputted. This file path is invalid. Check to ensure that the fasta path is corect.')
        exit_flag = True 


    if exit_flag == False: 
        cmd = 'python lib/consensus.py '+args.w+' '+args.f+' '+args.r
        os.system(cmd)
    else: 
        print('XenoFind [ERROR] - At least one file path not properly set. XenoFind exiting.')
        sys.exit()


#Guppy paths and stuff - add that to xr-parms
if args.subparsers == 'low_qual': 
    if os.path.exists(os.path.expanduser(args.f))==False:
        print('XenoFind [ERROR] - Inputted fast5/pod5 directory path not found or is not valid. Check to ensure path to fast5/pod5 directory is correct.')
        exit_flag = True 

    if os.path.exists(os.path.expanduser(args.r))==False:
        print('XenoFind [ERROR] - Inputted consensus fasta path is not valid. Check to ensure that the fasta path is correct')
        exit_flag = True 


    if exit_flag == False: 
        cmd = 'python lib/zf_low_qual.py '+args.w+' '+args.f+' '+args.r
        os.system(cmd)
    else: 
        print('XenoFind [ERROR] - At least one file path not properly set. Xemora basecaller exiting.')
        sys.exit()

