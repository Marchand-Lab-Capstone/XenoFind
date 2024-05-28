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

parser = argparse.ArgumentParser(
        usage='"python xenofind.py  [-h] {consensus, model_gen, find, full}',
        formatter_class=argparse.RawTextHelpFormatter,
	description= textwrap.dedent('''



______________________________________________________________________________


********** XenoFind: an XNA / modified base position detector *********


Xenofind command groups (additional help available within each command group):
	consensus Generate a consensus sequence in the form of a fasta file. 
	find Detect modification positions using consensus features and machine learning 
    model_gen Generate a new XenoFind model using consensus features 
    full Run both the consensus and find pipelines
         '''))


subparsers = parser.add_subparsers(dest='subparsers')



#Consensus Formation
parser_train = subparsers.add_parser('consensus', help='[-w working_dir] [-f raw_dir] [-r placeholder_fasta]')
parser_train.add_argument('-w',metavar = '[working_dir]', type=str,required = True, help='Path to output directory for storing outputs and intermediates')
parser_train.add_argument('-f',metavar ='[raw_dir]', type=str,required = True, help='Path to input directories containing multi-fast5 or pod5 folder')
parser_train.add_argument('-r',metavar = '[placeholder_fasta]', type=str, required = True, help='Path to FASTA (.fa, .fasta) file of sequence or sequences with barcodes for alignment and randomer placeholders for consensus formation. ')

#Consensus Formation
parser_train = subparsers.add_parser('model_gen', help='[-w working_dir] [-f raw_dir] [-r placeholder_fasta]')
parser_train.add_argument('-w',metavar = '[working_dir]', type=str,required = True, help='Path to output directory for storing outputs and intermediates')
parser_train.add_argument('-f',metavar ='[raw_dir]', type=str,required = True, help='Path to input directories containing multi-fast5 or pod5 folder')
parser_train.add_argument('-r',metavar = '[ref_fasta]', type=str, required = True, help='Path to FASTA (.fa, .fasta) file of sequence or sequences for model training. ')

#XNA Detection
parser_basecall = subparsers.add_parser('find', help='[-w working_dir] [-f raw_dir] [-r consensus_fasta] ')
parser_basecall.add_argument('-w',metavar = '[working_dir]', type=str,required = True, help='Path to output directory for storing outputs and intermediates')
parser_basecall.add_argument('-f',metavar ='[raw_dir]', type=str,required = True, help='Path to input directories containing multi-fast5 or pod5 files for consensus formed dataset')
parser_basecall.add_argument('-r',metavar = '[con_fasta]', type=str, required = True, help='Path to consensus fasta')

#Full Pipeline
parser_full = subparsers.add_parser('full', help='[-w working_dir] [-f raw_dir] [-r placeholder_fasta]')
parser_full.add_argument('-w',metavar = '[working_dir]', type=str,required = True, help='Path to output directory for storing outputs and intermediates')
parser_full.add_argument('-f',metavar ='[raw_dir]', type=str,required = True, help='Path to input directories containing multi-fast5 or pod5 folder')
parser_full.add_argument('-r',metavar = '[placeholder_fasta]', type=str, required = True, help='Path to FASTA (.fa, .fasta) file of sequence or sequences with barcodes for alignment and randomer placeholders for consensus formation. ')

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

    if os.path.exists(os.path.expanduser(args.r))==False:
        print('XenoFind [ERROR] - Consensus requires a dummy fasta to be inputted. This file path is invalid. Check to ensure that the fasta path is corect.')
        exit_flag = True 

    if exit_flag == False: 
        cmd = 'python lib/weighted_consensus_split.py '+args.w+' '+args.f+' '+args.r
        os.system(cmd)
    else: 
        print('XenoFind [ERROR] - At least one file path not properly set. XenoFind exiting.')
        sys.exit()

if args.subparsers == 'model_gen': 
    if os.path.exists(os.path.expanduser(args.f))==False:
        print('XenoFind [ERROR] - model_gen requires either a Pod5 or Fast5 directory. This file path is invalid. Check to ensure pod5/fast path is correct.')
        exit_flag = True 

    if os.path.exists(os.path.expanduser(args.r))==False:
        print('XenoFind [ERROR] - model_gen requires a forward and reverse sequence fasta files for training. At least one of these file path is invalid. Check to ensure fasta path(s) are set properly..')
        exit_flag = True 

    if exit_flag == False: 
        cmd = 'python lib/xf_model_training.py '+args.w+' '+args.f+' '+args.r
        os.system(cmd)
    else: 
        print('XenoFind [ERROR] - At least one file path not properly set. XenoFind exiting.')
        sys.exit()

if args.subparsers == 'find': 
    if os.path.exists(os.path.expanduser(args.f))==False:
        print('XenoFind [ERROR] - Inputted fast5/pod5 directory path not found or is not valid. Check to ensure path to fast5/pod5 directory is correct.')
        exit_flag = True 

    if os.path.exists(os.path.expanduser(args.r))==False:
        print('XenoFind [ERROR] - Inputted consensus fasta path is not valid. Check to ensure that the fasta path is correct')
        exit_flag = True 


    if exit_flag == False: 
        cmd = 'python lib/xf_xna_find.py '+args.w+' '+args.f+' '+args.r
        os.system(cmd)
    else: 
        print('XenoFind [ERROR] - At least one file path not properly set. Exiting.')
        sys.exit()

if args.subparsers == 'full': 
    if os.path.exists(os.path.expanduser(args.f))==False:
        print('XenoFind [ERROR] - Consensus requires either a Pod5 or Fast5 directory. This file path is invalid. Check to ensure pod5/fast path is correct.')
        exit_flag = True 

    if os.path.exists(os.path.expanduser(args.r))==False:
        print('XenoFind [ERROR] - Consensus requires a dummy fasta to be inputted. This file path is invalid. Check to ensure that the fasta path is corect.')
        exit_flag = True 

    if exit_flag == False: 
        cmd = 'python lib/weighted_consensus_split.py '+args.w+' '+args.f+' '+args.r
        os.system(cmd)
        '''
        doing a scuffed extract for now 
        '''
        consensus_fasta = os.path.join(args.w, 'consensus_files/consensus.fasta')
        cmd = 'python lib/xf_xna_find.py '+args.w+' '+args.f' '+consensus_fasta
        os.system(cmd)
    else: 
        print('XenoFind [ERROR] - At least one file path not properly set. XenoFind exiting.')
        sys.exit()
