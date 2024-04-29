'''
merge_consensus.py
J. Sumabat, N. Lai, S. Peck, Y. Huang
4/24/2024 -- Created by S. Peck
4/29/2024 -- Updated by S. Peck

merge_consensus.py is designed around importing a bamfile of reads, a consensensus fasta, and a pod5 file, 
and then combining the information contained into separate json files organized by consensus sequence and constituent reads.
This was penned with the hope of avoiding/mitigating memory leakage from a variety of sources. 

load_in_data() - uses pysam to load all the relevant information from a bamfile to a list of dictionaries.
consensus_formattter() - reads a consensus fasta and returns it as a list of dictionaries.
load_pod5_data() - reads a pod5 file and loads relevant information into a list of dictionaries.
merge_bam_reads_by_id() - merge data from load_in_data, consensus_formatter, and load_pod5_data by comparing
                          read id and consensus id.
map_signal_to_moves() - map a nanopore sequencing signal to the corresponding moves table extracted from that read's bamfile.
convert_signal_to_dict() - convert merged data's signal data to be mapped onto each read's corresponding moves.
save_by_consensus() - exports each individual consensus and its constituent reads to json files at a specified folder.
main() - request paths to bam, pod5, and reference fasta files, merge, and 
         save those data by consensus to individual json files
         at a specified output.


'''

# import statements
import os
import pysam
import remora
import sys
import pod5
import numpy as np
import pandas as pd
import math
import line_profiler as lprof
import json
from Bio import SeqIO
import tables
import h5py # not sure if nessecary, too afraid to check
import re
import numpy as np
import scipy.signal as sig
import gc



def load_in_data(bamfile_path):
    """
    load_in_data takes in a path to a basecalled and aligned consensus bamfile,
    and converts the data into a list of dicts format.

    Parameters:
    bam_path: path, as str, to the basecalled and aligned consensus bamfile.

    Returns:
    list of dict of relevant data. 
        [sequence id, sequence, quality, length, reference index]
    """

    # Set an empty list to hold the data
    data_list = []

    # Open the bamfile
    with pysam.AlignmentFile(bamfile_path, 'r') as bamfile:

        # Loop through each read in the bamfile
        for read in bamfile:

            # If the read is mapped, 
            if (not read.is_unmapped):

                # Get the sequence id, alignment sequence, cigar string, reference start position
                # and the quality, and put it into a dictionary
                sequence_id = read.query_name
                ref_name = read.reference_name
                seq = read.query_alignment_sequence
                cigar = read.cigartuples
                ref = read.reference_start
                quality = read.query_alignment_qualities
                rev = read.is_reverse
                tag_dict = dict(read.tags)
                read_dict = {'seq_id':str(sequence_id),
                             'ref_name':ref_name,
                             'sequence':seq,
                             'quality':list(quality),
                             'cigar':cigar,
                             'len':int(len(seq)),
                             'ref':int(ref),
                             'rev':rev,
                             'moves':tag_dict['mv'],
                             'sig_len': tag_dict['ns'],
                             'trim_ofs':tag_dict['ts']}
                #if (ref == 0):

                # Append the dictionary to the data list
                del(sequence_id)
                del(ref_name)
                del(seq)
                del(cigar)
                del(ref)
                del(quality)
                del(rev)
                del(tag_dict)
                data_list.append(read_dict)
            del(read)

    # convert the list of dicts into a dataframe
    #df = pd.DataFrame(data_list)

    # convert the qualities to list.
    #df['quality'] = df['quality'].apply(list)
    bamfile.close()

    # return the dataframe
    return data_list


def consensus_formatter(consensus_fasta_path):
    """
    consensus_formatter takes in a consensus fasta path, 
    and returns a list of dicts containing the ids and sequences of all valid samples.

    Parameters:
    consensus_fasta_path: a path, as str, to the consensus fasta
    
    Returns:
    a list of dicts containing id and sequence
    """

    # Set up list to hold consensus dictionaries
    valid_consensus_list = []

    # Search through the records in the fasta
    for record in SeqIO.parse(consensus_fasta_path, "fasta"):

        # For each record, get the header
        header = record.description

        # set the boundary conditions to invalid values before adding 'em
        cons_dict = {"start":None,"end":None, "id": None, "seq": None}

        # Check if there isn't a gap in the header
        if (not ('_Gap' in header)):

            # Get the boundary conditions by searching
            match = re.search(r'BC 1: (\d+), BC 2: (\d+)', header)
            
            # If the search has a match, add the indexes to the record
            if match:
                cons_dict["start"] = int(match.group(1))
                cons_dict["end"] = int(match.group(2))

            # add the id and sequence (with P swapped for G) to the dict
            cons_dict['id'] = str(record.id),
            cons_dict['seq'] = ''.join(record.seq.replace('P','G'))

            # add the dict to the valid consensus list
            valid_consensus_list.append(cons_dict)

    # convert the list of dicts to a dataframe
    df = pd.DataFrame(valid_consensus_list)
    
    # convert the ID from a weird touple to an actual value
    df[['id']] = pd.DataFrame(df['id'].tolist(), index=df.index)#3#######

    # return the dataframe as a list of dicts
    return df.to_dict('records')


def load_pod5_data(p5_path):
    '''
    load_pod5_data reads read information from a pod5 file and then exports it as
    a list of dicts.
    
    Parameters:
    p5_path: path to pod5 file, as a str.
    
    Returns:
    a list of dictionaries containing sequence id 'seq_id', signal 'signal', and sampling frequency 'freq'.
    '''
    # For some godforsaken reason, pod5 somehow consumes upwards of 6 GiB of memory when this runs, 
    # and it DOES NOT DISSAPEAR. ESLIT'rkjasvkajs;flksarejva;lkrvj
    data_list = []
    
    # use pod5's reader object to read through the pod5 file
    with pod5.Reader(p5_path) as pod5_file:
        
        # look at the reads
        reads = pod5_file.reads()
        
        # loop through each read
        for read in reads:
            
            # save the sequence id, signal, and frequency to a dict
            seq_id = read.read_id
            signal = read.signal_pa
            freq = read.run_info.sample_rate
            data_dict = {'seq_id': str(seq_id),
                         'signal': signal,
                         'freq': freq}
            
            # append the dict to the list
            data_list.append(data_dict)
            
            # force garbage collection
            del(seq_id)
            del(signal)
            del(freq)
            del(data_dict)
            # A GC collect here makes it run in circles
        del(read)
        gc.collect()
        
    pod5_file.close()
    del (pod5_file)
            
    # return the data list 
    return data_list

#https://pod5-file-format.readthedocs.io/en/latest/reference/api/pod5.reader.html#pod5.reader.ReadRecord


def merge_bam_reads_by_id(bam_list_dict, read_list_dict, consensus_list_dict):
    '''
    merge_bam_reads_by_id() takes in lists of dicts extracted from a bamfile, pod5 file, and consensus fasta,
    and then merges them into a single list of dictionaries.
    
    Parameters:
    bam_list_dict: list of dicts extracted using the load_in_data method.
    read_list_dict: list of dicts extracted using the load_pod5_data method.
    consensus_list_dict: list of dicts extracted using the consensus_formatter method.
    
    Returns:
    a list of dictionaries of all the data, merged by read sequence id.
    '''
    
    # this code can be made more time efficient and memory-chunking but that's not important rn
    # convert inputs to dataframe
    bld = pd.DataFrame(bam_list_dict)
    rld = pd.DataFrame(read_list_dict)
    cld = pd.DataFrame(consensus_list_dict)
    
    # merge the bam dataframe and read dataframe by sequence id
    merged_data_df = pd.merge(bld, rld, on='seq_id', how='left')
    
    # loop through the consensus fasta,
    for i in range(len(cld)):
        
        # combine the id and sequence if the reference name at the index is the same
        # filter the merged data by reference name
        indexes = merged_data_df[merged_data_df['ref_name'] == cld['id'][i]].index
        merged_data_df.loc[indexes,'ref_seq'] = cld['seq'][i]
        del(indexes)
    
    # convert the dataframe back to a list of dicts.
    data_dict = merged_data_df.to_dict('records')
    
    # force garbage collector
    del(bld)
    del(rld)
    del(cld)
    del(merged_data_df)
    gc.collect()

    return data_dict


def map_signal_to_moves(moves_table, offset, signal):
    '''
    map_signal_to_moves() takes in a moves table, base offeset, and a signal,
    and then maps the signal to each corresponding move/observation in the table.
    
    Parameters:
    moves_table: an array or list containing the stride as the first index, and the moves
                 table of a nanopore sequencing pysam read. 
    offset: the alignment base position, as an int. 
    signal: a list or array containing nanopore sequencing signal corresponding to the moves table passed.
    
    Returns:
    a list of lists of signals that correspond to each move in the moves table.
    '''
    if (str(signal) == 'nan'):
        return np.asarray([])
    else:
        # Generate an empty list to hold the true moves
        moves = []

        # get the raw moves, beginning after the first index
        raw_moves = list(moves_table[1:])

        # get the stride, beginning before the moves
        stride = moves_table[0]

        # Loop through the raw moves
        for i in range(len(raw_moves)):

            # if the raw move is zero
            if (raw_moves[i] == 0):

                # append n=stride zeroes
                for j in range(stride):
                    moves.append(0)

            # otherwise, append a 1 and then n-1 zeroes
            else:
                moves.append(1)
                for j in range(stride-1):
                    moves.append(0)

        # get the length difference between the signal+offset and the moves
        lendiff = (len(signal)+offset) - len(moves)

        # convert the moves to an array
        moves = np.asarray(moves)

        # pad the moves with zeroes for the offset and the length difference
        moves = np.pad(moves, (offset, lendiff), 'constant')

        # get the indexes where a stride step occurs
        stride_indicies = (np.where(moves == 1)[0])

        # add the last index as the length of the signal
        stride_indicies = np.append(stride_indicies, len(signal))

        # generate empty list to hold per-observation signals
        signal_list = []

        # loop through the stride indexes
        for i in range(len(stride_indicies)-1):

            # get beginning and end indexes
            begindex = stride_indicies[i]
            enddex = stride_indicies[i+1]

            # append the signal that falls between those indexes
            signal_list.append(signal[begindex:enddex])

        # force run the garbage collector
        del(raw_moves)
        del(stride)
        del(lendiff)
        del(moves)
        del(stride_indicies)
        
        # return the list of observed signal segments
        return signal_list
    

def convert_signal_from_dict(m_list_dict):
    '''
    convert_signal_from_dict takes in a list of dictionaries that were generated from merge_bam_reads_by_id,
    and then maps the signal of each read to the corresponding move table.
    
    Parameters:
    m_list_dict: list of dictionaries merged from the merge_bam_reads_by_id method.
    
    Returns:
    a list of dicts with the signal now mapped to the corresponding moves table.
    '''
    # convert to dataframe
    merged_data_df = pd.DataFrame(merged_list_dict)
    
    # map the signal to moves
    merged_data_df['signal'] = merged_data_df.apply(lambda x: map_signal_to_moves(x['moves'], x['trim_ofs'], x['signal']), axis=1)
    
    # convert back to a dict
    end_dict = merged_data_df.to_dict('records')
    
    # delete the dataframe
    del(merged_data_df)
    
    # run garbage collection
    
    # return the dict.
    return end_dict
    
    
def save_by_consensus(merged_dict, savefile_path):
    '''
    save_by_consensus() saves reads by consensus in their own individual json files.
    
    Parameters:
    merged_dict: list of dicts containing relevant data keys: [ref_name,
                                                               quality,
                                                               len,
                                                               ref,
                                                               rev,
                                                               moves,
                                                               sig_len,
                                                               trim_offs,
                                                               freq,
                                                               ref_seq]
    savefile_path: folder that the json files will be exported to, as a str. MAKE SURE IT ENDS WITH '/'.
    
    '''
    # First, make it a dataframe again
    df = pd.DataFrame(merged_signal_list_dict)
    
    # Make directory
    os.makedirs(savefile_path,exist_ok=True)
    
    # then, get the reference names
    consensus_ids = list(df['ref_name'].unique())
    
    # loop through the reference names:
    for consensus_id in consensus_ids:
        # filter the data accordingly:
        consensus = df[df['ref_name'] == consensus_id]
        
        # Condense the data types (Making sure no info is lost)
        consensus.loc[:,'quality']=consensus['quality'].apply(np.int16)
        consensus.loc[:,'len']=consensus['len'].apply(np.int16)
        consensus.loc[:,'ref']=consensus['ref'].apply(np.int16)
        consensus.loc[:,'rev']=consensus['rev'].apply(np.int16)
        consensus.loc[:,'moves']=consensus['moves'].apply(np.asarray)
        consensus.loc[:,'sig_len']=consensus['sig_len'].apply(np.int16)
        consensus.loc[:,'trim_ofs']=consensus['trim_ofs'].apply(np.int16)
        consensus = consensus.reset_index(drop=True)
        # sub method to convert datatypes in a list to a list
        
        def listconvert(x):
            # submethod listconvert is used to map subsects of data entries to lists.
            for i in range(len(x)):
                x[i]=list(x[i])
            return x
        
        consensus.loc[:,'signal']=consensus['signal'].apply(listconvert)
        consensus.loc[:,'freq']=consensus['freq'].apply(np.float16)
        
        # Get the reference name, sequence, and frequency - 
        # #### THIS IS ASSUMING THAT THE FREQUENCY IS THE SAME FOR ALL READS ####
        ref_name = consensus['ref_name'][0]
        ref_seq = consensus['ref_seq'][0]
        freq = consensus['freq'][0]
        
        # remove those columns
        del(consensus['ref_seq'])
        del(consensus['ref_name'])
        del(consensus['freq'])
        
        # generate the header and filename
        header = {'ref_seq':ref_seq,'freq':freq}
        filename = "cons_merge__{}.json".format(ref_name)
        
        # generate the filepath
        filepath = savefile_path + filename
        
        # generate the json object with the header and the dataframe
        json_frame = consensus.to_json() # conver dataframe to json
        json_frame = json.loads(json_frame) # convert json back to dict
        json_out = json.dumps([header, json_frame]) # convert list of dicts to json

        del(consensus)
        # if the path exists, overwrite it
        if os.path.exists(filepath):
            os.remove(filepath)
        
        # write the json object
        with open(filepath, 'a') as file:
            file.write(json_out)
        file.close()
        
        # run garbage collector.
        del(file)
        del(ref_name)
        del(ref_seq)
        del(freq)
        del(header)
        del(filename)
        del(filepath)
        del(json_frame)
        del(json_out)

    
    # return the savefile path
    return savefile_path    


if __name__ == "__main__":
    
    bam_path = input("Bam Path: ")
    pod5_path = input("pod5 path: ")
    ref_fasta = input("reference fasta path: ")
    output_path = input("output path ")
    
    # Use the shannon entropy methods to load the bam
    loaded_bam = load_in_data(bam_path)
    gc.collect()
    print("bam loaded")
    #Load the consensus reference fasta using shannon entropy methods
    consens = consensus_formatter(ref_fasta)
    gc.collect()

    print('fasta loaded')
    
    print('processing reads...')
    # Load the pod5 read data using load_pod5_data -- MEMORY LEAK FROM pod5!!!!
    dl = load_pod5_data(pod5_path)
    gc.collect()

    print('reads processed')
    
    print('merging...')
    # Merge the data
    merged_list_dict = merge_bam_reads_by_id(loaded_bam, dl, consens)
    gc.collect()
    print('merged')
    
    print('mapping signals to moves...')
    merged_signal_list_dict = convert_signal_from_dict(merged_list_dict)
    gc.collect()
    print('signals mapped')
    
    print('saving to {}...'.format(output_path))
    
    save_by_consensus(merged_signal_list_dict, output_path)
    gc.collect()
    print('Closing. ')
    sys.exit()