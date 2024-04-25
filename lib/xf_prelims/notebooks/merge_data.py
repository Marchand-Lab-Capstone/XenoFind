'''
Solely to avoid issues where TOO MUCH memory leaks, 
here's this to help subdivide that stuff down.

I pray it works.

-Sebastian
'''

import os
import pysam
import remora
import sys
sys.path.append('..//')
import shannon_entropies as sp
import pod5
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import math
import line_profiler as lprof


def load_pod5_data(p5_path):
    
    # For some godforsaken reason, pod5 somehow consumes upwards of 6 GiB of memory when this runs, 
    # and it DOES NOT DISSAPEAR. ESLIT'rkjasvkajs;flksarejva;lkrvj
    data_list = []
    with pod5.Reader(p5_path) as pod5_file:
        reads = pod5_file.reads()
        for read in reads:
            seq_id = read.read_id
            signal = read.signal_pa
            freq = read.run_info.sample_rate
            data_dict = {'seq_id': str(seq_id),
                         'signal': signal,
                         'freq': freq}
            data_list.append(data_dict)
    pod5_file.close()
            
    return data_list

#https://pod5-file-format.readthedocs.io/en/latest/reference/api/pod5.reader.html#pod5.reader.ReadRecord


def merge_bam_reads_by_id(bam_list_dict, read_list_dict, consensus_list_dict):
        
    # this code can be made more time efficient and memory-chunking but that's not important rn
    bld = pd.DataFrame(bam_list_dict)
    rld = pd.DataFrame(read_list_dict)
    cld = pd.DataFrame(consensus_list_dict)
    
    merged_data_df = pd.merge(bld, rld, on='seq_id', how='left')
    
    for i in range(len(cld)):
        #merged_data_df.at[]
        indexes = merged_data_df[merged_data_df['ref_name'] == cld['id'][i]].index
        merged_data_df.loc[indexes,'ref_seq'] = cld['seq'][i]
    
    
    data_dict = merged_data_df.to_dict('records')

    return data_dict


def map_signal_to_moves(moves_table, offset, signal):
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

        # return the list of observed signal segments
        return signal_list
    

def convert_signal_from_dict(m_list_dict):
    merged_data_df = pd.DataFrame(merged_list_dict)
    merged_data_df['signal'] = merged_data_df.apply(lambda x: map_signal_to_moves(x['moves'], x['trim_ofs'], x['signal']), axis=1)
    merged_data_df['signal'].apply(np.tolist)
    end_dict = merged_data_df.to_dict('records')
    del(merged_data_df)
    return end_dict
    

if __name__ == "__main__":
    
    bam_path = input("Bam Path: ")
    pod5_path = input("pod5 path: ")
    ref_fasta = input("reference fasta path: ")
    output_path = input("output path ")
    
    # Use the shannon entropy methods to load the bam
    loaded_bam = sp.load_in_data(bam_path)
    print("bam loaded")
    #Load the consensus reference fasta using shannon entropy methods
    consens = sp.consensus_formatter(ref_fasta)

    print('fasta loaded')
    
    print('processing reads...')
    # Load the pod5 read data using load_pod5_data -- MEMORY LEAK FROM pod5!!!!
    dl = load_pod5_data(pod5_path)

    print('reads processed')
    
    print('merging...')
    # Merge the data
    merged_list_dict = merge_bam_reads_by_id(loaded_bam, dl, consens)
    print('merged')
    
    print('mapping signals to moves...')
    merged_signal_list_dict = convert_signal_from_dict(merged_list_dict)
    print('signals mapped')
    
    print('saving to {}...'.format(output_path))
    
    out_df = pd.DataFrame(merged_signal_list_dict)
    out_df.groupby("ref_name").apply(lambda x: x.to_csv("{}/consensus_sequence_{}.csv".format(output_path,x.name)))
    print('{} Closing. '.format(result))
    sys.exit()