'''
aggregate_reads.py 
J. Sumabat, N. Lai, S. Peck, Y. Huang
5/21/2024 -- Created by S. Peck

aggregate_reads.py is designed around importing a bamfile of reads, a consensensus fasta, and a pod5 file, 
and then combining the information contained into separate json files organized by consensus sequence and constituent reads.
This was penned with the hope of avoiding/mitigating memory leakage from a variety of sources. 

Note: Multiprocessing can handle generators per https://docs.python.org/3/library/multiprocessing.html, have not implemented.

Classes:
Read - stores read id, signal, frequency of given read
JointBamReads - stores reference name and iterator of Reads for a reference
ConsensusReadsData - stores all information that could be saved to a json file, as well as the method to do so.

Methods:
load_bam() - uses pysam to load relevant information from a bamfile.
load_fasta() - reads a fasta of consensus sequences.
map_signal_to_moves() - map a nanopore sequencing signal to the corresponding moves table extracted from that read's bamfile.
iterable_reads() - creates an iterable generator of Read objects constructed from the pod5 file.
aggregate_reads_to_reference() - generates a list of JointBamReads based on a pod5 path and bam data.
iterable_merged_data() - create a generator of ConsensusReadsData objects from alignment data and the fasta data
main() - run all steps needed to convert bam, pod5, and fasta to several json at output dir.
'''
import os
import sys
import pod5
import pysam
import numpy as np
import pandas as pd
import math
import json
from Bio import SeqIO
import re
import multiprocessing
import datetime
import argparse
from pathlib import Path
from alive_progress import alive_bar
import feature_extraction as fe
import setup_methods
import shutil

VERBOSE = False
N_SAVED_TO_JSON = 0
N_TOTAL_CONSENSUS = 0
LAST_TIME = datetime.datetime.now()
TOTAL_ELAPSED = 0
store_json = False
store_parquet = False

class Read: 
    '''
    class Read stores the read id, signal, and frequency of a given read.
    '''
    def __init__(self, rid, sig, freq):
        '''
        assign passed variables from constructor to public variables for the read.
        Parameters:
        rid: read id, can be a string.
        sig: signal extracted from pod5 as array or list.
        freq: sampling frequency from pod5 read.
        
        Returns:
        N/A
        '''
        self.read_id = rid
        self.sig = np.array(list(sig), dtype=np.float32)
        self.freq = freq
    def __reduce__(self):
        return Read, (self.read_id, self.sig, self.freq,) 
          
            
class JointBamReads:
    '''
    class JointBamReads stores the reference name and read iterator for a given reference sequence.
    '''
    def __init__(self, ref_name, readiterator):
        '''
        assign the passed variables from constructor to public variables for the
        reference sequence.
        
        Parameters:
        ref_name: reference sequence the reads belong to, as str
        readiterator: a generator object of Read objects.
        
        Returns:
        N/A
        '''
        self.ref_name = ref_name
        self.read_iterator = readiterator
        global N_TOTAL_CONSENSUS
        N_TOTAL_CONSENSUS += 1
        
    def __reduce__(self):
        return JointBamReads, (self.ref_name, self.read_iterator,) 
    

class ConsensusReadsData:
    '''
    class ConsensusReadsData stores the reference name, sequence, frequency, and list of reads data.
    '''
    def __init__(self, ref_name, reference_sequence, freq, reads_data_list):
        '''
        Generates an object with publicly available reference name, sequence, frequency, and reads data
        Parameters:
        ref_name: reference name as str
        reference_sequence: reference sequence id as str
        freq: sampling frequency.
        reads_data_list: list of dictionaries of read information
        
        Returns:
        N/A
        '''
        self.ref_name = ref_name.replace(":", "#")
        self.ref_seq = reference_sequence
        self.freq = freq
        self.reads_data = reads_data_list
        self.batched_assembled_features = []
        
    def __reduce__(self):
        return ConsensusReadsData, (self.ref_name, self.ref_seq, self.freq, self.reads_data,) 

        
    def export_to_json(self, path):
        '''
        export_to_json takes in a directory path and saves the ConsensusReadsData object
        to a json file at that path.
        Parameters:
        path: path to directory, as str.
        
        Returns:
        json object that was saved at the path.
        '''
        
        
        # generate the header from the reference sequence and frequency
        header = {'ref_seq':self.ref_seq, 'freq':self.freq}
        
        # generate the filename from the reference name
        filename = "cons_merge__{}.json".format(self.ref_name)
        
        # print the location the file is going to be saved
        #if VERBOSE: print(filename, end='\r')
        
        # generate the filepath the json will be saved at
        filepath = "/"+path[1:] + filename
        
        # generate the dataframe to json of the reads data
        json_frame = pd.DataFrame(self.reads_data).to_json()
        
        # convert the json export of the dataframe back to a dataframe again in the format json likes
        json_frame = json.loads(json_frame)
        
        # convert the header and json frame to a new json object
        json_out = json.dumps([header, json_frame])

        # remove any existing json with that filepath
        if os.path.exists(filepath):
            os.remove(filepath)

        if len(self.reads_data) > 0: 
            # create and open a json with that filepath
            with open(filepath, 'a') as file:
                # write the contents of the json
                file.write(json_out)
            
            # close the json
            file.close()
        # estimate the runtime.
        global N_SAVED_TO_JSON 
        global LAST_TIME
        global TOTAL_ELAPSED
        N_SAVED_TO_JSON += 1
        curtime = datetime.datetime.now()
        duration = (curtime-LAST_TIME).total_seconds()
        TOTAL_ELAPSED += duration
        mean_dt = TOTAL_ELAPSED/N_SAVED_TO_JSON
        est_time = datetime.timedelta(seconds = mean_dt*N_TOTAL_CONSENSUS)
        
        if VERBOSE: print("[ aggregate_reads.py {} ] {}/{} saved. Estimated Time: {}. Elapsed: {}.".format(datetime.datetime.now(), N_SAVED_TO_JSON, N_TOTAL_CONSENSUS, est_time, datetime.timedelta(seconds=TOTAL_ELAPSED)), end='\r')
        LAST_TIME = curtime

        # return the json object that was saved
        return json_out
    
    
    def feature_extraction(self, batch_size):
        '''
        feature extraction runs the feature extraction methods on the
        data stored within the current object.
        
        Parameters:
        batch_size: n of reads per batched consensus, as an int.
        
        Returns:
        a list of the dictionaries of assembled features.
        '''
        # read the alignment data.
        read_alignment_data = fe.extract_alignment_from_dict(self.reads_data, self.ref_seq)
        self.reads_data = None
        # batch the data.
        batched_reads = fe.batch_read_alignment_data(read_alignment_data, batch_size)

        # iterate through the batches and accumulate the parameters to pass.
        prepped_batches = []
    
        for i in range(len(batched_reads)):
            batch = batched_reads[i]
            prepped_batches.append((batch, self.ref_seq, self.freq, self.ref_name))
        
        # iterate through the parameters and process the features.
        self.batched_assembled_features = []
        with multiprocessing.Pool(processes=len(batched_reads)) as feat_gen_pool:
            self.batched_assembled_features = feat_gen_pool.starmap(fe.extract_batch_features, prepped_batches)
        
        # return the features. 
        return self.batched_assembled_features
    
    def save_batches_to_parquet(self, output_path):
        '''
        save_batches_to_parquet takes in an output superdirectory,
        and saves the batched consensus features as parquets to that directory.
        
        Parameters:
        output_path: path, as str, to superdirectory for parquet files.
        
        Returns:
        N/A
        '''

        # get the total number of batches
        total_batches = len(self.batched_assembled_features)
        
        # loop through each batch
        for i in range(total_batches):
            
            # get the batch of features
            batched_consensus = self.batched_assembled_features[i]
            
            # create the filepath name
            output_filepath = output_path + self.ref_name + "/"+ "consensus{}.parquet".format(i)
            
            # check if the path to the file exists
            setup_methods.check_make_dir(output_path + self.ref_name + "/", False)
            # remove any existing parquet directory with that filepath
            if os.path.exists(output_filepath):
                os.remove(output_filepath)

            # save the batch to parquet.
            batched_consensus.to_parquet(output_filepath)
           
        # estimate the runtime ----
        global N_SAVED_TO_JSON 
        global LAST_TIME
        global TOTAL_ELAPSED
        N_SAVED_TO_JSON += 1
    
        
        curtime = datetime.datetime.now()
        duration = (curtime-LAST_TIME).total_seconds()
        TOTAL_ELAPSED += duration
        mean_dt = TOTAL_ELAPSED/N_SAVED_TO_JSON
        est_time = datetime.timedelta(seconds = mean_dt*N_TOTAL_CONSENSUS)
        
        if VERBOSE: print("[ aggregate_reads.py {} ] {}/{} saved to parquet batches. Estimated Time: {}. Elapsed: {}.".format(datetime.datetime.now(), N_SAVED_TO_JSON, N_TOTAL_CONSENSUS, est_time, datetime.timedelta(seconds=TOTAL_ELAPSED)), end='\r')
        LAST_TIME = curtime
        # end runtime estimation ----
    

def percenttobar(frac):
    
    """
    converts a fractional value to a loading progress bar string.
    """
    bar_str = "|"
    max_bars = 20
    perc = frac*2000
    n_bars = int(perc/100)
    for i in range(n_bars):
        bar_str += "="
    for i in range(max_bars-n_bars):
        bar_str += " "
    bar_str += "|  {}%           ".format(round(frac*100, 3))
    return bar_str
    
    
def load_bam(bamfile_path):
    """
    load_bam takes in a path to a basecalled and aligned consensus bamfile,
    and converts the data into a list of dicts format.

    Parameters:
    bam_path: path, as str, to the basecalled and aligned consensus bamfile.

    Returns:
    list of dict of relevant data. 
        [sequence id, sequence, quality, length, reference index]
    """

    # Set an empty list to hold the data
    data_list = []
    n_reads = 4858652 /2
    reads_read = 0
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
            reads_read += 1
            print(percenttobar(reads_read/n_reads), end='\r')

    # convert the list of dicts into a dataframe
    #df = pd.DataFrame(data_list)

    # convert the qualities to list.
    #df['quality'] = df['quality'].apply(list)
    bamfile.close()

    # return the list of dicts
    return data_list


def load_fasta(consensus_fasta_path):
    """
    load_fasta takes in a consensus fasta path, 
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
        bar()
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
        del(moves_table)
        del(offset)
        del(signal)

        # return the list of observed signal segments
        return signal_list

    
def iterable_reads(p5_path, valid_reads):
    '''
    iterable_reads uses a path to a pod5 file and a list of valid read ids,
    and returns an iterable generator of Read objects constructed from the pod5 file.
    Parameters:
    p5_path: path, as str, to pod5 file.
    valid_reads: a list of read ids that exist in the pod5 file, as str.
    
    Returns:
    generator of Read Objects based on the valid reads
    
    Yeilds:
    Read Object from the generator
    '''
    
    # open the pod5 file
    with pod5.Reader(Path(p5_path)) as pod5_file:
        
        # loop through the reads that are valid
        for read in pod5_file.reads(selection=valid_reads, preload=["samples"]):
            
            # generate a read object and return it through Yield.
            out_read = Read(read.read_id, read.signal_pa, read.run_info.sample_rate)
            yield out_read

            
def aggregate_reads_to_reference(p5_path, alignment_data):
    '''
    aggregate_reads_to_reference generates a list of JointBamReads based on a pod5 path and bam data.
    Parameters:
    p5_path: path to pod5 file, as str
    alignment_data: data extracted from bamfile using load_bam.
    
    Returns:
    a list of JointBamRead objects that are the intersection of the valid reads in the pod5 and bam.
    '''
    # convert the bam data to a dataframe for ease of processing
    temp_df = pd.DataFrame(alignment_data)
    
    # loop through the pod5 file to get the set of the read ids
    with pod5.Reader(Path(p5_path)) as pod5_file:
        pod5_read_ids = set(pod5_file.read_ids)

    # create an empty list of aggregate references
    aggregate_refs = []
    
    # create a list of the bam reference sequence ids
    bam_ref_ids = list({x['ref_name'] for x in alignment_data})
        
    # loop through the uniqe references
    with alive_bar(len(bam_ref_ids)) as bar:
        
        # for each unique reference in the bam_reference ids
        for unique_ref in bam_ref_ids[:]:
            #print(unique_ref, end='\r')
            # get the set of unique references 
            seqs = set(temp_df[temp_df['ref_name'] == unique_ref]['seq_id'])

            # get the intersection of valid reads in both the bam and pod5
            valid_reads_in_ref = list(pod5_read_ids.intersection(seqs))

            # append the JointBamRead object for the unique reference
            aggregate_refs.append(JointBamReads(unique_ref, iterable_reads(p5_path, valid_reads_in_ref)))
            bar()
        

    # return the list of aggregate references.
    return aggregate_refs
    

def iterable_merged_data(aggregate_refs, alignment_data, consensus_data):
    '''
    create an generator of ConsensusReadsData objects from alignment data and the fasta data, as well as a
    list of JointBamReads objects.
    
    Parameters:
    aggregate_refs: list of reads both in bam and pod5
    alignment_data: data extracted from bamfile from load_bam
    consensus_data: data extracted from fasta file from load_fasta
    
    Returns:
    generator object of ConsensusReadsData
    
    Yeilds:
    ConsensusReadsData for each reference sequence in the fasta.
    '''
    
    # create a dataframe of the bam and fasta data
    bam_df = pd.DataFrame(alignment_data)
    fasta_df = pd.DataFrame(consensus_data)
    
    # set a variable to store sampling frequency -- assumed same for every read in a reference
    freq = None
    
    # for each reference JointBamRead object in the aggregate references:
    for ref in (aggregate_refs):
        
        # get the reference name as str
        ref_name = str(ref.ref_name)
        
        # get ther eference sequence by filtering the fasta by the reference name
        reference_sequence = fasta_df[fasta_df['id'] == ref_name].iloc[0]['seq']
        
        # generate an empty list to hold the data
        dict_list = []
        
        # loop through each read in the readiterator of the jointbamread
        for read in ref.read_iterator:
            
            # get the sequence id as str
            seq_id = str(read.read_id)

            # get the bamfile's read info for the sequence (Assumes no duplicate reads in a reference)
            tmp_df = bam_df[bam_df['ref_name']==ref_name]
            bam_read_info = tmp_df[tmp_df['seq_id'] == seq_id].iloc[0]

            # get the sequence, quality, cigar string, length, reference,
            #reversibility, moves, signal length, trim offset, sampling frequency,
            # signal, positional offset to reference. 
            sequence = bam_read_info['sequence']
            quality = np.asarray(bam_read_info['quality'], dtype='int16')
            cigar = bam_read_info['cigar']
            length = np.int16(bam_read_info['len'])
            ref = np.int16(bam_read_info['ref'])
            rev = np.int16(bam_read_info['rev'])
            moves = np.asarray(bam_read_info['moves'], dtype='int16')
            sig_len = np.int16(bam_read_info['sig_len'])
            trim_ofs = np.int16(bam_read_info['trim_ofs'])
            freq = read.freq
            read_signal = read.sig
            offset = bam_read_info.trim_ofs
            
            # map the signal to the moves table. 
            mapped_signal = map_signal_to_moves(moves, offset, read_signal)

            # put all that data in a dictionary
            read_dict = {'seq_id':seq_id,
                         'sequence':sequence,
                         'quality':quality,
                         'cigar':cigar,
                         'len':length,
                         'ref':ref,
                         'rev':rev,
                         'moves':moves,
                         'sig_len':sig_len,
                         'trim_ofs':trim_ofs,
                         'signal': mapped_signal}
            
            #append it to the list of dictionaries. 
            dict_list.append(read_dict)
            
        # yeild a consensusreadsdata object built off of the reference sequence and constituent reads.
        yield ConsensusReadsData(ref_name, reference_sequence, freq, dict_list)
        


def main(p5_path, bam_path, fasta_path, out_path, export, batchsize = 100, verb = False):
    '''
    main() serves to act as a method for processing each step in the process of extracting reads to
    constituent json files.
    
    Parameters:
    p5_path: path, as str, to pod5 file
    bam_path: path, as str, to bamfile from basecalling
    fasta_path: path, as str, to reference fasta
    out_path: path, as str, to directory that will hold json files
    export: boolean (default false) for if the output should be saved. Currently does not work unless run as script, as no JSON or PARQUET is specified.
    batchsize: size of the batch for multiprocessing the read aggregation. defualt of 100.
    verb: flips verbosity variable if True
    
    Returns:
    generator object for the consensuses, or path to saved files if exporting.
    '''
    global VERBOSE
    if verb: 
        VERBOSE = True
    
    # keep user updated on process.
    print("[ aggregate_reads.py {} ] Process started.".format(datetime.datetime.now()))
    if VERBOSE: print("[ aggregate_reads.py {} ]  Reading Bam at {}...".format(datetime.datetime.now(),bam_path), end = '\r')
    
    # load the bam data
    bam_data = load_bam(bam_path)
    if VERBOSE: print("[ aggregate_reads.py {} ] Bam read.                                                                                                                           ".format(datetime.datetime.now()))
    
    #  read the fasta data.
    if VERBOSE: print("[ aggregate_reads.py {} ] Reading fasta at {}...".format(datetime.datetime.now(), fasta_path), end="\r")
    fasta_data = load_fasta(fasta_path)
    if VERBOSE: print("[ aggregate_reads.py {} ] Fasta read.                                                                                                                           ".format(datetime.datetime.now()))
    
    # aggregate the bam read data from the pod5 file.
    if VERBOSE: print("[ aggregate_reads.py {} ] Aggregating bam read data from pod5 file at {}...".format(datetime.datetime.now(), p5_path), end="\r")

    list_of_joint_bam_reads = aggregate_reads_to_reference(p5_path, bam_data) ######
    if VERBOSE: print("[ aggregate_reads.py {} ] Aggregation complete.                                                                                                                           ")
    
    # read the batched bam reads
    if VERBOSE: print("[ aggregate_reads.py {} ] Batching the data...".format(datetime.datetime.now()), end='\r')
    batched_bam_reads = [list_of_joint_bam_reads[i:i+batchsize] for i in range(0, len(list_of_joint_bam_reads), batchsize)]
    if VERBOSE: print("[ aggregate_reads.py {} ] Data batched.                            ".format(datetime.datetime.now()))
    
    # merge the data to an iterable
    if VERBOSE: print("[ aggregate_reads.py {} ] Generating Consensus Reads Data...".format(datetime.datetime.now()), end='\r')
    iterable_merged = iterable_merged_data(list_of_joint_bam_reads, bam_data, fasta_data)
    if VERBOSE: print("[ aggregate_reads.py {} ] Consensus Reads Data Generated.                                                                                                             ".format(datetime.datetime.now()))
    
    LAST_TIME = datetime.datetime.now()

    export_dir = os.listdir(out_path)
    if export:
        if VERBOSE: print("[ aggregate_reads.py {} ] Exporting to directory at {}...".format(LAST_TIME, out_path))
        for read in iterable_merged:
            
            if N_SAVED_TO_JSON == 0:
                LAST_TIME = datetime.datetime.now()
            
            if store_json:
                if read.ref_name + '.json' not in export_dir:
                    read.export_to_json(out_path)
                    
            elif store_parquet:
                if read.ref_name not in export_dir:
                    read.feature_extraction(batchsize)
                    read.save_batches_to_parquet(out_path)
        if VERBOSE: print("[ aggregate_reads.py {} ] Export complete. Have a nice day! :)                               ".format(datetime.datetime.now()))
        return out_path
    else:
        return iterable_merged
                
        
if __name__ == '__main__':
    
    # setup the args
    parser = argparse.ArgumentParser(description='Merges Pod5, Bam, and Fasta Files into consensus-organized read compilations in Json. ')
    parser.add_argument('-y', help="Confirm output directory initialization", action="store_true")
    parser.add_argument('-v', "--verbose", help="Make output verbose.", action="store_true")
    parser.add_argument('-bam', help="Bam filepath.")
    parser.add_argument('-pod5', help="Pod5 filepath.")
    parser.add_argument('-fasta', help="Fasta filepath.")
    parser.add_argument('-output', help="Output filepath.")
    parser.add_argument('-json', help="flags output as JSON type", action = "store_true")
    parser.add_argument('-parquet', help="flags output as parquet type", action = "store_true")
    parser.add_argument('-batch_size', help="number of reads per batched consensus. Default=100")
    
    
    # Parse the arguments
    args = parser.parse_args()
    
    # check the arguments and assign them
    wipe_output = args.y
    VERBOSE = args.verbose
    bam_path = args.bam
    pod5_path = args.pod5
    fasta_path = args.fasta
    output_path = args.output
    store_json = args.json
    store_parquet = args.parquet
    export = False
    batchsize = args.batch_size
    
    if store_json:
        export = True
    elif store_parquet:
        export = True
    else:
        export = False
        
    # double check the args
    if type(bam_path) == type(None): bam_path = input("[ aggregate_reads.py ] Bam Path: ")
    if type(pod5_path) == type(None): pod5_path = input("[ aggregate_reads.py ] Pod5 path: ")
    if type(fasta_path) == type(None): fasta_path = input("[ aggregate_reads.py ] Reference fasta path: ")
    if type(output_path) == type(None) and export: output_path = input("[ aggregate_reads.py ] Output path: ")
    if type(batch_size) == type(None): batchsize = 100

    # Check if output is ok to be removed
    if (wipe_output == False) and (export):
        resp = input("[ aggregate_reads.py ] files in {} will be removed. Continue? [Y/N]:".format(output_path))
        if str(resp) == 'Y':
            wipe_output = True

    # wipe the output directory if selected
    if (wipe_output == True) and (export):
        shutil.rmtree(output_path)
    
    # run the main method.
    main(pod5_path, bam_path, fasta_path, output_path, export, batchsize)
    print("[ aggregate_reads.py {} ] Closing.".format(datetime.datetime.now()))
    sys.exit()

    
        