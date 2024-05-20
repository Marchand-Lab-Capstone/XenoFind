'''
feature_extraction.py
J. Sumabat, N, Lai, S. Peck, Y. Huang
4/29/2024 -- Created by S. Peck
5/5/2024 -- Updated by S. Peck
5/16/2024 -- Updated by S. Peck

feature_extraction.py contains methods useful for extracting possible ml features from
a consensus of reads and their aggregate information as produced by merge_consensus.py.
This is designed to operate on one json file at a time, and is hitherto incomplete with respect to
actually extracting all features. Presently, this includes base-level and quality-score-level feature extraction,
with signal extraction in escrow. 

Methods:
load_info_from_json() - method to load data from a json exported from merge_consensus.py
extract_coordinate_data() - method to extract coordinate data from the loaded json
align_observation_ops() - method to align a mapped signal to the string operations
shift_to_alignment() - method to shift the alignment of signals, base operations, and qualities
parse_read_aignment_info() - method to prepare read data to be aligned to reference
extract_alignment_from_dict() - method to extract the reference-based alignment data
calc_base_probs() - method to calculate the base observation probabilities
shannon_entropies() - method to calculate shannon entropies 
get_first_base() - method to get the first base in a list of bases
check_base_match() - method to check if a base matches the reference
mismatch_prob() - method to check probability of mismatches at given bases
is_in_del() - method to check if s+omething is an insertion/deletion or not.
filter_by_indel() - method to filter data to not conain information corresponding to insertions or deletions
remove_in_dels() - method to remove insertions and deletions
getmean() - method to get the mean value of dataset
mean_std_qual_by_base() - method to get the mean and standard deviations by base
get_raw_median() - method to get the raw median (median of All values in a base position for quality, not just individual means)
convert_signals_to_objects() - method to convert signals in df to signal objects
convert_signal_list() - method to convert a list of signals to a list of signal objects
convert_signal_obs_to_groups() - method to convert a df of signal objects to signal object groups
convert_signal_obj_list() - method to convert a list of signal object list to signalgroups
twas_the_monster_mash() - method to convert signal observations to usable features and squish it all together.
feature_extraction() - method to extract features to a pandas dataframe from a json file.
batch_read_alignment_data() - method to generate a list of subdivided batches of reads
extract_batch_features() - method to extract the features of a given batch.
batched_feature_extraction() - method to combine batching and extraction of features.
main() - method to run feature extraction on a given json file.

Classes:
Signal - a class that contains relevant information for a signal sequence
GroupedSignals - a grouping of Signal objects and relevant information and methods regarding that data

'''

import os
import sys
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import math
import json
import gc
import datetime
import multiprocessing
import argparse
import warnings
warnings.filterwarnings("ignore")

VERBOSE = False

def load_info_from_json(path):
    '''
    load_info_from_json loads the information from a consensus json file that
    was exported from merge_consensus.py's save_by_consensus(),
    and provides them as a handful of variables.
    
    Parameters:
    path: path, as str, to the consensus json file.
    
    Returns:
    consensus_id, data_dict, ref_seq, freq,
        consensus_id - the id of the consensus sequence
        data_dict - list of dicts of all the reads' data in the file
        ref_seq - the reference sequence as a string
        freq - sampling frequency of the reads in the file
    '''
    
    # create a dummy variable to hold the data
    data = None
    
    # parse the consensus id from the filename
    consensus_id = path.split('__')[-1].split(".")[0]
    
    # open the file
    with open(path) as f:
        # load the data from json
        data = json.load(f)
        
    # close the file
    f.close()
    
    # get the reference sequence from the header
    ref_seq = data[0]['ref_seq']
    
    # get the frequency from the header
    freq = data[0]['freq']
    
    # get the data 
    data_dict = data[1]
    
    #gc.collect()
    # return the information extracted.
    return consensus_id, data_dict, ref_seq, freq


def extract_coordinate_data(alignment_data):
    '''
    extract_coordinate_data() takes in the alignment data (a pandas series 
    containing a list of three datasets per entry) and then extracts the data to individual
    pandas series.
    
    Parameters:
    alignment_data: a pandas series containing lists of three lists of dicts where
    keys are bases and values are, in order, bases, qualities, and signals. Note: This could probably
    be done in a more storage-efficient format.
    
    Returns:
    a list containing three dataframes of consensus-aligned bases, qualities, and signals
    '''
    
    # convert the alignment data to a dataframe
    df = pd.DataFrame(alignment_data.to_list(), columns=['bases','qualities','signals'])
    
    # assign each dict in the dataframe to its own dataframe
    base_coordinate_df = pd.DataFrame(list(df['bases']))
    qual_coordinate_df = pd.DataFrame(list(df['qualities']))
    signal_coordinate_df = pd.DataFrame(list(df['signals']))
    
    # memory stuff
    del(df)
    
    # create a list to contain the dataframes
    data_out = [base_coordinate_df, qual_coordinate_df, signal_coordinate_df]
    
    # more garbage collection
    del(base_coordinate_df)
    del(qual_coordinate_df)
    del(signal_coordinate_df)
    #gc.collect()
    
    # return the dataframe list.
    return data_out


def align_observation_ops(cig, sig):
    '''
    align_observation_ops takes in a cigar and list of mapped signals, and then
    aligns the signal by cigar string operations.
    
    Parameters:
    cig: cigar information (a list of cigar touples)
    sig: a list of move-mapped signal data.
    
    Returns: a dataframe containing two series - one for operation, and 
             one for the signal corresponding to the operation. 
    '''
    
    # set up a dummy list to hold observations
    observation_dict_list = []
    
    # set up an indexer to index the operations by
    observation_indexer = 0

    # Loop through each observaiton in the cigar
    for observation in cig:
        
        # Get the observation type and duration
        ob_type = observation[0]
        ob_len = observation[1]
    
        # loop through the observation duration
        for i in range(ob_len):
            
            # insert the observation and the signal corresponding to that observation
            try:
                if ob_type == 2:
                    # The number of deletions correspond to the difference in signal count, so 
                    # a deletion means there is no signal at that observation (somehow)
                    observation_dict = {'operation':str(ob_type),
                                        'signal':None}
                    sig.insert(i+observation_indexer, None)
                    
                    # if it isnt a deletion, append the signal and operation
                else:
                    observation_dict = {'operation':str(ob_type),
                                        'signal':sig[i+observation_indexer]}
                    
            # if the above fails, it's probably because the signal list is empty, so append none.
            except:
                observation_dict = {'operation':str(ob_type),
                                    'signal':None}
                sig.insert(i+observation_indexer, None)

            # append the observation to the dummy list
            observation_dict_list.append(observation_dict)
            
        # index the observation indexer
        observation_indexer += ob_len
        #gc.collect()
        
    # convert to a dataframe, and filter the not-so-valuable operations
    obs_df = pd.DataFrame(observation_dict_list)
    valued_operations = obs_df[(obs_df['operation'] != '4') |( obs_df['operation'] != '4')].reset_index(drop='true')
    

    # return the valued operations
    return valued_operations


def shift_to_alignment(ops, sigs, seq, quals, ref_len, r):
    '''
    shift_to_alignment() shifts the alignment of signals, base operations, and qualities to
    their proper corresponding alignment for a signle read. 
    
    Parameters:
    ops: a pandas series of operations extracted from a cigar string on a per-read basis
    sigs: a pandas series of signal lists corresponding to the operations
    seq: the raw read sequence, as a list of chars.
    quals: the qualities of each observed base/operation, as a list.
    ref_len: the length of the reference sequence the read is mapped to, as an int
    r: the reference position the read is aligned to on the reference, as an int
    
    Returns: a list of three dictionaries of reference-aligned bases/operations,
             qualities, and signals. 
    
    '''
    # per-read alignment shifter
    
    # set up a dummy variable to hold the current base position
    base_position = 0
    
    # set up a list of the base keys
    base_keys = np.arange(ref_len).tolist()
    
    # set up empty dictionaries indexed by base position in the reference
    base_dict = {key:[] for key in base_keys}
    qual_dict =  {key:[] for key in base_keys}
    sig_dict = {key:[] for key in base_keys}
    
    # loop through each operation
    for i in range(len(ops)):
        op = ops[i]
        sig = sigs[i]
        
        shifted_bp = base_position + r
        # if it's a 0, we put the corresponding base in that index
        if (op == str(0)):
            # get the base at the base position in the sequence and append it,
            # likewise do the same for quality and signal
            base_dict[shifted_bp].append(seq[base_position])
            qual_dict[shifted_bp].append(quals[base_position])
            sig_dict[shifted_bp].append(sig)
            
            #index the base position by 1
            base_position += 1
        
        # if it's a 1, we have an insertion at that index
        elif (op == str(1)):
            
            # append the insertion indexer into the previous base position
            base_dict[shifted_bp-1].append("+")
            
            # append the quality to the previous base position
            qual_dict[shifted_bp-1].append(quals[base_position])
            
            # append the signal to the previous base position
            sig_dict[shifted_bp-1].append(sig)
            
            # delete the inserted base and quality to not get indexing errors
            del(seq[base_position])
            del(quals[base_position])
        
        # if it's a deletion,
        elif(op == str(2)):
            
            # insert a deletion indexer into the sequence and signals
            seq.insert(base_position, "-")
            quals.insert(base_position, None)
            
            # append the signal
            sig_dict[shifted_bp].append(sig)
            
            #append the deletion base
            base_dict[shifted_bp].append("-")
            
            # index the base position
            base_position += 1
            
            
    return [base_dict, qual_dict, sig_dict]


def parse_read_alignment_info(signal, cigar, ref_pos, qualities, read_sequence, reference_sequence):
    '''
    parse_read_alignment_info is used to prepare read data to be aligned to a reference,
    and then align it and provide the resulting information.
    
    Parameters:
    signal: an array of move-mapped signal lists.
    cigar: a list of touples of corresponding cigar information
    ref_pos: reference position on a read as an int
    qualities: an array of qualities
    read_sequence: a string containing the read sequence of bases
    reference_fasta: a string containing the reference sequence the read was mapped to.
    
    Returns: a list of three dictionaries corresponding to the
             aligned base operations, qualities, and signals.
    '''
    # per-row reader for the alignment info
    
    # Make copies to avoid issues with overwriting data, and ensure datatypes are right.
    sig = list(signal.copy())
    c = cigar.copy()
    r = int(ref_pos)
    q = qualities.copy()
    seq = list(read_sequence) 
    ref_seq_len = len(reference_sequence)
    
    # Align the signals to observations using the eponymous method w/
    # the cigar string, and the signals
    ops_and_signals = align_observation_ops(c, sig)
    
    # isolate the signals and operations:
    o = ops_and_signals['operation'].copy()
    op_sig = ops_and_signals['signal'].copy()
    
    # get the alignment dictionaries as a list of dicts
    dicts = shift_to_alignment(o, op_sig, seq, q, ref_seq_len, r)
    del sig
    del c
    del r
    del q
    del seq
    del ref_seq_len
    
    return dicts


def extract_alignment_from_dict(m_list_dict, r_seq):
    '''
    extract_alignment_from_dict takes in a list of dictionaries of merged information,
    and the reference sequence for the reads in the list of dictionaries, and 
    then generates the reference-based alignment data.
    
    Parameters:
    m_list_dict: a list of dictionaries of merged read data, all corresponding to a single reference sequence.
                 each dict must have keys [signal, cigar, ref, quality, sequence].
    r_seq: the reference sequence as a str that all reads in the dict correspond to.
    
    Returns:
    a list of dicts containing the data and newly added alignment data as a new key.
    '''
    # extract the alignment info from the merged dict
    e_df = pd.DataFrame(m_list_dict)
    
    # apply parse_read_alignment_info
    e_df['alignment_data'] = e_df.apply(lambda x: parse_read_alignment_info(x['signal'],
                                                                            x['cigar'],
                                                                            x['ref'],
                                                                            x['quality'],
                                                                            x['sequence'],
                                                                            r_seq), axis = 1)
    # convert back to a dict
    out_dict = e_df.to_dict("records")
    
    # delete the dataframe
    del(e_df)
    
    # collect the garbage.
    #gc.collect()
    
    # return the dict.
    return out_dict


def calc_base_probs(sequence_coordinates):
    """
    base_probs takes in a dataframe of base sequence coordinates and calculates the
    numerical probability of each base at each position by the number of occurences.

    Parameters:
    sequence_coordinates: A pandas dataframe where each row is a read, and each
                          column is a base position. Thus, each cell should have
                          a single base, as a string.

    Returns:
    a dataframe containing statistical likelyhoods of each base at a given position,
    with base positons as columns and base id's as rows. Has an extra column, n,
    that provides number of sample measurements in that row.
    """
    # convert the sequence coordinate lists to strings
    seq_coords = sequence_coordinates.applymap(lambda x: ''.join(x))
    
    # Create a list to hold each option for a base
    base_options = ['A', 'G', 'C', 'T', '-', '+']

    # Generate an empty list to hold the probabilities
    list_of_prob_dicts = []
    list_of_obs = []

    # for each column,
    for column in list(seq_coords.columns):
        
        # convert column to a list with each char being unique, positionally
        charlist = list(''.join(seq_coords[column].to_list()))
        
        # convert the charlist to a pandas series
        char_series = pd.Series((x[0] for x in charlist))

        n_obs = len(char_series)
        list_of_obs.append(n_obs)
        
        # get the number of each value in the column as a dictionary
        counts_dict = char_series.value_counts().to_dict()

        # Loop through each option in the base options
        for option in base_options:

            # if the option isn't in the counts dictionary, set the value to zero
            if option not in counts_dict:
                counts_dict[option]=0

        # Sum up the amount of the base
        total = sum(counts_dict.values(), 0.0)

        probs_dict = {}

        # Add each base to the probability dictionaries
        for i in range(len(counts_dict)):
            key = list(counts_dict.keys())[i]
            value = list(counts_dict.values())[i]
            try:
                probs_dict[key] = value/total
            except:
                probs_dict[key] = None
            
        # append the probability dictionary of the bases to the list
        list_of_prob_dicts.append(probs_dict)

    # convert the list of dicts to a dataframe
    return_df = pd.DataFrame(list_of_prob_dicts)
    return_df['n']=pd.Series(list_of_obs)
    
    #gc.collect()
    # Return the dataframe
    return return_df


def shannon_entropies(base_probs):
    """
    shannon_entropies takes in a dataframe with base options as columns and
    the base positions as rows, with the values being probability of each 
    base option at each position, and reuturns the shannon entropies
    of each position in the sequence.

    Parameters:
    base_probs: dataframe containing probabilities as described above

    Returns:
    a series of entropy at each position.
    """

    # Create a duplicate of the base probabilities so no overlapping
    # transforms. 
    tmp_df = base_probs.copy()

    # Set up a list for the base options
    base_options = ['A', 'G', 'C', 'T', '-', '+']

    n_col = tmp_df['n']
    tmp_df = tmp_df.drop(columns=['n'])
    tmp_df = tmp_df.fillna(value=np.nan)

    # loop through each base option in the probabilities
    for column in list(tmp_df.columns):

        # get the series of the column
        base_ser = tmp_df[column]

        # Compute the individual shannon entropy of the column
        tmp_df[column] = base_ser*np.log2(base_ser)

    # Sum all of the shannon entropies for each row, and *-1 to make them positive.
    entropy = tmp_df.sum(axis=1)*-1
    e_max = entropy.max()
    e_min = entropy.min()
    n_entropy = (entropy - e_min) / (e_max-e_min)

    #gc.collect()
    # Return the entropy
    return entropy


def get_first_base(base_list):
    '''
    get the first base in a list of bases.
    '''
    if base_list != []:
        return base_list[0]
    else:
        return ""
    

def check_base_match(base, refbase):
    '''
    check if a base matches the reference base.
    '''
    if base != "":
        if (refbase == base):
            return 0
        else:
            return 1
    else:
        return None
    
    
def mismatch_prob(sequence_coordinates, ref_seq):
    '''
    mismatch_prob checks the aligned reads versus a reference sequence,
    and returns the probability of a mismatch at each base position.
    '''
    
    # get the first base of each observation of each read.
    seq_coords = sequence_coordinates.applymap(lambda x: get_first_base(x))
    
    # loop through each base in the reference sequence
    for i in range(len(ref_seq)):
        
        # get the corresponding column of read bases at that index
        refbase = ref_seq[i]
        
        # check if the bases match
        seq_coords[i] = seq_coords[i].apply(lambda x: check_base_match(x, refbase))
        
    # get the mean value of each, which is also the probability in this instance.
    means = seq_coords.mean()
    
    # return the probabilities.
    return means


def is_in_del(ls_bases):
    
    '''
    is_in_del() takes in a list of bases,
    and checks if the base is an insertion or deletion, 
    and returns a list of booleans corresponding to if the 
    base is either an insertion or deletion.
    
    Parameters:
    ls_bases: a list of bases or operations.
    
    Returns:
    list of booleans corresponding to if the base is an insertion/deletion or not.
    '''
    
    # dummy list to hold booleans
    ls_bools = []
    
    # loop through each base
    for base in ls_bases:
        
        # check if it's an insertion, deletion, or other
        if base == "+":
            ls_bools.append(True)
        elif base == "-":
            ls_bools.append(True)
        else:
            ls_bools.append(False)
            
    # return the list of booleans
    return ls_bools


def filter_by_indel(booleans, values):
    '''
    filter_by_indel() takes in a list of booleans and values (of the same length)
    and returns the values where there was not an insertion or deletion.
    
    Parameters:
    booleans: a list of bools corresponding to i/d or not.
    values: a list of values corresponding to each index in the booleans
    
    Returns: 
    a filtered list to contain only the values at positions that werent insertions or deletions.
    '''
    # generate dummy list for output
    out = []
    
    # loop through the booleans in the list
    for i in range(len(booleans)):
        
        # check if the boolean is not true
        if not booleans[i]:
            
            # append the value at that index
            out.append(values[i])
            
    return out


def remove_in_dels(base_coords, to_filter_coords):
    '''
    remove_in_dels takes in base coordinate dataframe, and a base coordinate dataframe of 
    data to be filtered and then filters that data by if it had insertions/deletions or not
    
    Parameters:
    base_coords: a base coordinate dataframe (columns are index of alignment sequence, rows are reads)
    to_filter_coords: a coordinate dataframe of unknown values to be filtered
    
    Returns:
    a coordinate dataframe with all values corresponding to insertions or deletions removed
    '''
    
    # apply the is_in_del function to each value in the df
    boolmap = base_coords.applymap(lambda x: is_in_del(x))
    
    # copy the filter df
    filter_df = to_filter_coords.copy()
    
    # loop through each column in the dataframe that will be filtered
    for column in to_filter_coords.columns:
        
        # get the transpose of the dataframe that is the combination of the current column of the
        # boolmap and filter_coords
        joined_map = pd.DataFrame([boolmap[column], to_filter_coords[column]]).T
        
        # apply the filteR_by_indel method to each row in this dataframe,
        # and set it to the filter_df's column
        joined_map = joined_map.set_axis(['m', 'q'], axis=1)
        filter_df[column] = joined_map.apply(lambda x: filter_by_indel(x['m'], x['q']), axis=1)
        
    # return the filter df.
    return filter_df
        
        
def getmean(lsval):
    '''
    get the mean value from a list of values,
    while checking if it isnt empty, either.
    '''
    if (lsval != []):
        mean = np.mean(np.asarray(lsval))
        return mean
    else:
        return None
    
    
def mean_std_qual_by_base(qual_coord):
    '''
    mean_std_qual_by_base() gets the mean and standard deviation quality value by
    the base position.
    
    Parameters:
    qual_coord: coordinate dataframe of qualities of each base in each read mapped to a reference
    
    Returns:
    (meanvals, stdvals, medians) - the means, standard deviations, and median qualities
    as pandas series for each base in the reference sequence
    '''
    
    # apply the getmean method to each cell
    qual_means = qual_coord.applymap(getmean)
    
    # get the mean of each column
    meanvals = qual_means.mean()
    
    # get the std of each column
    stdvals = qual_means.std()
    
    # get the median of each column
    medians = qual_means.median()
            
    # return the mean, std, medians
    return meanvals, stdvals, medians


def get_raw_median(qual_cord):
    '''
    get_raw_median takes in a quality coordinate dataframe,
    and then provides the median quality of each base in the reference.
    
    Parameters:
    qual_cord: quality coordinate dataframe
    
    Returns:
    pandas series of median quality where ALL qualities are involved, not averaging
    qual per base.
    '''
    # empty list for medians
    medians = []
    
    for column in qual_cord.columns:
        # get the series of qualities for the base
        base_qual_ser = qual_cord[column]
        
        # generate an empty list to hold the qualities
        base_qual_list = []
        for i in range(len(base_qual_ser)):
            
            # get the list of qualities for the current base for a specific read
            read_quals = base_qual_ser[i]
            
            for j in range(len(read_quals)):
                # append each observed base for the current read to the list of quals
                base_qual_list.append(read_quals[j])
        # convert to numpy array
        bq_array = np.asarray(base_qual_list)
        med = np.median(bq_array)
        medians.append(med)
        
    return pd.Series(medians)
        
                
def convert_signals_to_objects(seq_coords, f):
    '''
    convert_signals_to_obects converts each signal value in the sequence coordinates
    to a signal object using the convert_signal_list method.
    
    Parameters:
    seq_coords: read signal coordinate dataframe
    f: sampling frequency as an int. 
    '''
    
    converted = seq_coords.applymap(lambda x: convert_signal_list(x, f))
    return converted


# signal features
class Signal:
    '''
    Signal class is designed to contain and automatically run a handful of valuable pieces
    of information about a signal list, such as mean, std, etc.
    
    
    Parameters:
    signal: list of signal values in picoamperes
    freq: sampling frequency
    
    Methods:
        None
    '''
    def __init__(self, signal, freq):

        # get the mean, std, median, duration, signal, coulombic signal,
        # coulombs moved, mean coulombs, std coulombs, median coulombs
        signal_array = np.asarray(signal)
        self.mean = np.mean(signal_array)
        self.std = np.std(signal_array)
        self.median = np.median(signal_array)
        time_step = 1/freq
        time = np.full((signal_array.shape)[0], time_step)
        self.time = np.sum(time)
        self.signal = signal_array
        self.coulomb_signal = signal_array*time_step
        self.coulombs_moved = self.coulomb_signal.sum()
        self.mean_coulombs_per_step = np.mean(self.coulomb_signal) # per time step
        self.std_coulombs_per_step = np.std(self.coulomb_signal) # per time step
        self.median_coulombs_per_step = np.median(self.coulomb_signal)
        self.fft_peak = (np.abs(np.fft.fft(signal_array))**2).max()
    
    
class GroupedSignals:
    '''
    GroupedSignals is designed to contain useful methods for processing groups of signal objects.
    '''
    def __init__(self, sig_list):
        '''
        __init__ takes in a list of signal objects and loads the data to relevant variables
        '''
        
        if (type(sig_list) != type(None)):
            self.signal_list = sig_list

            # generate the features of the group.
            self._get_group_features()
            
        else:
            self.signal_list = None
            self.mean_of_means = None
            self.pooled_std = None
            self.total_signal = None
            self.total_signal_mean = None
            self.total_signal_std = None
            self.total_signal_median = None
            
        
        
    def _get_group_features(self):
        '''
        private method for calculating the group features.
        '''
        
        # empty lists for values
        sigmeans = []
        sigstds = []
        total_signal = []
        total_time = 0
        total_coulombs = 0 
        fft_peaks = []
        # loop through each signal in the list
        for signal in self.signal_list:
            total_time += signal.time
            total_coulombs += signal.coulombs_moved
            # get the mean, std, and total signal
            sigmeans.append(signal.mean)
            sigstds.append(np.square(signal.std))
            total_signal.extend(list(signal.signal))
            fft_peaks.append(signal.fft_peak)
            
        # get the mean of means
        mean_of_means = np.mean(np.asarray(sigmeans))
        
        # get pooled standard deviation
        pooled_std = np.sqrt(np.asarray(sigstds).sum()/(len(sigstds)))
        
        # get the total signal and get the mean, std, and median of the total signal
        total_signal = np.asarray(total_signal)
        total_sig_mean = np.mean(total_signal)
        total_sig_std = np.std(total_signal)
        total_sig_median = np.median(total_signal)
        fft_peaks = np.asarray(fft_peaks)
        mean_fft_peak = np.mean(np.asarray(fft_peaks))
        
        # assign to class variables -- should these be calculated on demand rather than stored?
        # yes, but that takes more time we don't have rn.
        self.mean_of_means = mean_of_means
        self.pooled_std = pooled_std
        self.total_signal = total_signal
        self.total_signal_mean = total_sig_mean
        self.total_signal_std = total_sig_std
        self.total_signal_median = total_sig_median
        self.total_fft_peak = (np.abs(np.fft.fft(total_signal))**2).max()
        self.mean_fft_peaks = mean_fft_peak
        self.residence_time = total_time
        
        
        
        
        #### MORE GO BELOW, DONT KNWO WHAT YET
        
        return None

    
def convert_signal_list(sig_list, f):
    '''
    convert_signal_list takes in a list of signals and a sampling frequency,
    and returns a list of signal objects.
    
    Parameters:
    sig_list: list of lists of signals
    f: sampling frequency, hz
    
    returns:
    list of signal objects,
    OR
    None if there are no signals in the list.
    '''
    # check if the list is not empty
    if (sig_list != []):
        return_list = []
        
        # loop through the signals
        for i in range(len(sig_list)):
            if type(sig_list[i]) != type(None):
                # create a signal object with that subsect of the signal and that frequency
                sig_ob = Signal(sig_list[i], f)
                return_list.append(sig_ob)
        return return_list
    else:
        return None
    

def convert_signal_obs_to_groups(sig_obj_coords):
    '''
    map the conversion of signal object lists to signal groups.
    '''
    converted = sig_obj_coords.applymap(lambda x: convert_signal_obj_list(x))
    return converted


def convert_signal_obj_list(x):
    if type(x) != type(None) and (len(x) > 0):
        gp = GroupedSignals(x)
        return gp
    else:
        return None
    

def twas_the_monster_mash(grouped_coords):
    # it was a graveyard bash
    
    def get_info_if_not_none(val):
        if type(val) != type(None):
            return [val.total_signal_mean, val.total_signal_std, val.total_signal_median, val.total_fft_peak, val.residence_time, val.mean_fft_peaks]
        else:
            return [None, None, None, None, None, None]
    base_dict = {}
    for column in grouped_coords:
        base_col = grouped_coords[column]
        grouped_vals = base_col.apply(lambda x: get_info_if_not_none(x))
        
        sep_vals = pd.DataFrame(grouped_vals.to_list(), columns=['sig_means', 'sig_stds', 'sig_meds', 'sig_peaks', 'sig_time', 'sig_peaks_means'])

        base_features = {}
        mean_vals = dict(sep_vals.mean())
        for key in mean_vals.keys():
            base_features['mean_{}'.format(key)] = mean_vals[key]
        
        std_vals = dict(sep_vals.std())
        for key in mean_vals.keys():
            base_features['std_{}'.format(key)] = std_vals[key]
            
        med_vals = dict(sep_vals.median())
        for key in mean_vals.keys():
            base_features['median_{}'.format(key)] = med_vals[key]
            
                         
        base_dict[column] = base_features
        
    return pd.DataFrame(base_dict)

    

def feature_extraction(json_path):
    '''
    feature_extraction takes a path to a json file, and loads all of the features in that dataset
    to a dataframe, with a column called 'XNA_PRESENT' representing the idx where the xna is present
    
    Parameters:
    json_path: path, as str, to a json file that is exported by merge_consensus.
    
    Returns:
    pandas Dataframe with all features and corresponding values.
    '''
    
    if VERBOSE: print(str(datetime.datetime.now()) + ' loading json... ')
    # yoink the consensus id, data dictionary, reference sequence, and frequency from the json path
    consensus_id, data_dict, ref_seq, freq = load_info_from_json(json_path)
    if VERBOSE: print(str(datetime.datetime.now()) + ' json loaded. ')
    
    if VERBOSE: print(str(datetime.datetime.now()) + ' reading alignment... ')
    # get the read_alignment_data
    read_alignment_data = extract_alignment_from_dict(data_dict, ref_seq)
    
    # get te read_alignment dataframe
    read_alignment_df = pd.DataFrame(read_alignment_data)
    if VERBOSE: print(str(datetime.datetime.now()) + ' alignment read. ')
    
    ## IF THERE IS A SPLIT IN DATA, SPLIT IT HERE AND MAKE THE FOLLOWING A PER-SPLIT METHOD
    
    if VERBOSE: print(str(datetime.datetime.now()) + ' generating data coordinates... ')
    # base-based features, base_probs is 6 but lets just look at insertions and deletions
    # get the alignment coordinate data
    data_coordinate_list = extract_coordinate_data(read_alignment_df['alignment_data'])

    if VERBOSE: print(str(datetime.datetime.now()) + ' data coordinates generated. ')
    if VERBOSE: print(str(datetime.datetime.now()) + ' calculating observation features... ')
    # get the base probabilites at each position -- THIS HAS INSERTION AND DELETION PERCENTAGE BUILT IN
    base_probs = calc_base_probs(data_coordinate_list[0])

    # get the probability of mismatch at each position
    mismatch_probs = mismatch_prob(data_coordinate_list[0], ref_seq)

    # Get the shannon entropy of each base position.
    shentropy_series = shannon_entropies(base_probs)
    
    if VERBOSE: print(str(datetime.datetime.now()) + ' observation features calculated.  ')
    if VERBOSE: print(str(datetime.datetime.now()) + ' calculating quality features...  ')
    # quality score features
    # get the qualities with no insertions or deletions
    no_in_del_qual = remove_in_dels(data_coordinate_list[0], data_coordinate_list[1])

    # get the mean, std, and median with no insertions or deletions
    no_mvals, no_stvals, no_meanmedian = mean_std_qual_by_base(no_in_del_qual)

    # get the mean, std, and median while including insertions or deletions
    mvals, stvals, meanmedian = mean_std_qual_by_base(data_coordinate_list[1])

    # get the raw median quality score of each base.
    rawmeds = get_raw_median(data_coordinate_list[1])
    
    if VERBOSE: print(str(datetime.datetime.now()) + ' quality features calculated.  ')
    if VERBOSE: print(str(datetime.datetime.now()) + ' calculating signal features...  ')
    # Signal features
    # convert the signals to objects
    sig_obj_coords = convert_signals_to_objects(data_coordinate_list[2], freq)

    # convert the signal lists to signal group objects
    sig_group_coords = convert_signal_obs_to_groups(sig_obj_coords)
    
    #get the signal objects with no insertions or deletions
    no_sig_obj_coords = convert_signals_to_objects(remove_in_dels(data_coordinate_list[0], data_coordinate_list[2]), freq)

    # get signal groups with no insertions or deletions
    no_sig_group_coords = convert_signal_obs_to_groups(no_sig_obj_coords)
    
    # get the mess of signal features
    mess_with_in_del = twas_the_monster_mash(sig_group_coords)
    mess_wo_in_del = twas_the_monster_mash(no_sig_group_coords)
    
    if VERBOSE: print(str(datetime.datetime.now()) + ' signal features calculated.  ')
    if VERBOSE: print(str(datetime.datetime.now()) + ' assembling features...  ')
    # combine all of the features as well as the xna position into one cumulative dataframe
    mismatch_probs = pd.DataFrame(mismatch_probs, columns=['mm_prob'])
    shentropy_col = pd.DataFrame(shentropy_series, columns=['shentropy'])
    no_quals = pd.concat([no_mvals, no_stvals, no_meanmedian], axis=1)
    no_quals.columns = ['n_qmean', 'n_qst', 'n_qmed']
    id_quals = pd.concat([mvals, stvals, meanmedian], axis=1)
    id_quals.columns = ['i-d_qmean', 'i-d_qst', 'i-d_qmed']
    rawmeds = pd.DataFrame(rawmeds, columns = ['r_qmed'])
    in_del_sigs = mess_with_in_del.T.add_suffix('_i-d')
    wo_in_del_sigs = mess_wo_in_del.T.add_suffix('_w/o')

    xna_idx = int(consensus_id.split(':')[-1].split(']')[0])
    xna_df = pd.DataFrame(mismatch_probs.index, columns = ['XNA_PRESENT'])
    xna_df['XNA_PRESENT'] = 0
    xna_df['XNA_PRESENT'][xna_idx] = 1
    
    assembled_features = pd.concat([xna_df, base_probs, mismatch_probs, shentropy_col, no_quals, id_quals, rawmeds, in_del_sigs, wo_in_del_sigs], axis=1)
    if VERBOSE: print(str(datetime.datetime.now()) + ' features assembled.  ')
    
    return assembled_features
    
    
def batch_read_alignment_data(read_alignment_data_list, batch_size):
    '''
    batch_read_alignment_data takes in a list of read alignment data, and a batch size,
    and generates a list of subdivided batches of the size of the batch. if there arent enough to fill the last,
    it populates it with what it can - occasionally resulting in smaller batches.
    
    Parameters:
    read_alignment_data_list: a list of dicts containing info as exported by extract_alignment_from_dict.
    batch_size: a positive integer representing the size of each batch, in reads.
    
    Returns:
    a list of batched reads with the specified batchsize.
    '''

    # generate the list
    batched_reads = [read_alignment_data_list[i:i+batch_size] for i in range(0, len(read_alignment_data_list), batch_size)]
    
    # return the batched reads
    return batched_reads
    

def extract_batch_features(read_batch, ref_seq, freq, consensus_id):
    '''
    extract_batch_features takes in a read batch, and the reference data,
    and then generates the features for that batch as if it were an entire consensus.
    
    Parameters:
    read_batch: a list of dicts of read data
    ref_seq: a reference sequence corresponding to the reads, as str
    freq: a frequency value, as float or int, corresponding to read sampling frequency
    consensus_id: the id of the consensus sequence, containing xna location.
    
    Returns:
    a conjunction of the features of the batch, as a pandas dataframe.
    '''
    read_alignment_df = pd.DataFrame(read_batch)
    data_coordinate_list = extract_coordinate_data(read_alignment_df['alignment_data'])
    # get the base probabilites at each position -- THIS HAS INSERTION AND DELETION PERCENTAGE BUILT IN
    base_probs = calc_base_probs(data_coordinate_list[0])

    # get the probability of mismatch at each position
    mismatch_probs = mismatch_prob(data_coordinate_list[0], ref_seq)

    # Get the shannon entropy of each base position.
    shentropy_series = shannon_entropies(base_probs)
    # quality score features
    # get the qualities with no insertions or deletions
    no_in_del_qual = remove_in_dels(data_coordinate_list[0], data_coordinate_list[1])

    # get the mean, std, and median with no insertions or deletions
    no_mvals, no_stvals, no_meanmedian = mean_std_qual_by_base(no_in_del_qual)

    # get the mean, std, and median while including insertions or deletions
    mvals, stvals, meanmedian = mean_std_qual_by_base(data_coordinate_list[1])

    # get the raw median quality score of each base.
    rawmeds = get_raw_median(data_coordinate_list[1])
    
    # Signal features
    # convert the signals to objects
    sig_obj_coords = convert_signals_to_objects(data_coordinate_list[2], freq)

    # convert the signal lists to signal group objects
    sig_group_coords = convert_signal_obs_to_groups(sig_obj_coords)
    
    #get the signal objects with no insertions or deletions
    no_sig_obj_coords = convert_signals_to_objects(remove_in_dels(data_coordinate_list[0], data_coordinate_list[2]), freq)

    # get signal groups with no insertions or deletions
    no_sig_group_coords = convert_signal_obs_to_groups(no_sig_obj_coords)
    
    # get the mess of signal features
    mess_with_in_del = twas_the_monster_mash(sig_group_coords)
    mess_wo_in_del = twas_the_monster_mash(no_sig_group_coords)
    
     # combine all of the features as well as the xna position into one cumulative dataframe
    mismatch_probs = pd.DataFrame(mismatch_probs, columns=['mm_prob'])
    shentropy_col = pd.DataFrame(shentropy_series, columns=['shentropy'])
    no_quals = pd.concat([no_mvals, no_stvals, no_meanmedian], axis=1)
    no_quals.columns = ['n_qmean', 'n_qst', 'n_qmed']
    id_quals = pd.concat([mvals, stvals, meanmedian], axis=1)
    id_quals.columns = ['i-d_qmean', 'i-d_qst', 'i-d_qmed']
    rawmeds = pd.DataFrame(rawmeds, columns = ['r_qmed'])
    in_del_sigs = mess_with_in_del.T.add_suffix('_i-d')
    wo_in_del_sigs = mess_wo_in_del.T.add_suffix('_w/o')

    xna_idx = int(consensus_id.split(':')[-1].split(']')[0])
    xna_df = pd.DataFrame(mismatch_probs.index, columns = ['XNA_PRESENT'])
    xna_df['XNA_PRESENT'] = 0
    xna_df['XNA_PRESENT'][xna_idx] = 1
    
    assembled_features = pd.concat([xna_df, base_probs, mismatch_probs, shentropy_col, no_quals, id_quals, rawmeds, in_del_sigs, wo_in_del_sigs], axis=1)
    
    return assembled_features


def batched_feature_extraction(json_path, batch_size):
    '''
    batched_feature_extraction takes in a json filepath and a batch size, and returns
    a list of the assembled features as batches.
    
    Parameters:
    json_path: the path, as str, to a json file of data
    batch_size: positive int representing size of each batch of reads
    '''
    
    # extract the cnsensus id, data dict, reference sequence, and frequency from the json
    consensus_id, data_dict, ref_seq, freq = load_info_from_json(json_path)
    
    # read the alignment data
    read_alignment_data =extract_alignment_from_dict(data_dict, ref_seq)
    
    # convert the reads to batches of reads
    batched_reads = batch_read_alignment_data(read_alignment_data, batch_size)
    
    # create an empty list to hold the touples of reads
    prepped_batches = []
    
    # loop through each batch of reads and convert it to a touple with the reference info
    for i in range(len(batched_reads)):
        batch = batched_reads[i]
        prepped_batches.append((batch, ref_seq, freq, consensus_id))
        
    # create a list of the assembled features
    batched_assembled_features = []
    
    # use the multiprocessing to run all the batches at once. 
    with multiprocessing.Pool(processes=len(batched_reads)) as feat_gen_pool:
        batched_assembled_features = feat_gen_pool.starmap(extract_batch_features, prepped_batches)
        
    # return the batched features.
    return batched_assembled_features

    
def main(json_path, batch_size=None):
    '''
    main takes in a json path and batch size, and returns a list of assembled features by batch.
    
    Parameters:
    json_path: path, as str, to json file
    batch_size: default None, size of batch of reads.
    '''
    
    if VERBOSE: print('[ feature_extraction.py ] Batch Size: {}'.format(batch_size))
    # create an empty list to populate with the feature batches
    assembled_features = []
    
    # check batch size, get the features accordingly
    if type(batch_size) == type(None):
        assembled_features = [feature_extraction(json_path)]
    else:
        assembled_features = batched_feature_extraction(json_path, int(batch_size))
    
    return assembled_features


if __name__ == '__main__':
    
    parser = argparse.ArgumentParser(description='Convert json consensus files to feature information.')
    parser.add_argument('-v',"--verbose", help="Make output verbose.", action="store_true")
    parser.add_argument('-batch_size', help="Size of batches of reads.")
    parser.add_argument('-json_path', help="Json filepath.")
    parser.add_argument('-output', help="Output filepath.")
    
        # Parse the arguments
    args = parser.parse_args()
    
    # check the arguments and assign them
    VERBOSE = args.verbose
    batch_size = args.batch_size
    json_path = args.json_path
    output_path = args.output

    if type(json_path) == type(None): bam_path = input("[ feature_extraction.py ] Json Path: ")
    if type(output_path) == type(None): pod5_path = input("[ feature_extraction.py ] Output path: ")
    
    assembled_features = main(json_path, batch_size)
    
    for i in range(len(assembled_features)):
        batched_consensus = assembled_features[i]
        output_filepath = output_path + "consensus{}.parquet".format(i)
        print(output_filepath)
        batched_consensus.to_parquet(output_filepath)
        
    sys.exit()
    
        
    