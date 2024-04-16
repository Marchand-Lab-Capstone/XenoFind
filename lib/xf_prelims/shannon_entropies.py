"""
shannon_entropies.py
J. Sumabat, N. Lai, S. Peck, Y. Huang
4/10/2024 -- Updated by S. Peck

shannon_entropies.py contains methods involved with generating
'shannon entropy' values for consensus data and the corresponding reads.

consensus_formatter() - generates dataframe of ids and
                        sequences of all valid samples in 
                        a consensus fasta file.
get_record_indexes() - gets the indexes of the 'zone of interest'
                       that are in the consensus fasta
load_in_data() - converts information in a bamfile into a dataframe
shift_sequence() - uses a cigar string to get bases and operations by base index
format_sequences() - uses the shift_sequence method to format a dataframe of reads
calc_base_probs() - gets the probability (n instances/m options) of
                    each base or operation per read.
shannon_entropies() - uses base probabilities to generate shannon entropies.
get_top_n_reads() - get the top n alignments by number of reads aligned
wrapper() - generate entropy dataframe from consensus fasta and from a bamfile
generate_frw_figure() - produces a figure of the forward, reverse,
                        and all alignment entropies
"""

# Import relevant libraries
import pysam
import os
from Bio import SeqIO
import pandas as pd
import tables
import h5py # not sure if nessecary, too afraid to check
import re
import matplotlib.pyplot as plt
import numpy as np
import scipy.signal as sig

def consensus_formatter(consensus_fasta_path):
    """
    consensus_formatter takes in a consensus fasta path, 
    and returns a dataframe containing the ids and sequences of all valid samples.

    Parameters:
    consensus_fasta_path: a path, as str, to the consensus fasta
    
    Returns:
    a pandas dataframe containing id and sequence.
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

    # return the dataframe
    return df


def get_record_indexes(consensus_path):
    """
    get_record_indexes takes in a path to a labeled consensus file, and returns a list
    containing dictionaries of the start and end indicies of the unknown region of the consensus.

    Parameters:
    consensus_path: path to consensus fasta (typically labeled_consensus.fasta) as a string

    Returns:
    list of dicts of indexes.
    
    """

    # Generate an empty list to hold the indexes 
    list_record_indexes = []

    # Loop through each record in the fasta file
    for record in SeqIO.parse(consensus_path, "fasta"):

        # For each record, get the header
        header = record.description

        # Get the boundary conditions by searching
        match = re.search(r'BC 1: (\d+), BC 2: (\d+)', header)

        # set the boundary conditions to invalid values before adding 'em
        idxs = {"start":-1,"end":-1}

        # If the search has a match, add the indexes to the record
        if match:
            idxs["start"] = int(match.group(1))
            idxs["end"] = int(match.group(2))
        list_record_indexes.append(idxs)

    # return the list of the record indexes
    return list_record_indexes


def load_in_data(bamfile_path):
    """
    load_in_data takes in a path to a basecalled and aligned consensus bamfile,
    and converts the data into a dataframe format.

    Parameters:
    bam_path: path, as str, to the basecalled and aligned consensus bamfile.

    Returns:
    pandas dataframe of the relevant data. [sequence id, sequence, quality, length, reference index]
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
                read_dict = {'seq_id':sequence_id,
                             'ref_name':ref_name,
                             'sequence':seq,
                             'quality':quality,
                             'cigar':cigar,
                             'len':len(seq),
                             'ref':ref,
                             'rev':rev,
                             'moves':tag_dict['mv'],
                             'sig_len': tag_dict['ns'],
                             'trim_ofs':tag_dict['ts']}
                #if (ref == 0):

                # Append the dictionary to the data list
                data_list.append(read_dict)

    # convert the list of dicts into a dataframe
    df = pd.DataFrame(data_list)

    # convert the qualities to list.
    df['quality'] = df['quality'].apply(list)

    # return the dataframe
    return df


def shift_sequence(sequence, ref, cigar):
    """
    shift_sequences takes in a sequence, corresponding quality, and reference position, and
    shifts it to align it with the proper positioning.

    Parameters:
    sequence: DNA read sequence as str.
    quality: DNA read qualities as list, same length of str.
    ref: reference position on the reference read, as int.
    cigar: the cigar string of the read, as str.

    returns:
    a dict of base positions with corresponding bases or insertions/deletions.
    """

    # Get the sequence as a list
    seq = list(sequence)

    # set current index to zero,
    # make a list to hold all of the operations.
    current_index = 0
    operations = []

    # set up an empty dict for the true sequence and qualities
    true_seq = {}
    true_quals = {}

    # for each index of an offset of a reference,
    # append 'S' to the operation sequence
    for i in range(ref):
        operations.append("S")

    # For each operation in the cigar,
    for operation in cigar:

        # Get the operation type and length
        op_type = operation[0]
        op_length = operation[1]

        # Loop through the operations
        for i in range(op_length):

            # If the operation isn't four or 5, 
            if (op_type != 4 and op_type != 5):

                # Append the operation to the list.
                operations.append(str(op_type))

    #print(''.join(seq))
    #print(''.join(operations))
    base_indexer = 0
    dels = 0
    for i in range(len(operations)-1):

        if (operations[i] == str(0)):
            #print(''.join(true_seq.values()))
            #print("Base Index: {}, Operation index: {}, Operation: {}".format(base_indexer, i, operations[i]))
            if (base_indexer in true_seq):
                true_seq[base_indexer] += seq[base_indexer]
            else:
                true_seq[base_indexer] = seq[base_indexer]
            base_indexer += 1
            
        elif (operations[i] == str(1)):
            if (base_indexer in true_seq):
                true_seq[base_indexer] += "+"
            else:
                true_seq[base_indexer] = "+"
                
        elif (operations[i] == str(2)):
            seq.insert(base_indexer, "-")
            if (base_indexer in true_seq):
                true_seq[base_indexer] += "-"
            else:
                true_seq[base_indexer] = "-"
            base_indexer += 1

        elif (operations[i] == "S"):
            seq.insert(base_indexer, "N")
            if (base_indexer in true_seq):
                true_seq[base_indexer] += ""
            else:
                true_seq[base_indexer] = ""
            base_indexer += 1
        
                
        # Delete the inserted base to shift everything over
        #del seq[index]
    
        # 3 - skipped region from reference
        # 4 - soft clipping
        # 5 - hard clipping
        # 6 - padding
        # 7 - sequence match
        # 8 - sequence mismatch

    # Return the sequence and the qualities
    return true_seq


def format_sequences(reads_df, indexes):
    """
    format_sequences takes in a dataframe of data previously exported from a sam file, and 
    indexes of the start and end of an unknown region, and returns the reads and quality scores for those segments.

    Parameters:
    reads_df: dataframe containing [sequence id, sequence, quality, length, reference index]
    indexes: a dictionary {"start", "end"} that contains the start and end indexes of an unknown region.

    Returns:
    a dataframe with rows indexed by read and columns by base position, containing strings
    of the corresponding bases at that location. 
    """

    #step 1: remedying insertions, deletions and alignment

    # Duplicate the reads dataframe
    l_df = reads_df.copy()

    # Shift each read using the shift_sequence method
    l_df['shift_dict'] = l_df.apply(lambda x: shift_sequence(x['sequence'],
                                                        x['ref'],
                                                        x['cigar']),
                               axis = 1)

    # Create new columns and split the shift touple column to two individual columns

    shift_bases_df = pd.DataFrame(l_df['shift_dict'].to_list())

    
    #Step 2, trimming off the ends

    # Trim the sequences and qualities
    if (indexes["start"] == None):
        indexes["start"] = 0
    if (indexes["end"] == None):
        indexes["end"] = len(shift_bases_df.columns)
    shift_bases_df = shift_bases_df.iloc[:, indexes["start"]:indexes["end"]]
    shift_bases_df.fillna('', inplace=True)
    

    # return the dataframe
    return shift_bases_df


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

    # Create a list to hold each option for a base
    base_options = ['A', 'G', 'C', 'T', '-', '+']

    # Generate an empty list to hold the probabilities
    list_of_prob_dicts = []
    list_of_obs = []

    # for each column,
    for column in list(sequence_coordinates.columns):
        
        # convert column to a list with each char being unique, positionally
        charlist = list(''.join(sequence_coordinates[column].to_list()))
        
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

    # Return the entropy
    return entropy


def get_top_n_reads(cons_ids, read_ids, n=10):
    """
    get_top_n_reads takes in a pandas series of alignment id's, and a pandas series of
    id's that reads aligned to (from a loaded bam or sam file),
    and returns a list of the ids of the top n alignments with the most corresponding
    reads.
    
    Parameters:
    cons_ids: pandas series of consensus sequence id's
    read_ids: pandas series of the ids a read aligned to
    n: int value (MUST BE LOWER THAN CONS_IDS LENGTH) of the top 
        number of reads desired. Default n = 10
        
    Returns:
    a list of the top n id's sorted by number of aligned reads
    """
    
    # Get the number of instances of each alignment id in the reads
    count_series = read_ids.value_counts()
    
    # filter to make sure the reads in counts correspond to the consensus id's
    counts_in_alignment = count_series[count_series.index.isin(cons_ids)]
    
    # Get the top n id's
    n_counts = counts_in_alignment[:n].index
    
    # return as a list
    return list(n_counts)


def wrapper(c_f, s_p, rev = None, n=10, verbose = False):
    """
    wrapper takes in the path to a consensus fasta,
    the path to a bamfile, 
    and returns a dataframe of the shannon entropies of 
    the top n alignments by read.
    
    Parameters:
    c_f: consensus fasta path, as str.
    s_p: bam file path, as str.
    Optional: 
    rev: (default) None - return all entropies
         True - return entropies of reads that were reversed
         False - return entropies of reads that were not reversed
    n: int value lower than the number of alignments, representing the top n
       alignments by read
    debug: boolean on wether or not to show steps
       
    Returns:
    pandas dataframe containing shannon entropies 'sh_ent',
    and alignment id 'alignment_id'.
    """
    
    # use the consensus formatter to format the consensus data
    form_consens = consensus_formatter(c_f)
    
    # Load in the bamfile data 
    loaded_df = load_in_data(s_p)
    
    # Get the top n alignments by read
    top_read_ids = get_top_n_reads(form_consens['id'], loaded_df['ref_name'], n=n)
    
    # Generate adummy list to hold the entropy data
    ent_data = []
    
    # loop through the top reads, 
    for i in range(len(top_read_ids)):
        
        # if verbose, print progress
        if verbose:
            print("{}/{} completed".format(i,len(top_read_ids)), end="\r")
        
        # filter the loaded data to data with just the id of the current top read
        subset_df = loaded_df[loaded_df['ref_name'] == top_read_ids[i]]
        
        # if the length of the top read is greater than zero:
        if (len(subset_df) >0):
            
            # Get the read indexes from the consensus fasta
            read_indexes = {'start':form_consens['start'][i], 'end':form_consens['end'][i]}
            
            # Get the series of booleans for if reversed
            reindexed_reverse = pd.Series(list(subset_df['rev'])) 
            
            # format the dna sequence 
            formatted_df = format_sequences(subset_df, read_indexes)
            
            # Check the reverse parameter
            if (rev is not None):
                # filter to either reversed or not
                formatted_df = formatted_df[reindexed_reverse==rev]
            
            # get the base probabilities of the formatted data
            base_probs_df = calc_base_probs(formatted_df)

            # Get the shannon entropies from the probabilities
            shentropy_series = shannon_entropies(base_probs_df)

            # put the id and entropy into a dict
            ent_dict = {'alignment_id':top_read_ids[i],
                        'sh_ent':list(shentropy_series)}
            
            # append the dict to the list
            ent_data.append(ent_dict)

    # convert to a dataframe
    ent_df = pd.DataFrame(ent_data)
    
    # return the dataframe of id's and entropies
    return ent_df


def generate_frw_figure(bam_path_list, cons_path_list, out_path, n=10, verbose=False):
    """
    Generate a figure that has the forward, reverse, and both shannon entropies
    for the top n reads to a specified path.
    
    Parameters:
    bam_path_list: list of paths to bam files, in order of forward, reverse, both.
                   Paths themselves are strings.
    cons_path_list: list of paths to consensus fastas, in order of forward, reverse, both.
                    Paths themselves are strings.
    out_path: path the figure will be output to.
    Optional:
    n: (default 10) number of top alignments by read count.
    verbose: (defualt False) display runtime information as a boolean
    
    Returns:
    string path to saved figure.
    """
    
    bam_for = bam_path_list[0]
    bam_rev = bam_path_list[1]
    bam_both = bam_path_list[2]
    
    cons_for  = cons_path_list[0]
    cons_rev = cons_path_list[1]
    cons_both = cons_path_list[2]
    
    ent_df_for = wrapper(cons_for, bam_for, False, n, verbose)
    ent_df_rev = wrapper(cons_rev, bam_rev, True, n, verbose)
    ent_df_both = wrapper(cons_both, bam_both, None, n, verbose)
    
    
    fig, ax = plt.subplots(len(ent_df_both),1, figsize = (20, 5*(len(ent_df_for))))

    for i in range(len(ent_df_both)):
    
    
        plotting_ser_for = pd.Series(ent_df_for['sh_ent'][i])
        base_indexes_for = list(plotting_ser_for.index)

        for j in base_indexes_for:
            base_indexes_for[j] = str(base_indexes_for[j])

        plotting_ser_rev = pd.Series(ent_df_rev['sh_ent'][i])
        base_indexes_rev = list(plotting_ser_rev.index)

        for j in base_indexes_rev:
            base_indexes_rev[j] = str(base_indexes_rev[j])

        plotting_ser_both = pd.Series(ent_df_both['sh_ent'][i])
        base_indexes_both = list(plotting_ser_both.index)

        for j in base_indexes_both:
            base_indexes_both[j] = str(base_indexes_both[j])

        ax[i].scatter(x=base_indexes_for, y=plotting_ser_for, s=10, c='tab:blue')
        ax[i].scatter(x=base_indexes_rev, y=plotting_ser_rev, s=10, c='tab:red')
        ax[i].scatter(x=base_indexes_both, y=plotting_ser_both, s=10, c='magenta')
        ax[i].grid(axis="x")
        ax[i].scatter(x=[0, 0], y=[0, 0], marker = "^", c='tab:orange')
        #ax[i].set_title("Shannon entropy per base")
        ax[i].set_xlabel("Base Index")
        ax[i].set_ylabel("Shannon Entropy")
        ax[i].set_ylim(0, 3.2)
        ax[i].tick_params(axis='x',labelrotation=-90)
        ax[i].legend(['forward', 'reverse','both'])

    outfile_path = '{}.png'.format(out_path)
    fig.savefig(outfile_path)
    
    return outfile_path
