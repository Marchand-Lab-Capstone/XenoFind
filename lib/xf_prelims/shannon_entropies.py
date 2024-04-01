import pysam
import os
from Bio import SeqIO
import pandas as pd
import tables
import h5py 
import re
import matplotlib.pyplot as plt
import numpy as np


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


def load_in_data(samfile_path):
    """
    load_in_data takes in a path to a basecalled and aligned consensus samfile,
    and converts the data into a dataframe format.

    Parameters:
    sam_path: path, as str, to the basecalled and aligned consensus samfile.

    Returns:
    pandas dataframe of the relevant data. [sequence id, sequence, quality, length, reference index]
    """

    # Set an empty list to hold the data
    data_list = []

    # Open the samfile
    with pysam.AlignmentFile(samfile_path, 'r') as samfile:

        # Loop through each read in the samfile
        for read in samfile:

            # If the read is mapped, 
            if (not read.is_unmapped):

                # Get the sequence id, alignment sequence, cigar string, reference start position
                # and the quality, and put it into a dictionary
                sequence_id = read.query_name
                seq = read.query_alignment_sequence
                cigar = read.cigartuples
                ref = read.reference_start
                quality = read.query_alignment_qualities
                read_dict = {'seq_id':sequence_id,
                             'sequence':seq,
                             'quality':quality,
                             'cigar':cigar,
                             'len':len(seq),
                             'ref':ref}
                #if (ref == 0):

                # Append the dictionary to the data list
                data_list.append(read_dict)

    # convert the list of dicts into a dataframe
    df = pd.DataFrame(data_list)

    # convert the qualities to list.
    df['quality'] = df['quality'].apply(list)

    # return the dataframe
    return df


def shift_sequence(sequence, qualities, ref, cigar):
    """
    shift_sequences takes in a sequence, corresponding quality, and reference position, and
    shifts it to align it with the proper positioning.

    Parameters:
    sequence: DNA read sequence as str.
    quality: DNA read qualities as list, same length of str.
    ref: reference position on the reference read, as int.
    cigar: the cigar string of the read, as str.

    returns:
    a touple of shifted sequence and shifted qualities.
    """

    # Get the sequence as a list
    seq = list(sequence)
    quals = qualities

    # set current index to zero,
    # make a list to hold all of the operations.
    current_index = 0
    operations = []

    # set up an empty list for the true sequence and qualities
    true_seq = []
    true_quals = []

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
                
    #print(''.join(operations))

    # Get the indexes of the operations with an offset
    index_minus_1 = [i for i, x in enumerate(operations) if x == "S"]

    # For each index in the offset index list,
    for index in index_minus_1:

        # insert a "deletion" at each index where offset.
        seq.insert(index, "-")
        quals.insert(index, -1)

    # get each index in the operations that has a deletion
    index_2 = [i for i, x in enumerate(operations) if x == "2"]

    # for each deletion index,
    for index in index_2:

        # insert a deletion at the index of deletion,
        # and a quality of -1 at the deletion site. 
        seq.insert(index, "-")
        quals.insert(index, -1)

    # for each index that has an insertion,
    index_1 = [i for i, x in enumerate(operations) if x == "1"]

    # Loop through the indexes that have insertions
    for index in index_1:

        # Delete the inserted base to shift everything over
        del seq[index]

        # set an addition symbol at the previous sequence
        seq[index-1] = "+"

    #print(''.join(seq))

        # 3 - skipped region from reference
        # 4 - soft clipping
        # 5 - hard clipping
        # 6 - padding
        # 7 - sequence match
        # 8 - sequence mismatch

    # Return the sequence and the qualities
    return seq, quals


def trim(vals, i, j):
    """
    trim takes in a starting and ending index, as well as a sequence of values,
    and then trims it to just be the values between those inexes.

    Parameters:
    vals: values as a list or a string to be trimmed
    i: starting index, as an int
    j: ending index, as an int

    Returns:
    either a list or a string of the sliced values, depending on what was passed
    """

    # Check if the type is a string, if it is, perform it and return as a string
    if type(vals) == type(""):
        return ''.join(list(vals)[i:j])

    # else if it is a list, return it as a list
    elif type(vals) == type([]):
        return vals[i:j]


def format_sequences(reads_df, indexes):
    """
    format_sequences takes in a dataframe of data previously exported from a sam file, and 
    indexes of the start and end of an unknown region, and returns the reads and quality scores for those segments.

    Parameters:
    reads_df: dataframe containing [sequence id, sequence, quality, length, reference index]
    indexes: a dictionary {"start", "end"} that contains the start and end indexes of an unknown region.
    """

    #step 1: remedying insertions, deletions and alignment

    # Duplicate the reads dataframe
    l_df = reads_df.copy()

    # Shift each read using the shift_sequence method
    l_df['shift'] = l_df.apply(lambda x: shift_sequence(x['sequence'],
                                                        x['quality'],
                                                        x['ref'],
                                                        x['cigar']),
                               axis = 1)

    # Create new columns and split the shift touple column to two individual columns
    new_col_list = ['shifted_seq','shifted_qual']
    for n,col in enumerate(new_col_list):
        l_df[col] = l_df['shift'].apply(lambda x: x[n])

    # Delete the shift column
    l_df = l_df.drop('shift',axis=1)

    #Step 2, trimming off the ends

    # Trim the sequences and qualities using the trimming method
    trim_sequence_series = l_df.apply(lambda x: trim(x['shifted_seq'],
                                                     indexes["start"],
                                                     indexes["end"]),
                                      axis = 1)
    trim_qual_series = l_df.apply(lambda x: trim(x['shifted_qual'],
                                                 indexes["start"],
                                                 indexes["end"]),
                                  axis = 1)

    #Step 3: formatting to an puput dataframe

    # Generate the output dataframe columns and assign the columns with the data
    out_df = pd.DataFrame(columns = ['id', 'seq', 'qual'])
    out_df['id']=l_df['seq_id']
    out_df['seq']=trim_sequence_series
    out_df['qual']=trim_qual_series

    # return the dataframe
    return out_df


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
    with base positons as columns and base id's as rows. 
    """

    # Create a list to hold each option for a base
    base_options = ['A', 'G', 'C', 'T', '-', '+']

    # Generate an empty list to hold the probabilities
    list_of_prob_dicts = []

    # for each column,
    for column in list(sequence_coordinates.columns):

        # get the number of each value in the column as a dictionary
        counts_dict = sequence_coordinates[column].value_counts().to_dict()

        # Loop through each option in the base options
        for option in base_options:

            # if the option isn't in the counts dictionary, set the value to zero
            if option not in counts_dict:
                counts_dict[option]=0

        # Sum up the amount of the base
        total = sum(counts_dict.values(), 0.0)

        #print(total)

        # Add each base to the probability dictionaries
        probs_dict = {key: value / total for key, value in counts_dict.items()}

        # append the probability dictionary of the bases to the list
        list_of_prob_dicts.append(probs_dict)

    # convert the list of dicts to a dataframe
    return_df = pd.DataFrame(list_of_prob_dicts)
    
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

    # loop through each base option in the probabilities
    for column in list(tmp_df.columns):

        # get the series of the column
        base_ser = tmp_df[column]

        # Compute the individual shannon entropy of the column
        tmp_df[column] = base_ser*np.log2(base_ser)

    # Sum all of the shannon entropies for each row, and *-1 to make them positive.
    entropy = tmp_df.sum(axis=1)*-1

    # Return the entropy
    return entropy
    

def file_shannon_entropy(sam_path, consensus_fasta = None, trimmed="N"):
    """
    file_shannon_entropy takes in sam path, and optional value for
    trimming and a sam reference,
    and provides the shannon entropy for the basecalled sequence.

    Parameters:
    sam_path: path to aligned samfile with basecalled reads, as a str.
    consensus_fasta: path to a fasta file containing labelling
                     that is output by consensus_methods.py.
                     default value of None.
    trimmed: flag for if trimming should happen to restrict the shannon entropy
             solely to the region of interest. Default value of "N" ('N'-No, 'Y'-yes)

    Returns: 
    dataframe of shannon entropy at each position along the reads.
    """

    loaded_df = None
    try:
        # Load in the data
        loaded_df = load_in_data(sam_path)
    except IOError:
        print("Samfile read failed. Did you pass the correct path?")

    # set up variable for read_indexes
    read_indexes = None

    # Check if all the info is available for trimming
    if (consensus_fasta != None) and (trimmed == "Y"):

        try:
            # if it is, get the trim indexes
            read_indexes = get_record_indexes(consensus_fasta)[0]
        except IOError:
            print("Consensus fasta index read failed. Is it the correctly indexed type or path?")
        read_indexes['start']=read_indexes['start']-1
        read_indexes['end']=read_indexes['end']+2
    else:
        # Otherwise, go from the beginning to the end.
        read_indexes = {'start':0,'end':loaded_df['sequence'].apply(len).max()}

    
    # get the formatted and trimmed data
    formatted_df = format_sequences(loaded_df, read_indexes)

    # convert the formatted sequence and quality data to per-base position columns
    seq_coordinate_df = pd.DataFrame(formatted_df['seq'].apply(list).to_list())
    quality_coordinate_df = pd.DataFrame(formatted_df['qual'].to_list())

    # calculate the base probabilities of the sequence coordinates:
    base_probs_df = calc_base_probs(seq_coordinate_df)

    # Get the series of the shannon entropy
    shentropy_series = shannon_entropies(base_probs_df)

    #return the shannon entropy series
    return shentropy_series


def generate_shannon_plot(series, ax, indexes=None):
    """
    generate_shannon_plot takes in a series of shannon entropies,
    a matplotlib ax, and plots the shannon entropies.

    Parameters:
    series: shannon entropies per base as a pandas series
    ax: maplotlib axes object
    indexes: a list of indexes that are the edge of the series of interest. Default None.

    Returns:
    a matplotlib axes of the shannon entropy plot
    """

    # Get the base indexes for x values
    base_indexes = list(series.index)
    for i in base_indexes:
        base_indexes[i] = str(base_indexes[i])

    ax.scatter(x=base_indexes, y=series)
    if (indexes != None):
        ax.scatter(x=indexes, y=[0,0])
    ax.grid(axis='x')

    return ax


def main():
    samfile = input("Please provide path to basecalled samfile: ")
    conspath = "E"
    trimQ = "E"

    consensusQ = "E"
    while (list(consensusQ)[0] != "N") and (list(consensusQ)[0] != "Y"):
        #print(list(consensusQ)[0] != "Y")
        consensusQ = input("Do you have a consensus file with index headers? (Y/N): ")

    if (list(consensusQ)[0] == 'Y'):
        conspath = input("Please provide valid path to consensus file: ")

        while (list(trimQ)[0] != "N") and (list(trimQ)[0] != "Y"):
            trimQ = input("Should the unknown region be isolated? ")

    return file_shannon_entropy(samfile, conspath, trimQ)


if __name__ == '__main__':
    main()
    
        


    

    
    