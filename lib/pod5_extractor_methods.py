import sys
import numpy as np
from alive_progress import alive_bar
import pod5
import multiprocessing
import psutil
import datetime
import os
import gc

def yoink_from_read(read):
    '''
    yoink_from_read takes in a pod5.Reader.Read object,
    and returns a dict of relevant info:
    
    Parameters:
    read: the Read Object to be parsed
    
    Returns:
    dictionary with seq_id, signal, and freq as keys.
    '''
    
    # save the sequence id, signal, and frequency to a dict
    seq_id = read.read_id
    signal = read.signal_pa
    freq = read.run_info.sample_rate
    data_dict = (str(seq_id),
                 np.array(list(signal), dtype=np.float32),
                 np.float16(freq/1000))
    return data_dict

def subprocess(batch, path):
    '''
    subprocess takes in a batch and a path to a pod5 file, and extracts the
    information at the indexes in that batch.
    
    Parameters:
    batch: a list of values corresponding to an index range in a pod5 file.
          the initial value in the list will be used as the first value in the slice,
          the final will be the last value in the slice: [batch[0]:batch[-1]]
    path: path, as str, to a pod5 file.
    
    Returns:
    a list of dicts of the reads at the indexes in the defined range.
    
    '''
    # get the initial and final indexes
    initial = batch[0]
    final = batch[-1]
    
    # generate an empty sublist for storing the read data
    subdata_list = []
    
    # loop through each read in the index range in the pod5 file
    for read in list(pod5.Reader(path).reads())[initial:final]:
        
        # append the read data dict to the subdata list.
        subdata_list.append(yoink_from_read(read))
        
    
    # return the subdata list
    return subdata_list

def flatten_list(list_of_lists):
    '''
    flatten_list takes a list of lists and makes it a list of all the values in that list.
    
    Parameters: 
    list_of_lists: a list of lists.
    
    Returns: 
    a list with all the items in the sublist becoming the items in the superlist. 
    '''
    flattened_list = [item for sublist in list_of_lists for item in sublist]
    return flattened_list

def mem_check(data_list):
    print('check called')
    sum_val = 0
    for j in range(len(data_list)):
        if (j%10000 == 0):print(j)
        a=sys.getsizeof(data_list[j][0])
        b=sys.getsizeof(data_list[j][2])
        c=0
        for i in range(len(data_list[j][1])):
            c+=sys.getsizeof(data_list[j][1][i])
        sum_val +=(a+b+c)/1073741824
    return sum_val

def extract_pod5_dict(p5_path, index_batch_size, multiprocess_pool_size):
    '''
    Extract_pod5_dict takes in a path, an index batch size, and multiprocessing pool size, 
    and then extracts the information from a pod5 file without consuming an excess of memory.
    
    Parameters:
    p5_path: path to a pod5 file, as str
    index_batch_size: size of the number of reads to be extracted within a pool, as an int.
    multiprocess_pool_size: number of pools to be used, as an int.
    
    Returns:
    a list of dicts of the relevant read information extracted from the pod5_file. 
    '''
    # set up an empty list to contain the data dictionaries
    data_list = []
    
    # get the length of the pod5 file by counting the reads
    pod5_length = len(list(pod5.Reader(p5_path).reads()))
    
    
    # get the indexes of each read by arranging a scale from 0 to pod5_length
    pod5_indexes = np.arange(0, pod5_length, 1)
    
    # take the indexes and convert them into batches of indexes, in accordance with the batchsize that was passed
    batched_indexes = [pod5_indexes[i:i+index_batch_size] for i in range(0, len(pod5_indexes), index_batch_size)][:6]
    
    index = 0
    # start a loading bar
    with alive_bar(len(batched_indexes)) as pod5_loading_progress:

        # get the batch size of the pools
        batched_batch_size = multiprocess_pool_size
        
        # generate batches of each batch of index based on pool size
        batch_of_batches = [batched_indexes[i:i+batched_batch_size] for i in range(0, len(batched_indexes), batched_batch_size)]

        # loop through each sub batch in the multiprocessing batches
        for sub_batch in batch_of_batches:
            

            # generate a list the same length of the sub batch containing just the pod5 path
            repeating_path_list = list(np.full(fill_value=p5_path, shape=len(sub_batch)))
            
            # create a list of touples containing the sub)batch and the path 
            sub_batch = list(zip(sub_batch, repeating_path_list))

            # use multiprocessing to get the batch sizes
            with multiprocessing.Pool(batched_batch_size) as pool:

                # extend the data list with the results from the pooling of the subprocess
                data_list.extend(flatten_list(pool.starmap(subprocess, sub_batch)))
                gc.collect()
            # progress updates
            for i in range(batched_batch_size):
                pod5_loading_progress()
                index += 1
            print('{}/{}'.format(index, len(batched_indexes)))
            #print('{} gb'.format(mem_check(data_list)))
            memory_report()
            gc.collect()

                
    # return the data list. 
    return data_list


def memory_report():
    '''
    memory_report is a method of reporting the memory use statistics at any given point in the process,
    through printing a string of information.
    
    Parameters:
    None
    Returns:
    None
    '''
    
    # use psutil to get the process memory info in Gb
    proc = psutil.Process(os.getpid())
    rss = proc.memory_info().rss/(1*10**9)
    vms= proc.memory_info().vms/(1*10**9)
    shared = proc.memory_info().shared/(1*10**9)
    data = proc.memory_info().data/(1*10**9)
    print("[ Memory Report {} ]------------------------------------------------".format(datetime.datetime.now()))
    print("RSS: {} Gb, VMS: {} Gb, SHARED: {} Gb, DATA: {} Gb".format(rss, vms, shared, data))
    print("--------------------------------------------------------------------------------------------")


if __name__ == '__main__':
    p5_path = "/media/sebastian/Slepnir/xenofind_datasets/Working_directory_PZ/model_training/merged_pod5/merged.pod5"

    extract_pod5_dict(p5_path, 20000, 3)