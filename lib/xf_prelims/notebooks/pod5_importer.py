import os
import sys
import pod5
import gc
import multiprocessing
import pandas as pd

def load_pod5_data(p5_path):
    print('pod5 load called...')
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
    
    
    # ----------------------------------------------- Segment that hogs 7gigs. 
    # use pod5's reader object to read through the pod5 file
    with pod5.Reader(p5_path) as pod5_file:
        
        for read in pod5_file.reads():
            data_list.append(yoink(read))
        pod5_file.close()
        
    print('saving to pod5.pickle...')
    pd.DataFrame(data_list).to_pickle('pod5.pickle')
    
    return pd.DataFrame(data_list)
            
    # return the data list 
    #return data_list



def yoink(read):
    seq_id = read.read_id
    signal = read.signal_pa
    freq = read.run_info.sample_rate
    data_dict = {'seq_id': str(seq_id),
                 'signal': signal,
                 'freq': freq}
    return data_dict

        

if __name__ == "__main__":
    import argparse
    print("running...")
    parser = argparse.ArgumentParser(description='pod5 loading')
    parser.add_argument('path', metavar='path', help='Path to Pod5 file.')
    args = parser.parse_args()
    sys.exit(load_pod5_data(args.path))
    