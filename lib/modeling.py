'''
modeling.py
S. Peck

modeling.py contains methods critical for running an evaluation of XNA location using 
the machine learning methods from a given PyTorchClassifier. 

run_check_device() - checks the current device pytorch is running on, or can run on
percenttobar() - converts a percentage value to a progress bar string
activ_func_parser() - parses an activity function string.
loss_func_parser() - parses a loss function from parameters
get_model_params() - generates the model parameter tuple from a directory.
load_models() - loads the pytorch classifier models from the directories passed.
load_pcas() - loads the window model and the base model PCA as exported from model training. NOT GOOD CODE OR ML TODO
window_detection() - uses a pytorcclassifier to detect XNA windows from a dataframe of features
windowed_base_detection() - uses a pytorchclassifier to detect XNA within a window from a dataframe of features.
'''

# First, we import all the relevant packages:
# Base packages for data processing and cleanliness
import os
import sys
import numpy as np
import itertools
import pandas as pd
import random
import datetime

# adjust the path to include the directory for our own python files
sys.path.append('..//..//model_gen')
import feature_extraction as fe

# Then,the machine learning packages
import torch
from torch import nn
from sklearn.model_selection import train_test_split
from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA
from sklearn.metrics import confusion_matrix
import seaborn as sns


# -------------- END MASS IMPORT OF JUNK

import PyTorchClassifer as ptc


# METHODS FOLLOW

import ast
import pickle as pk
DEVICE=(
    "cuda"
    if torch.cuda.is_available()
    else "mps"
    if torch.backends.mps.is_available()
    else "cpu"
)

def run_check_device():
    '''
    run_check_device() checks which devices are available to run pytorch on and returns it as a string.
    Parameters:
    N/A
    
    Returns:
    string indicating torch capable device.
    '''
    d = (
    "cuda"
    if torch.cuda.is_available()
    else "mps"
    if torch.backends.mps.is_available()
    else "cpu")
    
    DEVICE = d
    return d


def percenttobar(frac):
    '''
    Percenttobar takes in a fraction and reutrns a progress bar.
    Parameters:
    frac: fraction, as a float
    
    Returns:
    string representation of progress out of 100% to get that fraction.
    '''
    bar_str = "|"
    max_bars = 20
    perc = frac*2000
    n_bars = int(perc/100)
    for i in range(n_bars):
        bar_str += "="
    for i in range(max_bars-n_bars):
        bar_str += " "
    bar_str += "|  {}%                ".format(round(frac*100, 3))
    return bar_str


def activ_func_parser(act_str):
    '''
    parses an activation function string. Currently only works with nn.tanh()
    
    Parameters:
    act_str: string containing the activation function from params
    
    Returns:
    currently only returns nn.Tanh() function
    
    TODO: Add more activity function options
    '''
    if act_str == 'Tanh()':
        return nn.Tanh()
    else:
        return nn.Tanh()

    
def loss_func_parser(loss_str, wt, device):
    '''
    parses a loss function string with weight and device. Currently only uses
    nn.CrossEntropyLoss()
    
    Parameters: 
    loss_str: the string representation of a loss function
    wt: list containing weights of classes
    device: device to pass tensors to
    
    Returns:
    if the weight array has more than 1 class, nn.CrossEntropyLoss with weights
    otherwise, nn.CrossEntropyLoss()
    
    TODO: 
    Add more loss functions
    '''
    if len(wt) > 1:
        wt = torch.tensor(wt).type(torch.float).to(device)
        if loss_str == 'CrossEntropyLoss()':
            return nn.CrossEntropyLoss(weight = wt)
    else:
        return nn.CrossEntropyLoss()


def get_model_params(model_dir):
    '''
    generates the model paramater touple from a given model directory.
    assumes model directory contains 'params.txt'.
    
    Parameters:
    model_dir: directory to model parameters
    
    Returns: a tuple of the parameters for the model
    '''
    DEVICE = run_check_device()
    param_string = None
    with open(model_dir+"params.txt", 'r') as f:
        param_string = f.readline()

    param_list = param_string.split("|")

    # Check if the device is matching the current device and update accordingly.

    if param_list[0] == DEVICE:
        DEVICE = param_list[0]
    else:
        param_list[0] = DEVICE

    param_list[9] = ast.literal_eval(param_list[9]) # weights notation

    param_list[1] = int(param_list[1]) # PCA features
    param_list[2] = int(param_list[2]) # Classes
    param_list[3] = ast.literal_eval(param_list[3]) # Network list
    param_list[4] = activ_func_parser(param_list[4]) # activation function
    param_list[5] = loss_func_parser(param_list[5], param_list[9], DEVICE) # loss function
    param_list[6] = int(param_list[6]) # epochs
    param_list[7] = float(param_list[7]) # learning rate
    param_list[8] = int(param_list[8]) # random state


    param_list = param_list[:-1]

    return tuple(param_list)


def load_models(window_model_dir, base_model_dir):
    '''
    this is very hastily written code and should be revised to be better. 
    forgive me :( -S
    
    load_models takes in a window model and base model directory,
    and then loads the parameters to generate PyTorchClassifiers with them
    
    Parameters:
    window_model_dir: path, as str, to window model directory
    base_model_dir: path, as str, to base_model directory
    
    Returns:
    pytorchclassifiers of( window model, base model) as tuple
    
    TODO: 
    Make individual so that it doesnt require both base and window.
    
    '''
    run_check_device()
    # open params:
    window_model_params = get_model_params(window_model_dir)
    base_model_params = get_model_params(base_model_dir)
    
    #generate the models
    loaded_window_model = ptc.PyTorchClassifier(*window_model_params)
    loaded_base_model = ptc.PyTorchClassifier(*base_model_params)
    
    #load the model states
    loaded_window_model.load_model_state(window_model_dir + 'state.pt')
    loaded_base_model.load_model_state(base_model_dir + 'state.pt')
    
    return loaded_window_model, loaded_base_model


def load_pcas(window_model_dir, base_model_dir):
    '''
    Dave, I know this is bad code and bad ML but genuinely i am writing this on wednesday at 8:10 pm -S
    
    load_pcas tkaes in a window model directory, a base model directory,
    and then loads the PCA for each that was exported.
    
    Parameters:
    window_model_dir: path to window model directory as str
    base_model_dir: path to base model directory as str
    
    Returns:
    tuple of windowed pca, based pca.
    
    TODO:
    REMOVE THIS, PCA SHOULD BE RUN ON EACH INDIVIDUAL DATASET
    '''
    pca_window = pk.load(open(window_model_dir+"pca.pkl",'rb'))
    pca_base = pk.load(open(base_model_dir+"pca.pkl", 'rb'))
    return pca_window, pca_base


def window_detection(window_model, read_feature_df, pca):
    '''
    uses a window model to detect xna windows from the data
    
    Parameters:
    window_model: pytorchclassifier of the window model
    read_feature_df: dataframe of all the read features for each base
    pca: the pca used when training the window model
    
    Returns:
    a list of tuples representing the windows XNA were detected in using the model.
    '''
    DEVICE = run_check_device()
    # step 1: split the data into windows of size 7 ----------------
    window_size = 7
    windows = [read_feature_df[i:int(i+window_size)] for i in range(0, len(read_feature_df), int(window_size))]
    
    window_classes = []
    window_features = []
    
    # for each window minus the potentially inconsistently sized one,
    for base_window in windows[:-1]:
        
        # generate an empty list for the window's 1d features
        window_sub_features = []
        
        # check if there's an XNA, append a 1 to the classes if so. When testing, there should never be anything but zeroes here.
        if len(base_window[base_window['XNA_PRESENT'] > 0]) > 0:
            window_classes.append(1)
        else:
            window_classes.append(0)
            
        # loop through each base features, and extend the features to the 1d feature array
        for base in base_window.drop(columns=['XNA_PRESENT']).values:
            window_sub_features.extend(base)
            
        # append the window's 1d features to the list.
        window_features.append(window_sub_features)
    # step 2: scale the data.
    
    scaled_features = StandardScaler().fit_transform(window_features)
    
     # step 3as;lkvskvna;lkf: generate PCA - THIS IS BAD FORM AND SHOULD NOT BE DONE THIS WAY BUT TIME HAS FORCED MY HAND
    #n_pca_features = window_model.n_features
    #pca = PCA(n_components=n_pca_features, random_state = window_model.random_state)
    #pca = PCA.fit_transform(scaled_features)
    
    # step 4: generate the feature tensor and label tensor
    feature_tensor = torch.tensor(pca.transform(scaled_features)).type(torch.float).to(DEVICE)
    # --> unuused because it is only for training rn:  label_tensor = torch.tensor(np.asarray(window_classes)).type(torch.LongTensor).to(DEVICE)
    
    # step 5: run the model!
    window_logits = window_model.model(feature_tensor)
    
    # step 6: convert the model logits to predictions
    window_predictions = torch.softmax(window_logits, dim=1).argmax(dim=1)
    
    # step 7: extract the predicted windows:
    #convert the predictions to a list
    predictions = window_predictions.tolist()

    # convert the predictions to a series
    pred_series=pd.Series(predictions)
    
    # get the lsit of windows predicted
    idd_windows = list(zip(np.asarray(pred_series[pred_series> 0].index.tolist())*(window_size), (np.asarray(pred_series[pred_series> 0].index.tolist())+1)*(window_size)))
    
    return idd_windows
    

def windowed_base_detection(base_model, window, read_feature_df, pca):
    '''
    uses a base model to detect bases from a window
    
    parameters:
    base_model: pytorchclassifier trained on windowed bases
    window: tuple of lower and upper bound of window of detection
    read_feature_df: dataframe containing the read features for all bases in the reference
    pca: pca the base_model was trained on
    
    Returns:
    list of bases detected as XNA in the given window
    
    TODO:
    REMOVE PCA FROM PARAMS AND INCLUDE INSIDE
    '''
    DEVICE = run_check_device()
    
    # step 1: generate the detection region from the window
    window_bases = read_feature_df[window[0]:window[1]]
    
    # step 2: get the features and classes
    window_base_features = window_bases.loc[:, window_bases.columns != 'XNA_PRESENT']
    window_base_features = window_bases['XNA_PRESENT']


    # step 3: scale the data:
    base_scaler = StandardScaler()
    
    # For some heavens forsaken reason, this is the only way this runs and it wont work separated
    base_scale = StandardScaler().fit_transform(window_bases.loc[:, window_bases.columns != 'XNA_PRESENT'])

        
    # step 4: generate PCA - THIS IS BAD FORM AND SHOULD NOT BE DONE THIS WAY BUT TIME HAS FORCED MY HAND
    #n_pca_features = base_model.n_features
    #pca = PCA(n_components=n_pca_features, random_state = window_model.random_state)
    #pca.fit(scaled_features)
    
    # step 5: generate the feature and label tensors
    region_features = torch.tensor(pca.transform(base_scale)).type(torch.float).to(DEVICE)
    # UNUSED    region_labels = torch.tensor(np.asarray(xna_region_classes)).type(torch.LongTensor).to(device)
        
    # step 6: Run the model! 
    region_logits = base_model.model(region_features)
        
    # step 7: get the predicted base classes
    xna_predictions = torch.softmax(region_logits, dim=1).argmax(dim=1)
        
    # step 8: get the list of predictions
    xna = xna_predictions.tolist()

    # get the series of predicted bases
    pred_xna_series = pd.Series(xna)
    
    # return a list of the bases id'd as XNA. 
    return (np.asarray(pred_xna_series[pred_xna_series >0].index.tolist()) + window[0]).tolist()
                     
