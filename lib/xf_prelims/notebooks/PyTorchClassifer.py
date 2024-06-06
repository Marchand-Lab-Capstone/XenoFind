# ------------------------ MASS IMPORT OF JUNK THAT I MAY OR MAY NOT NEED 
# First, we import all the relevant packages:
# Base packages for data processing and cleanliness
import os
import sys
import numpy as np
import subprocess
import multiprocessing
import itertools
import matplotlib.pyplot as plt
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


class PyTorchClassifier:
    '''
    Originally written by S. Peck and A. Mahmoud, Adapted for this project by S. Peck, and is licensed for this sole implementation. 
    
    PyTorchClassifier - Sebastian Peck, 12/1/2023
    PyTorchClassifier is an sklearn-styled class to aid in
    the classification of data using a pytorch neural network.
    Initialization requires the device to run pytorch on,
    number of input features, and number of classes.
    
    Updated 6/4/2024 - S. Peck
    UPDATED 6/5/2024 - S. PECK
    
    '''
    def __init__(self,
                 device,
                 n_features:int,
                 out_features:int,
                 nl_list:list=[5, 5, 5],
                 activation_function=nn.ReLU(),
                 loss_function=nn.CrossEntropyLoss(),
                 n_epochs:int=1000,
                 learning_rate:float=.1,
                 random_state:int=42,
                 class_weights:list=[0]):
        '''
        initializing PyTorchClassifyer sets pytoch's random
        seed and generates the loss function, optimizer,
        and model of neural network based on passed parameters.
        
        Parameters:
        device: the device running pytorch
        n_features: number of features in data as int
        out_features: number of classes in the data as int
        nl_list: a list, with number of indexes as layers, 
                 and values representing neurons at that layer(as ints) Default = [5, 5, 5]
        activation_function: the torch.nn activation function of choice. Default = ReLU()
        loss_function: the torch.nn loss function of choice. Defualt = CrossEntropyLoss()
        n_epochs: int of number of epochs to perform when training. Default= 1000
        learning_rate: float representing learning rate. Defualt = .1
        random_state: the random state to use. Default = 42
        class_weights: List of class weights USED FOR BOOKKEEPING AND DOES NOT ACTUALLY PASS TO LOSS FXN.

        '''
        
        # Generate a string of the parameters
        self._param_string = (str(device) +'|'+ str(n_features) +'|'+ str(out_features) +'|'+ str(nl_list) +'|'+ str(activation_function) +'|'+str(loss_function)+'|'+ str(n_epochs) +'|'+ str(learning_rate) +'|'+ str(random_state) + '|' + str(class_weights))
        
        # assign the random seed to pytorch
        torch.manual_seed(random_state)
        
        # set the class variables
        self.epochs = n_epochs
        self.activation_function = activation_function
        self.loss_function = loss_function
        self.device = device
        self.n_features = n_features
        self.random_state = random_state
        
        # Generate the model using the local ClassifierNetwork class
        self.model = (self.ClassifierNetwork(input_features=n_features,output_features=out_features,nl_list = nl_list, activation_function = self.activation_function)).to(device)
        
        # Optimize the model using SGD
        self.optimizer = torch.optim.SGD(params = self.model.parameters(), lr=learning_rate)
        

    
    def train_model(self, dataloader, epochs=None):
        if type(epochs) == type(None):epochs=self.epochs
        '''
        train_model takes in a torch formatted dataset of features,
        and the corresponding torch formatted dataset of classes,
        and trains the current model version on that data.
        
        Parameters:
        dataloader: a TensorDataset DataLoader object containing the data to be processed.
        epochs: an int representing number of epochs to train this dataset. Default is the passed epochs when creating the model.
        
        Returns:
        The trained model.
        '''
        
        # Set up the model to train
        self.model.train()
        
        # repeat for every epoch
        for epoch in range(epochs):
            
            
            # load a subset from the dataloader
            for i, (x, y) in enumerate(dataloader):
                
                # generate the features and classes
                features, classes = x.to(self.device), y.to(self.device)

                # Forward Pass - pass the features and
                # convert to class logits
                class_logits = self.model(features)

                # Convert the logits to probabilities to labels using softmax
                class_predictions = torch.softmax(class_logits, dim=1).argmax(dim=1)

                # Calculate loss
                # passing the logits and the training classes
                loss = self.loss_function(class_logits, classes)

                # Reset the optimizer gradient
                self.optimizer.zero_grad()

                # set backwards loss training
                loss.backward()

                # step the optimizer by one
                self.optimizer.step()
            print("{}  {}/{}         ".format(percenttobar(epoch/epochs), epoch, epochs), end="\r")
            
        # return the model
        return self.model
    
    
    def save_model_state(self, path, filename):
        '''
        save_model saves the model state to the given path.
        
        Parameters:
        path: path, as str, to the save location directory, must end with '/'
        filename: name of the file to be saved, as str
        
        Returns:
        path to saved file, as str
        '''
        
        # generate savefile strings for the state and parameters
        save_path_str = '{}{}/state.pt'.format(path, filename)
        save_path_str_2 = '{}{}/params.txt'.format(path, filename)
        
        if not os.path.isfile('{}{}'.format(path, filename)):
            os.mkdir('{}{}'.format(path, filename))
        
        # save the state
        torch.save(self.model.state_dict(), save_path_str)
        
        # save the parameters
        with open(save_path_str_2, 'w') as f:
            f.write(self._param_string)
        
        # return the directory
        return save_path_str
    
    
    def load_model_state(self, path):
        '''
        load_model_state loads the model state from a given path.
        
        Parameters:
        path: Path, as str, to the saved parameters file
        
        Returns:
        the loaded model, in evaluation state
        '''
        
        # load the model from the path
        self.model.load_state_dict(torch.load(path))

        # set the model into evaluation mode
        self.model.eval()
        
        # return the model
        return self.model
            
    
    def accuracy_score(self, y_true, y_pred):
        '''
        accuracy_score takes in two torch tensors of true and predicted values
        and returns the accuracy as a percentage fraction.
        
        Parameters:
        y_true: pytorch tensor of LongTensors representing classes of each point
        y_pred: pytorch tesnor of Longtensors representing predicted classses of each point
        
        Returns:
        the accuracy fraction
        '''
        
        # get number of the true predicted values to get the valid ones using torch.eq()
        valid = torch.eq(y_true, y_pred).sum().item()
        
        # Generate the precentage fraction
        acc = (valid/len(y_pred))
        
        # return the accuracy
        return acc
    
    
    def test_model(self, features, classes):
        '''
        test_model takes in a set of features and corresponding classes to
        test how the model performs with that dataset.
        
        Parameters:
        features: a pytorch tensor of floats representing the data for each feature
        classes: a pytorch tensor of LongTensors representing the classes for each point
        
        Returns:
        a touple containing the loss and the accuracy.
        
        '''
        
        # generate the logit values from the model by passing the test features
        test_logits = self.model(features)
        
        # generate the predictions using torch.softmax() and the test logits
        test_predictions = torch.softmax(test_logits, dim=1).argmax(dim=1)
        
        # Generate the loss using the class's loss function
        test_loss = self.loss_function(test_logits,
                                       classes)
        
        # Generate the accuracy score using the class accuracy score function
        test_acc = self.accuracy_score(classes,
                                       test_predictions)
        
        # Generate a string to hold the results for pretty-printing if needed
        #results = "Loss: {}, Accuracy Score: {}".format(test_loss, test_acc)
        #print(results)
        
        # return a touple of the loss value and the test accuracy.
        return (test_loss.item(), test_acc)
    
        
    class ClassifierNetwork(nn.Module):
        '''
        ClassifierNetwork is a subclass of the Pytorch nn.Module,
        which is meant to be used as a classifier.
        '''
        def __init__(self,
                     input_features:int,
                     output_features:int,
                     nl_list:list=[5, 5, 5],
                     activation_function=nn.ReLU()):
            '''
            Extends nn.Module.
            
            initializer sets up the model stack in sequential order according to the passed parameters.
            Uses linear layer stack of linear nn layers.
            
            Parameters:
            input_features: an int representing the number of features in the dataset
            output_features: an int representing the number of classes in the dataset
            nl_list: list representing number of neurons per layer,
                     with n indexes as the number of layers, Default = [5, 5, 5]
            activation_function: activation function of choice between layers. Default=nn.ReLU()
            
            Returns:
            N/A
            '''
            
            # Initialize the superclass
            super().__init__()
            
            # generate the starting neurons as the first index in the list
            starting_neurons = nl_list[0]
            
            # set up a sequential stack
            self.linear_layer_stack = nn.Sequential()
            
            # add the first layer where the input features is the first number of neurons, and the output features are the starting neurons
            self.linear_layer_stack.append(nn.Linear(in_features = input_features, out_features = starting_neurons))
            
            
            # loop through each layer of the nl_list.
            for i in range(len(nl_list)-1):
                # append the activation function to the linear layer stack
                self.linear_layer_stack.append(activation_function)

                # append a new layer where the input features is the current index in the nl_list, and the subsequent is the next.
                self.linear_layer_stack.append(nn.Linear(in_features=nl_list[i], out_features = nl_list[i+1]))
                
            # add the output layer
            self.linear_layer_stack.append(activation_function)
            self.linear_layer_stack.append(nn.Linear(in_features=nl_list[-1], out_features = output_features))


        def forward(self, x):
            '''
            forward takes in a set of features and applies it to the linear layer stack.
            
            Paramters:
            x: set of features. 
            
            Returns:
            the linear layer stack with the passed features.
            '''
            return self.linear_layer_stack(x)
        
        
    def __str__(self):
        '''
        export self as string.
        '''
        return self._param_string