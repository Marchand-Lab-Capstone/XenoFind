import os
import sys
import torch
import torch.nn as nn
import pandas as pd
import numpy as np
import warnings
from alive_progress import alive_bar
from sklearn.model_selection import train_test_split
import setup_methods as setup  
import feature_extraction as fe
import argparse
import datetime
import json

os.environ["OMP_NUM_THREADS"] = "1"  # Set OpenMP threads to 1 to reduce OpenBLAS warnings

class CustomPlusCoolDownScheduler:
    """
    Custom learning rate scheduler with a cooldown period after a specified number of epochs.

    Args:
        opt (Optimizer): Optimizer to be used.
        initial_lr (float): Initial learning rate.
        lr_scheduler_str (str): Learning rate scheduler type.
        lr_scheduler_kwargs (dict): Arguments for the learning rate scheduler.
        epochs (int): Number of epochs.
        cool_down_epochs (int): Number of cooldown epochs.
        cool_down_lr (float): Learning rate during cooldown.
    """
    def __init__(self, opt, initial_lr, lr_scheduler_str, lr_scheduler_kwargs, epochs, cool_down_epochs, cool_down_lr):
        self.last_epoch = -1
        self.opt = opt
        self.initial_lr = initial_lr
        self.num_epochs = epochs
        self.num_cool_down_epochs = cool_down_epochs
        self.cool_down_lr = cool_down_lr
        self.total_epochs = epochs + cool_down_epochs
        self.custom_scheduler = getattr(torch.optim.lr_scheduler, lr_scheduler_str)(
            self.opt,
            **{lr_key: lr_val for lr_key, lr_val in lr_scheduler_kwargs.items()}
        )

    def get_last_lr(self):
        """
        Get the last computed learning rate.

        Returns:
            list: List of learning rates.
        """
        if self.last_epoch < self.num_epochs - 1:
            return self.custom_scheduler.get_last_lr()
        else:
            return [self.cool_down_lr for _ in self.opt.param_groups]

    def step(self, metrics):
        """
        Update the learning rate at each step (epoch).

        Args:
            metrics (float): The value of the monitored metric (e.g., validation loss).
        """
        self.last_epoch += 1
        if self.last_epoch < self.num_epochs - 1:
            return self.custom_scheduler.step(metrics)
        else:
            for pg in self.opt.param_groups:
                pg["lr"] = self.cool_down_lr

class ConvLSTMCell(nn.Module):
    """
    Convolutional LSTM cell.

    Args:
        input_channels (int): Number of input channels.
        hidden_channels (int): Number of hidden channels.
        kernel_size (int): Size of the convolutional kernel.
        bias (bool): Whether to include a bias term in the convolution.
    """
    def __init__(self, input_channels, hidden_channels, kernel_size, bias=True):
        super(ConvLSTMCell, self).__init__()
        self.input_channels = input_channels
        self.hidden_channels = hidden_channels
        self.kernel_size = kernel_size
        self.padding = kernel_size // 2  # Ensures output has the same spatial dimensions as input
        self.bias = bias
        # Convolutional layer that combines input and previous states
        self.conv = nn.Conv2d(in_channels=input_channels + hidden_channels,
                              out_channels=4 * hidden_channels,
                              kernel_size=self.kernel_size,
                              padding=self.padding,
                              bias=self.bias)

    def forward(self, input_tensor, cur_state):
        """
        Forward pass of the ConvLSTMCell.

        Args:
            input_tensor (torch.Tensor): Input tensor at a time step (batch, channels, height, width).
            cur_state (torch.Tensor): Current hidden and cell state concatenated.

        Returns:
            torch.Tensor: Next hidden and cell state concatenated.
        """
        h_cur, c_cur = torch.split(cur_state, self.hidden_channels, dim=1)  # Split current states
        combined = torch.cat([input_tensor, h_cur], dim=1)  # Concatenate input and hidden state
        combined_conv = self.conv(combined)  # Apply convolution

        # Split the convolution output into four parts for gates
        cc_i, cc_f, cc_o, cc_g = torch.split(combined_conv, self.hidden_channels, dim=1)
        i = torch.sigmoid(cc_i)  # Input gate
        f = torch.sigmoid(cc_f)  # Forget gate
        o = torch.sigmoid(cc_o)  # Output gate
        g = torch.tanh(cc_g)  # Candidate cell state

        c_next = f * c_cur + i * g  # Update cell state
        h_next = o * torch.tanh(c_next)  # Update hidden state

        return torch.cat([h_next, c_next], dim=1)  # Concatenate next hidden and cell state

    def init_hidden(self, batch_size, shape):
        """
        Initialize hidden and cell states with zeros.

        Args:
            batch_size (int): Size of the batch.
            shape (tuple): Spatial dimensions (height, width) of the input tensor.

        Returns:
            torch.Tensor: Initialized hidden and cell states concatenated.
        """
        height, width = shape
        h = torch.zeros(batch_size, self.hidden_channels, height, width, device=self.conv.weight.device)
        c = torch.zeros(batch_size, self.hidden_channels, height, width, device=self.conv.weight.device)
        return torch.cat([h, c], dim=1)

class ConvLSTM(nn.Module):
    """
    Convolutional LSTM Network.

    Args:
        input_channels (int): Number of input channels.
        hidden_channels (int): Number of hidden channels.
        kernel_size (int): Size of the convolutional kernel.
        num_layers (int): Number of ConvLSTM layers.
        bias (bool): Whether to include a bias term in the convolutions.
    """
    def __init__(self, input_channels, hidden_channels, kernel_size, num_layers, bias=True):
        super(ConvLSTM, self).__init__()
        self.input_channels = input_channels
        self.hidden_channels = hidden_channels
        self.kernel_size = kernel_size
        self.num_layers = num_layers
        self.bias = bias
        # Create a list of ConvLSTM cells for each layer
        cells = nn.ModuleList([ConvLSTMCell(input_channels if i == 0 else hidden_channels,
                                            hidden_channels, kernel_size, bias) for i in range(num_layers)])
        self.cells = cells
        self.fc = nn.Linear(hidden_channels, 1)  # Fully connected layer for binary classification

    def forward(self, input_tensor, hidden_state=None):
        """
        Forward pass of the ConvLSTM layer.

        Args:
            input_tensor (torch.Tensor): Input tensor (batch, sequence_length, channels, height, width).
            hidden_state (torch.Tensor, optional): Initial hidden states for each layer.

        Returns:
            tuple: Output of the last layer and the last hidden states.
        """
        layer_output_list = []
        last_state_list = []
        seq_len = input_tensor.size(1)  # Sequence length
        cur_layer_input = input_tensor  # Initialize current layer input

        for layer_idx, cell in enumerate(self.cells):
            # Initialize hidden states if not provided
            if hidden_state is None:
                cur_state = cell.init_hidden(input_tensor.size(0), (input_tensor.size(3), input_tensor.size(4)))
            else:
                cur_state = hidden_state[layer_idx]
            output_inner = []

            for t in range(seq_len):
                cur_state = cell(cur_layer_input[:, t, :, :, :], cur_state)
                h_next, c_next = torch.split(cur_state, self.hidden_channels, dim=1)
                output_inner.append(h_next)

            layer_output = torch.stack(output_inner, dim=1)  # Stack outputs
            cur_layer_input = layer_output  # Set current layer input to output of this layer

            layer_output_list.append(layer_output)  # Collect outputs
            last_state_list.append(cur_state)  # Collect last states

        # Use the last output from the last layer
        output = layer_output_list[-1][:, -1, :, :, :]
        output = torch.mean(output, dim=(2, 3))  # Global average pooling
        output = self.fc(output)  # Fully connected layer
        return torch.sigmoid(output), last_state_list  # Sigmoid activation for binary classification

    def init_hidden(self, batch_size, shape):
        """
        Initialize hidden states for all layers.

        Args:
            batch_size (int): Size of the batch.
            shape (tuple): Spatial dimensions (height, width) of the input tensor.

        Returns:
            list: Initialized hidden states for each layer.
        """
        init_states = []
        for cell in self.cells:
            init_states.append(cell.init_hidden(batch_size, shape))
        return init_states

def pad_or_truncate_features(df_list, max_columns):
    """
    Pads or truncates the feature matrices to have the same number of columns.
    
    Args:
        df_list (list of pd.DataFrame): List of dataframes to be padded or truncated.
        max_columns (int): Maximum number of columns to pad or truncate to.

    Returns:
        list of pd.DataFrame: List of padded or truncated dataframes.
    """
    padded_dfs = []
    for df in df_list:
        if df.shape[1] < max_columns:
            while df.shape[1] < max_columns:
                padding_col = pd.Series([0] * df.shape[0])
                df.insert(0, df.shape[1], padding_col)  # Add a column to the start
                if df.shape[1] < max_columns:
                    df[df.shape[1]] = 0  # Add a column to the end
            padded_df = df
        elif df.shape[1] > max_columns:
            while df.shape[1] > max_columns:
                df = df.drop(df.columns[0], axis=1)  # Remove one column from the 0 index
                if df.shape[1] > max_columns:
                    df = df.drop(df.columns[-1], axis=1)  # Remove one column from the last index
            padded_df = df
        else:
            padded_df = df
        padded_dfs.append(padded_df)
    return padded_dfs

def prepare_data_and_labels(dataframes, max_columns):
    """
    Prepare data and labels by converting a list of dataframes to tensors.

    Args:
        dataframes (list of pd.DataFrame): List of dataframes containing the data.
        max_columns (int): Maximum number of columns for padding/truncating.

    Returns:
        tuple: A tuple containing the data tensor and the labels tensor.
    """
    dataframes = pad_or_truncate_features(dataframes, max_columns)  # Pad or truncate the features
    data_tensors = []
    labels = []

    for df in dataframes:
        if df.empty:
            continue  # Skip empty dataframes

        df = df.reset_index(drop=True)  # Reset index to ensure the first row has index 0
        label = df.iloc[0, :].values[0]  # Extract the label from the first row
        labels.append(label)

        features_df = df.drop(index=0)
        tensor = torch.tensor(features_df.values, dtype=torch.float32)  # Convert dataframe to tensor
        data_tensors.append(tensor.unsqueeze(0))  # Add channel dimension

    data_tensor = torch.stack(data_tensors, dim=0)  # Stack tensors to create a batch
    data_tensor = data_tensor.unsqueeze(1)  # Add sequence dimension
    labels_tensor = torch.tensor(labels, dtype=torch.float32).unsqueeze(1)  # Convert labels to tensor and add dimension

    return data_tensor, labels_tensor

def consensus_features(json_dir):
    """
    Calculate consensus features for each JSON file.

    Args:
        json_dir (str): Directory containing JSON files.

    Returns:
        list of pd.DataFrame: List of dataframes containing consensus level features.
    """
    warnings.filterwarnings("ignore")  # stops a warning from spamming your output
    sys.path.append('..//')  # path to directory holding feature_extraction
    json_file_names = os.listdir(json_dir)
    cons_features_list = []
    
    with alive_bar(len(json_file_names), title="Processing JSON files") as bar:
        for i in range(len(json_file_names)):
            json_file_path = os.path.join(json_dir, json_file_names[i])
            consensus_features = fe.feature_extraction(json_file_path)
            cons_features_list.append(consensus_features.T)
            bar()  # Update the progress bar
    return cons_features_list

def update_model_with_new_data(model, new_dataframes, optimizer, criterion, max_columns, epochs=10):
    """
    Updates the model with new training data.
    
    Args:
        model (nn.Module): The ConvLSTM model to be updated.
        new_dataframes (list of pd.DataFrame): List of new dataframes for training.
        optimizer (torch.optim.Optimizer): Optimizer for the model.
        criterion (nn.Module): Loss function.
        max_columns (int): Maximum number of columns for padding/truncating.
        epochs (int): Number of epochs for training.
    """
    new_data, new_labels = prepare_data_and_labels(new_dataframes, max_columns)  # Prepare new data and labels
    
    for epoch in range(epochs):
        model.train()  # Set model to training mode
        optimizer.zero_grad()  # Zero the gradients
        output, _ = model(new_data)  # Forward pass

        loss = criterion(output, new_labels)  # Compute loss
        loss.backward()  # Backward pass
        optimizer.step()  # Update parameters
        print(f'Update Epoch {epoch+1}, Loss: {loss.item()}')  # Print loss

def calculate_accuracy(output, labels):
    """
    Calculate the accuracy of the model.

    Args:
        output (torch.Tensor): The output predictions of the model.
        labels (torch.Tensor): The true labels.

    Returns:
        float: The accuracy as a percentage.
    """
    predictions = (output > 0.5).float()  # Convert probabilities to binary predictions
    correct = (predictions == labels).float().sum()  # Count correct predictions
    accuracy = correct / labels.numel()  # Calculate accuracy
    return accuracy.item() * 100

def export_model_torchscript(model, save_filename):
    """
    Export the trained model to TorchScript format.

    Args:
        model (nn.Module): The trained model.
        save_filename (str): The filename to save the model.
    """
    model.eval()  # Set the model to evaluation mode
    scripted_model = torch.jit.script(model)  # Convert the model to TorchScript
    torch.jit.save(scripted_model, save_filename)  # Save the TorchScript model
    print(f'Model exported to {save_filename}')

if __name__ == "__main__":
    working_dir = '/home/marchandlab/github/jay/capstone/XenoFind/xenofind_test/240521_convlstm_test'
    directories_list = setup.setup_directory_system(working_dir)
    json_dir = '/home/marchandlab/github/jay/capstone/datasets/pz_xem_lib_json_subsets'
    model_dir = directories_list[6]
    new_json_dir = None  # Define new_json_dir to avoid NameError

    # Generate the list of dataframes using consensus features
    dataframes = consensus_features(json_dir)
    max_columns = max(df.shape[1] for df in dataframes)

    # Prepare data and labels
    data, labels = prepare_data_and_labels(dataframes, max_columns)
    train_data, val_data, train_labels, val_labels = train_test_split(data, labels, test_size=0.2, random_state=42)

    print(f"Training data tensor shape: {train_data.shape}")
    print(f"Validation data tensor shape: {val_data.shape}")

    # Initialize the model
    model = ConvLSTM(input_channels=1, hidden_channels=16, kernel_size=3, num_layers=2)

    # Example training loop
    optimizer = torch.optim.Adam(model.parameters(), lr=0.001)  # Adam optimizer
    criterion = nn.BCELoss()  # Binary cross-entropy loss
    scheduler = CustomPlusCoolDownScheduler(
        optimizer, 
        initial_lr=0.001, 
        lr_scheduler_str="ReduceLROnPlateau", 
        lr_scheduler_kwargs={"mode": "min", "patience": 3, "factor": 0.1}, 
        epochs=50, 
        cool_down_epochs=5, 
        cool_down_lr=1e-5
    )

    best_accuracy = 0.0
    epochs_since_improvement = 0

    for epoch in range(50):
        model.train()  # Set the model to training mode
        optimizer.zero_grad()  # Zero the gradients
        output, _ = model(train_data)  # Forward pass

        print(f"Training output values (min, max): {output.min().item()}, {output.max().item()}")

        if torch.isnan(output).any():
            print("NaNs detected in training output. Stopping training.")
            break

        loss = criterion(output, train_labels)  # Compute loss
        loss.backward()  # Backward pass

        torch.nn.utils.clip_grad_norm_(model.parameters(), max_norm=1.0)  # Gradient clipping

        optimizer.step()  # Update parameters

        train_accuracy = calculate_accuracy(output, train_labels)  # Calculate training accuracy
        print(f'Epoch {epoch+1}, Training Loss: {loss.item()}, Training Accuracy: {train_accuracy:.2f}%')

        model.eval()  # Set the model to evaluation mode
        with torch.no_grad():  # Disable gradient calculation
            val_output, _ = model(val_data)  # Forward pass
            val_loss = criterion(val_output, val_labels)  # Compute validation loss

            print(f"Validation output values (min, max): {val_output.min().item()}, {val_output.max().item()}")

            val_accuracy = calculate_accuracy(val_output, val_labels)  # Calculate validation accuracy
            print(f'Epoch {epoch+1}, Validation Loss: {val_loss.item()}, Validation Accuracy: {val_accuracy:.2f}%')

        scheduler.step(val_loss)  # Update learning rate based on validation loss

        if val_accuracy > best_accuracy:
            best_accuracy = val_accuracy
            epochs_since_improvement = 0
        else:
            epochs_since_improvement += 1

        if epochs_since_improvement >= 5:
            print("No improvement in validation accuracy after 5 epochs. Stopping training.")
            break

    # Save the trained model
    model_path = os.path.join(model_dir, 'convlstm_model.pt')
    torch.save(model.state_dict(), model_path)
    print(f'Model saved to {model_path}')

    # If there is a directory for new training data, update the model
    if new_json_dir and os.path.exists(new_json_dir):
        new_dataframes = consensus_features(new_json_dir)
        update_model_with_new_data(model, new_dataframes, optimizer, criterion, max_columns, epochs=10)
        torch.save(model.state_dict(), model_path)
        print(f'Updated model saved to {model_path}')

    # Export the trained model to TorchScript format
    export_model_torchscript(model, os.path.join(model_dir, 'convlstm_model_scripted.pt'))

