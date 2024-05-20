import torch
import torch.nn as nn
import pandas as pd
import numpy as np

class ConvLSTMCell(nn.Module):
    """
    Convolutional LSTM cell.

    Args:
        input_channels (int): Number of input channels.
        hidden_channels (int): Number of hidden channels.
        kernel_size (int): Size of the convolutional kernel.
        bias (bool): Whether to include a bias term in the convolution.

    Methods:
        forward(input_tensor, cur_state): Computes the forward pass of the cell.
        init_hidden(batch_size, shape): Initializes the hidden and cell states.
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
            cur_state (tuple): Current hidden state and cell state (h_cur, c_cur).

        Returns:
            tuple: Next hidden state and cell state (h_next, c_next).
        """
        h_cur, c_cur = cur_state  # Unpack current states
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

        return h_next, c_next

    def init_hidden(self, batch_size, shape):
        """
        Initialize hidden and cell states with zeros.

        Args:
            batch_size (int): Size of the batch.
            shape (tuple): Spatial dimensions (height, width) of the input tensor.

        Returns:
            tuple: Initialized hidden and cell states (h, c).
        """
        height, width = shape
        return (torch.zeros(batch_size, self.hidden_channels, height, width, device=self.conv.weight.device),
                torch.zeros(batch_size, self.hidden_channels, height, width, device=self.conv.weight.device))


class ConvLSTM(nn.Module):
    """
    Convolutional LSTM Network.

    Args:
        input_channels (int): Number of input channels.
        hidden_channels (int): Number of hidden channels.
        kernel_size (int): Size of the convolutional kernel.
        num_layers (int): Number of ConvLSTM layers.
        bias (bool): Whether to include a bias term in the convolutions.

    Methods:
        forward(input_tensor, hidden_state): Computes the forward pass of the network.
        init_hidden(batch_size, shape): Initializes hidden states for all layers.
    """
    
    def __init__(self, input_channels, hidden_channels, kernel_size, num_layers, bias=True):
        super(ConvLSTM, self).__init__()
        self.input_channels = input_channels
        self.hidden_channels = hidden_channels
        self.kernel_size = kernel_size
        self.num_layers = num_layers
        self.bias = bias

        # Create a list of ConvLSTM cells for each layer
        cells = []
        for i in range(self.num_layers):
            in_dim = self.input_channels if i == 0 else self.hidden_channels
            cells.append(ConvLSTMCell(in_dim, self.hidden_channels, self.kernel_size, self.bias))

        self.cells = nn.ModuleList(cells)  # Register the cells as a module list
        self.fc = nn.Linear(hidden_channels, 1)  # Fully connected layer for binary classification

    def forward(self, input_tensor, hidden_state=None):
        """
        Forward pass of the ConvLSTM layer.

        Args:
            input_tensor (torch.Tensor): Input tensor (batch, sequence_length, channels, height, width).
            hidden_state (list, optional): Initial hidden states for each layer.

        Returns:
            tuple: Output of the last layer and the last hidden states.
        """
        layer_output_list = []
        last_state_list = []

        seq_len = input_tensor.size(1)  # Sequence length
        cur_layer_input = input_tensor  # Initialize current layer input

        # Iterate through each layer
        for layer_idx in range(self.num_layers):
            # Initialize hidden states if not provided
            h, c = self.cells[layer_idx].init_hidden(input_tensor.size(0), (input_tensor.size(3), input_tensor.size(4))) if hidden_state is None else hidden_state[layer_idx]
            output_inner = []

            # Iterate through each time step
            for t in range(seq_len):
                h, c = self.cells[layer_idx](cur_layer_input[:, t, :, :, :], (h, c))
                output_inner.append(h)

            layer_output = torch.stack(output_inner, dim=1)  # Stack outputs
            cur_layer_input = layer_output  # Set current layer input to output of this layer

            layer_output_list.append(layer_output)  # Collect outputs
            last_state_list.append((h, c))  # Collect last states

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
        for i in range(self.num_layers):
            init_states.append(self.cells[i].init_hidden(batch_size, shape))
        return init_states


def prepare_data(dataframes):
    """
    Prepare data by converting a list of dataframes to a tensor.

    Args:
        dataframes (list of pd.DataFrame): List of dataframes containing the data.

    Returns:
        torch.Tensor: Tensor of shape (batch_size, sequence_length, channels, height, width).
    """
    tensors = []
    for df in dataframes:
        tensor = torch.tensor(df.values, dtype=torch.float32)  # Convert dataframe to tensor
        tensors.append(tensor.unsqueeze(0))  # Add channel dimension

    return torch.stack(tensors, dim=0)  # Stack tensors to create a batch


# Create a list of dataframes as an example
dataframes = [pd.DataFrame(np.random.rand(53, 133)) for _ in range(10)]
data = prepare_data(dataframes)  # Prepare data

# Initialize the model
model = ConvLSTM(input_channels=1, hidden_channels=16, kernel_size=3, num_layers=2)

# Forward pass
output, _ = model(data)
print(output)

# Example training loop
optimizer = torch.optim.Adam(model.parameters(), lr=0.001)  # Adam optimizer
criterion = nn.BCELoss()  # Binary cross-entropy loss

# Assuming you have labels for each dataframe in your batch
labels = torch.tensor([0, 1, 0, 1, 0, 1, 0, 1, 0, 1], dtype=torch.float32).unsqueeze(1)  # Binary labels

for epoch in range(10):
    optimizer.zero_grad()  # Zero the gradients
    output, _ = model(data)  # Forward pass
    loss = criterion(output, labels)  # Compute loss
    loss.backward()  # Backward pass
    optimizer.step()  # Update parameters
    print(f'Epoch {epoch+1}, Loss: {loss.item()}')  # Print loss

