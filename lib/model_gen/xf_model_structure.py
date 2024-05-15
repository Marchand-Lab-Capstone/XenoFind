import torch
import torch.nn as nn

class ConvLSTMCell(nn.Module):
    """
    Implements a Convolutional LSTM cell.

    This cell is a basic building block for a Convolutional LSTM layer, which processes
    input data incorporating both spatial and temporal dimensions.

    Attributes:
        input_channels (int): Number of input channels.
        hidden_channels (int): Number of output channels from the hidden state.
        kernel_size (tuple): Size of the convolutional kernel.
        padding (int): Padding added to the input for the convolution operation.
        bias (bool): If True, adds a learnable bias to the output.
        conv (nn.Module): Convolutional layer that performs combined input and hidden state convolutions.
    """

    def __init__(self, input_channels, hidden_channels, kernel_size, bias=True):
        super(ConvLSTMCell, self).__init__()
        self.input_channels = input_channels
        self.hidden_channels = hidden_channels
        self.kernel_size = kernel_size
        self.padding = kernel_size // 2
        self.bias = bias

        # The convolution layer that combines input and previous states
        self.conv = nn.Conv2d(in_channels=input_channels + hidden_channels,
                              out_channels=4 * hidden_channels,
                              kernel_size=self.kernel_size,
                              padding=self.padding,
                              bias=self.bias)

    def forward(self, input_tensor, cur_state):
        """
        Forward pass of the ConvLSTMCell.

        Args:
            input_tensor (torch.Tensor): The input tensor at a time step with shape [batch, channels, height, width].
            cur_state (tuple): A tuple (h, c) containing the current hidden state and the cell state.

        Returns:
            tuple: A tuple (h_next, c_next) containing the next hidden state and the cell state.
        """
        h_cur, c_cur = cur_state
        combined = torch.cat([input_tensor, h_cur], dim=1)  # combine along channel dimension
        combined_conv = self.conv(combined)

        # Split the combined convolution output into four parts for gates
        cc_i, cc_f, cc_o, cc_g = torch.split(combined_conv, self.hidden_channels, dim=1)
        i = torch.sigmoid(cc_i)  # input gate
        f = torch.sigmoid(cc_f)  # forget gate
        o = torch.sigmoid(cc_o)  # output gate
        g = torch.tanh(cc_g)  # new candidate

        c_next = f * c_cur + i * g  # cell state update
        h_next = o * torch.tanh(c_next)  # hidden state update

        return h_next, c_next

    def init_hidden(self, batch_size, shape):
        """
        Initialize hidden and cell states.

        Args:
            batch_size (int): The size of the batch.
            shape (tuple): The spatial dimensions (height, width) of the input tensor.

        Returns:
            tuple: Initialized hidden and cell states, each with dimensions [batch, channels, height, width].
        """
        height, width = shape
        return (torch.zeros(batch_size, self.hidden_channels, height, width, device=self.conv.weight.device),
                torch.zeros(batch_size, self.hidden_channels, height, width, device=self.conv.weight.device))


class ConvLSTM(nn.Module):
    """
    A Convolutional LSTM Network layer capable of processing sequences of images or spatial data.

    Attributes:
        input_channels (int): Number of input channels.
        hidden_channels (int): Number of channels in hidden states.
        kernel_size (tuple): Convolutional kernel size.
        num_layers (int): Number of ConvLSTM layers.
        bias (bool): If True, includes a bias term in convolution operations.
        cells (nn.ModuleList): List of ConvLSTM cells.
    """

    def __init__(self, input_channels, hidden_channels, kernel_size, num_layers, bias=True):
        super(ConvLSTM, self).__init__()
        self.input_channels = input_channels
        self.hidden_channels = hidden_channels
        self.kernel_size = kernel_size
        self.num_layers = num_layers
        self.bias = bias

        # Creating a list of ConvLSTM cells for each layer
        cells = []
        for i in range(self.num_layers):
            in_dim = self.input_channels if i == 0 else self.hidden_channels
            cells.append(ConvLSTMCell(in_dim, self.hidden_channels, self.kernel_size, self.bias))

        self.cells = nn.ModuleList(cells)

    def forward(self, input_tensor, hidden_state=None):
        """
        Forward pass of the ConvLSTM layer.

        Args:
            input_tensor (torch.Tensor): Input tensor of shape [batch, sequence_length, channels, height, width].
            hidden_state (list, optional): Initial hidden states for each layer.

        Returns:
            tuple: Outputs and last states from all layers.
        """
        layer_output_list = []
        last_state_list = []

        seq_len = input_tensor.size(1)
        cur_layer_input = input_tensor

        for layer_idx in range(self.num_layers):
            # Initialize hidden states if not provided
            h, c = self.cells[layer_idx].init_hidden(input_tensor.size(0), (input_tensor.size(3), input_tensor.size(4))) if hidden_state is None else hidden_state[layer_idx]
            output_inner = []
            for t in range(seq_len):
                h, c = self.cells[layer_idx](cur_layer_input[:, t, :, :, :], (h, c))
                output_inner.append(h)

            layer_output = torch.stack(output_inner, dim=1)
            cur_layer_input = layer_output

            layer_output_list.append(layer_output)
            last_state_list.append((h, c))

        return layer_output_list, last_state_list

    def init_hidden(self, batch_size, shape):
        """
        Initialize hidden states for all layers.

        Args:
            batch_size (int): The size of the batch.
            shape (tuple): The spatial dimensions (height, width) of the input tensor.

        Returns:
            list: A list of hidden states for each layer.
        """
        init_states = []
        for i in range(self.num_layers):
            init_states.append(self.cells[i].init_hidden(batch_size, shape))
        return init_states

