# XenoFind 

## Overview 

XenoFind is a modified/synthetic base detection pipeline stemming from Chem E 546/547 from the University of Washington. The purpose of this pipeline is to identify regions within DNA where non-standard bases may be present. This information will 

### Tentative Objectives 

Winter quarter:

Primary Objective: Benchmark and compare various modified base detection strategies

• Literature review for current state-of-the-field modified base detection strategies that can be applied to XNA world

• Develop pipelines to extract features from nanopore signal data formats (such as .FAST5, .POD5)

• Determine benchmarking metrics that will be used to evaluate and compare different strategies. This will be instrumental for the approach to measure performance and progress throughout the project.


Spring quarter:

Primary Objective: Create a modified basecalling pipeline, implementing one or multiple of the strategies for de novo modified base detection from winter quarter. (With or without ML)

• Tier 1 – Develop processing pipeline that will take raw signal data and output probability that a given signal does not belong to DNA.

o Low accuracy, Low precision

•Tier 2 – Develop an end-to-end basecalling pipeline (ML or not) that can be used accurately identify regions of DNA that contain an XNA.

o High accuracy, low precision

• Tier 3 – Develop an end-to-end basecalling pipeline that uses a (PyTorch model, such as a LSTM RNN or other) to identify position and identity of XNA within a sequence from a raw data input.

o High accuracy, high precision
