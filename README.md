# XenoFind: A de novo XNA detection pipeline

## About 

XenoFind is a modified/synthetic base detection pipeline stemming from Chem E 546/547 in collaboration with the Marchand Group at the University of Washington. XenoFind was motivated by the lack of de novo unnatural base-pairing xenonucleic acid (ubpXNAs or XNAs) basecalling pipelines and aims to bridge the gap between canonical & modified basecalling. The purpose of this pipeline is to identify regions within DNA where non-standard bases may be present and output that information in bioinformatic relevant file formats to be used in downstream analysis. 

## Dependencies 

## XenoFind command groups 

### Tentative Objectives 

Winter quarter:

Overall Objective: Benchmark and compare various modified base detection strategies

• Literature review for current state-of-the-field modified base detection strategies that can be applied to XNA world

• Develop pipelines to extract features from nanopore signal data formats (such as .FAST5, .POD5)

• Determine benchmarking metrics that will be used to evaluate and compare different strategies. This will be instrumental for the approach to measure performance and progress throughout the project.


Spring quarter:

Primary Objective: Create a modified basecalling pipeline, implementing one or multiple of the strategies for de novo modified base detection from winter quarter. (With or without ML)

• Tier 1 – Develop processing pipeline that will take raw signal data and output probability that a given signal does not belong to DNA. Low accuracy, Low precision

•Tier 2 – Develop an end-to-end basecalling pipeline (ML or not) that can be used accurately identify regions of DNA that contain an XNA. High accuracy, low precision

• Tier 3 – Develop an end-to-end basecalling pipeline that uses a (PyTorch model, such as a LSTM RNN or other) to identify position and identity of XNA within a sequence from a raw data input. High accuracy, high precision
