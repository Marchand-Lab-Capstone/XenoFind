# XenoFind: A de novo XNA detection pipeline

## About 

XenoFind is a modified/synthetic base detection pipeline stemming from Chem E 546/547 in collaboration with the Marchand Group at the University of Washington. XenoFind was motivated by the lack of de novo basecalling pipelines for unnatural base-pairing xenonucleic acid (ubpXNAs or XNAs)  and aims to bridge the gap between canonical & modified basecalling. The purpose of this pipeline is to identify regions within DNA where non-standard bases may be present and output that information in bioinformatic relevant file formats to be used in downstream analysis. 

## Dependencies 

## XenoFind command groups 

Xenofind has two primary subcommands to call: consensus or low_qual. 

Consensus requires the following inputs: a desired working directory, a pod5 or fast5 dataset, and a temporary reference file containing the following structure. 

    Barcode (known) - Region of Interest (NNNNN in reference) - Barcode (known)
Consensus formation is a mandatory preprocessing step for data with unidentified regions (the region of interest above). 
To run consensus formation, call the following command in terminal 

    python xenofind.py consensus -w [working directory] -f [pod5 or fast5 directory] -r [temporary reference fasta]

low_qual performs XNA detection by analyzing and performing statistical tests on a per base level and per consensus level.Regions that have several areas with low quality cores compared to the mean quality are logged as potential XNA positions. To call low quality score analysis, run the following command. 

	python xenofind.py low_qual -w [working directory] -f [pod5 or fast5 directory] -r [consensus fasta from consensus]


### Winter Objectives 

Overall Objective: Benchmark and compare various modified base detection strategies

• Literature review for current state-of-the-field modified base detection strategies that can be applied to XNA world

• Develop pipelines to extract features from nanopore signal data formats (such as .FAST5, .POD5)

• Create a consensus reference file using the nanopore sequencing data starting with a placeholder reference, only containing a known barcode region and an unknown region of interest. 

• Detect candidate XNA position using extracted quality scores and statistical methods. Identify these positions in the header of the consensus fasta file. 
