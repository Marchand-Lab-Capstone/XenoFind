# XenoFind: A de novo XNA consensus formation & detection pipeline for XNA aptamers 

## About 

XenoFind is a modified/synthetic base detection pipeline stemming from Chem E 546/547 in collaboration with the Marchand Group at the University of Washington. XenoFind was motivated by the lack of de novo basecalling pipelines for unnatural base-pairing xenonucleic acid (ubpXNAs or XNAs)  and aims to bridge the gap between canonical & modified basecalling in the context of XNA containing aptamer therapeutics. The XenoFind pipeline has two main functionalities: generating accurate consensus sequences semi-de novo for aptamer sequences and detecting the presence of XNAs in the consensus sequences.  

## Libraries

Within lib on the Main branch, there are currently four directories.


consensus_formation - scripts for consensus sequence generation (FASTA output)

model_gen - scripts for XenoFind PyTorch model generation (model output)

xna_finder - scripts for XNA detection within a consensus FASTA (xFASTA output)

xf_prelims - contains original drafts of the current tools and methods for xenofind. May or may not be functional. 

## Dependencies 
XenoFind requires ONT tools (pod5, Dorado), several bioinformatics packages, and various python packages. A full list of dependencies can be found in the xenofind-env.yml document, with the exception of dorado which requires manual installation through ONT. XenoFind was built and tested on Ubuntu 20.04.To use conda for installing dependencies, simply load xenofind-env.yml into a new environment with the following command: 

    conda env create -f xenofind-env.yml
To enter the XenoFind conda environment, then use:

	conda activate xenofind

XenoFind requires users to have the basecaller [dorado](https://github.com/nanoporetech/dorado) installed in a user's home directory. 
This can be achieved with the folllowing commands: 

	cd 
 	wget -qO- https://github.com/nanoporetech/dorado | tar xvz


## XenoFind command groups 

Xenofind has three subcommands to call: Consensus Formation, Find, and Model Generation. 

### Consensus Formation
Consensus Formation primarily uses VSEARCH to create accurate consensus sequences from the dataset structure above. These consensus sequences are decoupled for both the forward and reverse strands to reduce end specific nanopore sequencing error and are labelled accordingly. The consensus sequences are outputted in the FASTA file format in the desired working directory as "final_merged_consensus.fasta". This FASTA file will be a mandatory input in "Find". 


Consensus Formation requires the following inputs: a desired working directory, a POD5 or FAST5 dataset, and a reference FASTA file containing the following structure. 

    Barcode (known) - Region of Interest (NNNNN in reference) - Barcode (known)
Consensus Formation is a mandatory preprocessing step for data with unidentified regions (the region of interest above). 
To run consensus formation, call the following command in terminal 

    python xenofind.py consensus -w [desired working directory] -f [pod5 or fast5 directory] -r [NNN reference fasta]

### XNA Find
XNA Find performs XNA detection by using either the included XenoFind model or one generated by the user. XNA Find creates 'consensus-level features' from the inputted dataset using both raw and sequence space data. Consensus-level features are features generated by grouping reads based on alignment and generates features based on the groups of aligned reads (e.g., the average current signal, average residence time, Shannon Entropy, etc.). This differs significantly from other detection pipelines, which takes on a per-read level (e.g., the raw signal from a single read).

XNA Find requires the following inputs: a desired working directory, the POD5 or FAST5 dataset used for consensus formation, the consensus fasta generated from consensus formation, and a XenoFind model (provided or can be self generated). To run XNA Find, call the following in terminal.

	python xenofind.py find -w [working directory] -f [pod5 or fast5 directory] -r [consensus fasta from consensus formation] -m [XenoFind model]

XNA Find will output an updated consensus sequence reference file in the xFASTA format to be used with an XNA basecaller such as Xemora or Xenomorph. For more information on the xFASTA format, a section describing the format is provided in the Model Generation section below. 

### Model Generation 
XenoFind models are PyTorch based models that are trained using the aforementioned consensus-level features. XenoFind has a pre-trained model under the models directory but users can also generate their own models if desired.

Model Generation requires the following inputs: a desired working directory, a POD5 or FAST5 dataset, and a reference FASTA file containing the following structure.

	...NNNXNNN...

Where the '...' and 'NNN' regions are ATGC bases and X is an XNA present in lib/xf_params. New XNA base pairs can be added to xf_params by updating the xna_base_pairs and confounding_pairs variables. Model Generation will generate xFASTA files for both the forward and reverse strands. 

#### xFASTA format (adapted from [Xenomorph] https://github.com/xenobiolab/xenomorph)
Many tools used to read and manipulate nucleic acid sequences are not built to handle non ATGC bases. The xFASTA file format (.fa) stores XNA positional information in the header of each sequence. xFASTA files can be generated from standard Fasta files that contain non-ATGC bases (e.g. BSPZJVKX) in the sequence. xFASTA files are automatically generated in the standard xenomorph preprocessing workflow. The fasta2x command is provided for utility, but generally not required. Note that XNA bases in the input sequence will be replaced for a standard ATGC. Signal matching to sequence is highly sensitive what base is chosen as the replacement base. As default, XNA bases are replaced as followed: B>A, S>A, P>G, Z>C. Base substitution settings can be modified in lib/xf_params.py by changing the paired base in the confounding_base variable. 


        Standard FASTA format with XNA in sequence
                >header_for_a_read
                ATGGCAACAGGATGAPAAGGACGTA

        xFASTA format with XNA information stored in header. (P replaced with G in sequence)
                >header_for_a_read+X_POS[P:18]
                ATGGCAACAGGATGAGAAGGACGTA
