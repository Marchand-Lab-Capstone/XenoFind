# Technical Specifications
## 1. Introduction 

### 1.1 Purpose 
This document will provide a comprehensive overview of XenoFind, a modified base detection pipeline. In this document, the scope of the project, overall project description, hardware & system requirements, and tentative timeline for Winter 2024. 

### 1.2 Project Description
XenoFind is Python based pipeline that will be used for the detection of 

### 1.3 Definitions, Acronyms, and Abbreviations
### 1.3.1 Modified base - 
### 1.3.2 Xenonucleic Acid (XNA) - Chemically synthesized non-standard nucleobase
	1.3.2 
## 2. Requirements 

## 2.1 Hardware Requirements 
	### 2.1.1 Storage
Filesize should be minimal enough to not exceed that of a dataset it is trained on (no more than, say, one gigabyte)
	### 2.1.2 Memory
		Dependent on system, should run on at least 4 gigabytes of RAM
	### 2.1.3 CPU
		32 bit system bare minimum, multithreading is a stretch goal, so MVP is one core.
	### 2.1.4 GPU
GPU resources are recommended but not required. Recommended GPU: Nvidia RTX 3060, 12 GB VM. 
### 2.2 Software Requirements 
### 2.2.1 Operating System
Ubuntu 18.04 and 20.04
### 2.2.2 Python Version
	Python version 3.8.x
### 2.2.3 Software packages 
	<TODO:List Required Package Dependancies here>
	Pytorch
	Tensorflow
	etc….
### 2.2.4  External Programs 
Dorado (https://github.com/nanoporetech/dorado) 

## 3. Software Architecture 

### 3.1 Software Design 

### 3.2 Software Components
	3.2.1 xf_pipe.py 
		Python script that takes in user inputs such as datasets, output directories, desired method, etc. This script will be run if using terminal for inputs is undesired.  
### 3.3 Data Flow Diagram

## 7. Milestones and Timeline

## 8. Appendices 

### 8.1 References 
