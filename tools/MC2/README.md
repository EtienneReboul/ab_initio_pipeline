# MoRFchibi 2.0
Molecular Recognition Features (MoRFs) are a subset of protein-binding IDR that undergo a disorder-to-order transition upon binding. MoRFchibi 2.0 is a specialized prediction tool designed to identify the locations of MoRFs within protein sequences. Output scores are generated using an ensemble of convolutional neural network logistic regression models, followed by a reverse Bayes Rule to adjust for priors in the training data. These scores reflect MoRF probabilities normalized for the priors in the training data, making them individually interpretable and compatible with other tools utilizing the same scoring framework, such as IPA.

Training and test data are available at: https://doi.org/10.5281/zenodo.20299225 

#### Please reference the following publications:

Malhis N, Gsponer J. "Predicting molecular recognition features in protein sequences with MoRFchibi 2.0" *bioRxiv* (2025). [doi.org/10.1101/2025.01.31.635962.] (https://doi.org/10.1101/2025.01.31.635962) 

## Minimum Hardware Requirements

OS: Linux (tested on Ubuntu and CentOS7).

RAM: 8 GB minimum, 16 GB recommended.

CPU: Multicore with 4+ cores recommended.

## To install:
```bash
# Clone MC2
git clone https://github.com/NawarMalhis/MC2.git
# Clone the "annotated fasta format" library:	
git clone https://github.com/NawarMalhis/AFF.git
# Change directory:	
cd MC2
# Update the path to the AFF (annotated fasta format) folder in param.py
aff_path = '/xxx/xxx/AFF/'
# Create an mc2 environment:
conda env create -f mc2.yml
```

## To run:
```bash
# Activate the mc2 environment:
conda activate mc2_cpu
```
Input sequences are in the input.fasta, output will be in ./output/

Run MoRFchibi 2.0:
```bash
python3 MC2.py
```
## HTML Interface
[https://mc2.msl.ubc.ca]
