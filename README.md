# FASTA Sequence Analyzer

A Python-based bioinformatics project for analyzing DNA sequences
stored in FASTA format.

## Features

- Reads DNA sequences from a FASTA file
- Calculates sequence length
- Calculates GC content
- Calculates AT content
- Counts A, T, G and C nucleotides
- Validates DNA sequences
- Detects duplicate sequence IDs
- Calculates sequence statistics using NumPy
- Saves results in CSV format
- Saves analysis summary in a text file
- Generates GC content graph
- Generates sequence length graph
- Generates nucleotide composition graph

## Project Structure
fasta-sequence-analyzer/
│
├── data/
│   └── sequence.fasta
│
├── results/
│   ├── sequence_analysis.csv
│   ├── summary.txt
│   ├── gc_content.png
│   ├── sequence_lengths.png
│   └── nucleotide_composition.png
│
├── analyze.py
└── README.md
## Technologies Used

- Python
- Biopython
- NumPy
- Matplotlib
- CSV

## How to Run

Run the following command from the project folder:

python analyze.py

## Input

The program takes a multi-FASTA file as input:

data/sequence.fasta

## Output

The analysis results are saved inside the `results` folder.

### CSV

sequence_analysis.csv

Contains:

- Sequence ID
- Sequence length
- GC content
- AT content
- A, T, G and C counts

### Summary

summary.txt

Contains overall statistics of the FASTA dataset.

### Graphs

The program generates:

- GC content graph
- Sequence length graph
- Nucleotide composition graph
