## FASTA Sequence Analyzer

A Python-based bioinformatics tool that reads DNA sequences from a FASTA file and generates a full statistical report that includes per-sequence metrics, dataset-level summary statistics, and visualizations.

# Features
- Reads DNA sequences from a multi-FASTA file
- Calculates sequence length
- Calculates GC content and AT content
- Counts A, T, G, and C nucleotides per sequence
- Validates that sequences contain only standard DNA bases
- Detects duplicate sequence IDs
- Calculates dataset-level statistics (mean, min, max, standard deviation)    using NumPy
- Saves per-sequence results to a CSV file
- Saves a dataset-level summary to a text file
- Generates GC content, sequence length, and nucleotide composition plots

<<<<<<< HEAD
## Project Structure
                  FASTA SEQUENCE ANALYZER
                           │
                    sequence.fasta
                           │
                           ▼
                    ┌─────────────┐
                    │  analyze.py │
                    └──────┬──────┘
                           │
             ┌─────────────┼─────────────┐
             ▼             ▼             ▼
        GC Content    Sequence Length   Nucleotide
                                         Composition
             │             │             │
             └─────────────┼─────────────┘
                           ▼
                       RESULTS
                           │
             ┌─────────────┼─────────────┐
             ▼             ▼             ▼
           CSV           Summary        Plots
## Technologies Used
=======
# Workflow
              FASTA SEQUENCE ANALYZER
                       │
                sequence.fasta
                       │
                       ▼
                ┌─────────────┐
                │  analyze.py │
                └──────┬──────┘
                       │
         ┌─────────────┼─────────────┐
         ▼             ▼             ▼
    GC Content    Sequence Length   Nucleotide
                                     Composition
         │             │             │
         └─────────────┼─────────────┘
                       ▼
                   RESULTS
                       │
         ┌─────────────┼─────────────┐
         ▼             ▼             ▼
       CSV           Summary        Plots

# Project Structure
fasta-sequence-analyzer/
├── analyze.py
├── requirements.txt
├── README.md
├── data/
│   └── sequence.fasta
└── results/            (created automatically when you run the script)
    ├── sequence_analysis.csv
    ├── summary.txt
    ├── gc_content.png
    ├── sequence_lengths.png
    └── nucleotide_composition.png

# Technologies Used
Python 3.12.4
Biopython
NumPy
Matplotlib
>>>>>>> b76652e (Documents edited)

# Installation
Clone the repository:
   git clone https://github.com/anchal-bio/fasta-sequence-analyzer.git
   cd fasta-sequence-analyzer

Install the dependencies:
   pip install -r requirements.txt

# How to Run

Place your FASTA file at data/sequence.fasta, then run:

python analyze.py

The results/ folder is created automatically if it doesn't already exist.

# Input

A multi-FASTA file at data/sequence.fasta.

# Output

All results are saved inside the results/ folder.

- sequence_analysis.csv — per-sequence data:

* Sequence ID
* Sequence length
* GC content
* AT content
* A, T, G, and C counts

- summary.txt — dataset-level statistics: total sequence count, average/min/max/standard deviation of length and GC content, shortest and longest sequence, duplicate IDs, and invalid sequences.

- Plots:

* gc_content.png — GC content per sequence
* sequence_lengths.png — length per sequence
* nucleotide_composition.png — A/T/G/C composition per sequence

## License

Distributed under the MIT License.

## Author

<<<<<<< HEAD
The program generates:

- GC content graph
- Sequence length graph
- Nucleotide composition graph
=======
- **Anchal Joshi** - [GitHub Profile](https://github.com/anchal-bio)
>>>>>>> b76652e (Documents edited)
