"""
FASTA Sequence Analyzer

Reads a multi-FASTA file and computes per-sequence and dataset-level
statistics: length, GC/AT content, nucleotide composition, and
sequence validity. Results are saved as a CSV, a text summary, and plots.
"""

import os
import csv
from collections import Counter

import numpy as np
import matplotlib.pyplot as plt
from Bio import SeqIO


INPUT_FILE = "data/sequence.fasta"
OUTPUT_DIR = "results"


def calculate_gc_content(sequence):
    """Return GC content (%) of a sequence. Returns 0 for empty sequences."""
    if len(sequence) == 0:
        return 0.0
    gc_count = sequence.count("G") + sequence.count("C")
    return (gc_count / len(sequence)) * 100


def calculate_at_content(sequence):
    """Return AT content (%) of a sequence. Returns 0 for empty sequences."""
    if len(sequence) == 0:
        return 0.0
    at_count = sequence.count("A") + sequence.count("T")
    return (at_count / len(sequence)) * 100


def count_nucleotides(sequence):
    """Return counts of A, T, G, C in a sequence."""
    return (
        sequence.count("A"),
        sequence.count("T"),
        sequence.count("G"),
        sequence.count("C"),
    )


def validate_sequence(sequence):
    """Check whether a sequence contains only standard DNA bases (A, T, G, C)."""
    valid_bases = set("ATGC")
    return set(sequence).issubset(valid_bases)


def read_sequences(fasta_path):
    """
    Parse the FASTA file and return a list of dicts, one per sequence,
    containing its id, sequence, length, GC/AT content, nucleotide
    counts, and validity flag.
    """
    records = []

    for record in SeqIO.parse(fasta_path, "fasta"):
        # Uppercase so lowercase/soft-masked bases are still counted correctly
        sequence = str(record.seq).upper()

        a, t, g, c = count_nucleotides(sequence)

        records.append({
            "id": record.id,
            "length": len(sequence),
            "gc_content": calculate_gc_content(sequence),
            "at_content": calculate_at_content(sequence),
            "a": a, "t": t, "g": g, "c": c,
            "is_valid": validate_sequence(sequence),
        })

    return records


def find_duplicate_ids(records):
    """Return a list of sequence IDs that appear more than once."""
    id_counts = Counter(r["id"] for r in records)
    return [seq_id for seq_id, count in id_counts.items() if count > 1]


def write_csv(records, output_path):
    """Write per-sequence results to a CSV file."""
    with open(output_path, "w", newline="") as file:
        writer = csv.writer(file)
        writer.writerow(
            ["Sequence ID", "Length", "GC Content", "AT Content", "A", "T", "G", "C"]
        )

        for r in records:
            writer.writerow([
                r["id"], r["length"],
                round(r["gc_content"], 2), round(r["at_content"], 2),
                r["a"], r["t"], r["g"], r["c"],
            ])


def compute_statistics(records):
    """Compute dataset-level length and GC-content statistics using NumPy."""
    lengths = np.array([r["length"] for r in records])
    gc_contents = np.array([r["gc_content"] for r in records])

    shortest = min(records, key=lambda r: r["length"])
    longest = max(records, key=lambda r: r["length"])

    return {
        "count": len(records),
        "avg_length": lengths.mean(),
        "min_length": lengths.min(),
        "max_length": lengths.max(),
        "std_length": lengths.std(),
        "shortest_id": shortest["id"],
        "longest_id": longest["id"],
        "avg_gc": gc_contents.mean(),
        "min_gc": gc_contents.min(),
        "max_gc": gc_contents.max(),
        "std_gc": gc_contents.std(),
    }


def print_summary(stats, duplicate_ids, invalid_ids):
    """Print the dataset summary to the console."""
    print("Total Sequences:", stats["count"])
    print("Average Length:", round(stats["avg_length"], 2))
    print("Minimum Length:", stats["min_length"])
    print("Maximum Length:", stats["max_length"])
    print("Standard Deviation of Length:", round(stats["std_length"], 2))
    print("Shortest Sequence:", stats["shortest_id"], "-", stats["min_length"])
    print("Longest Sequence:", stats["longest_id"], "-", stats["max_length"])
    print("Average GC Content:", round(stats["avg_gc"], 2), "%")
    print("Minimum GC Content:", round(stats["min_gc"], 2), "%")
    print("Maximum GC Content:", round(stats["max_gc"], 2), "%")
    print("Standard Deviation of GC Content:", round(stats["std_gc"], 2), "%")
    print("Duplicate IDs:", duplicate_ids if duplicate_ids else "None")
    print("Invalid Sequences:", invalid_ids if invalid_ids else "None")


def save_summary(stats, duplicate_ids, invalid_ids, output_path):
    """Save the dataset summary to a text file."""
    with open(output_path, "w") as file:
        file.write("FASTA SEQUENCE ANALYSIS SUMMARY\n")
        file.write("===============================\n\n")
        file.write(f"Total Sequences: {stats['count']}\n")
        file.write(f"Average Length: {stats['avg_length']:.2f}\n")
        file.write(f"Minimum Length: {stats['min_length']}\n")
        file.write(f"Maximum Length: {stats['max_length']}\n")
        file.write(f"Standard Deviation of Length: {stats['std_length']:.2f}\n\n")
        file.write(f"Shortest Sequence: {stats['shortest_id']} - {stats['min_length']}\n")
        file.write(f"Longest Sequence: {stats['longest_id']} - {stats['max_length']}\n\n")
        file.write(f"Average GC Content: {stats['avg_gc']:.2f}%\n")
        file.write(f"Minimum GC Content: {stats['min_gc']:.2f}%\n")
        file.write(f"Maximum GC Content: {stats['max_gc']:.2f}%\n")
        file.write(f"Standard Deviation of GC Content: {stats['std_gc']:.2f}%\n\n")
        file.write(f"Duplicate IDs: {duplicate_ids if duplicate_ids else 'None'}\n")
        file.write(f"Invalid Sequences: {invalid_ids if invalid_ids else 'None'}\n")


def plot_gc_content(records, output_path):
    """Save a bar chart of GC content per sequence."""
    ids = [r["id"] for r in records]
    values = [r["gc_content"] for r in records]

    plt.figure()
    plt.bar(ids, values)
    plt.xlabel("Sequence ID")
    plt.ylabel("GC Content (%)")
    plt.title("GC Content of DNA Sequences")
    plt.xticks(rotation=45)
    plt.tight_layout()
    plt.savefig(output_path)
    plt.close()


def plot_sequence_lengths(records, output_path):
    """Save a bar chart of sequence lengths."""
    ids = [r["id"] for r in records]
    values = [r["length"] for r in records]

    plt.figure()
    plt.bar(ids, values)
    plt.xlabel("Sequence ID")
    plt.ylabel("Sequence Length (bp)")
    plt.title("Length of DNA Sequences")
    plt.xticks(rotation=45)
    plt.tight_layout()
    plt.savefig(output_path)
    plt.close()


def plot_nucleotide_composition(records, output_path):
    """Save a grouped bar chart of A/T/G/C counts per sequence."""
    ids = [r["id"] for r in records]
    a_counts = [r["a"] for r in records]
    t_counts = [r["t"] for r in records]
    g_counts = [r["g"] for r in records]
    c_counts = [r["c"] for r in records]

    x = np.arange(len(ids))
    width = 0.2

    plt.figure()
    plt.bar(x - 1.5 * width, a_counts, width, label="A")
    plt.bar(x - 0.5 * width, t_counts, width, label="T")
    plt.bar(x + 0.5 * width, g_counts, width, label="G")
    plt.bar(x + 1.5 * width, c_counts, width, label="C")
    plt.xlabel("Sequence ID")
    plt.ylabel("Nucleotide Count")
    plt.title("Nucleotide Composition of DNA Sequences")
    plt.xticks(x, ids, rotation=45)
    plt.legend()
    plt.tight_layout()
    plt.savefig(output_path)
    plt.close()


def main():
    if not os.path.exists(INPUT_FILE):
        print(f"Error: input file not found at '{INPUT_FILE}'")
        return

    os.makedirs(OUTPUT_DIR, exist_ok=True)

    records = read_sequences(INPUT_FILE)

    if not records:
        print("Error: no sequences found in the FASTA file.")
        return

    invalid_ids = [r["id"] for r in records if not r["is_valid"]]
    duplicate_ids = find_duplicate_ids(records)

    write_csv(records, os.path.join(OUTPUT_DIR, "sequence_analysis.csv"))

    stats = compute_statistics(records)
    print_summary(stats, duplicate_ids, invalid_ids)
    save_summary(stats, duplicate_ids, invalid_ids, os.path.join(OUTPUT_DIR, "summary.txt"))

    plot_gc_content(records, os.path.join(OUTPUT_DIR, "gc_content.png"))
    plot_sequence_lengths(records, os.path.join(OUTPUT_DIR, "sequence_lengths.png"))
    plot_nucleotide_composition(records, os.path.join(OUTPUT_DIR, "nucleotide_composition.png"))

    print(f"\nResults saved in '{OUTPUT_DIR}/' folder.")


if __name__ == "__main__":
    main()