from Bio import SeqIO
import csv
import numpy as np
import matplotlib.pyplot as plt


# Calculate GC content
def calculate_gc_content(sequence):
    gc_count = sequence.count("G") + sequence.count("C")
    gc_content = (gc_count / len(sequence)) * 100
    return gc_content


# Calculate AT content
def calculate_at_content(sequence):
    at_count = sequence.count("A") + sequence.count("T")
    at_content = (at_count / len(sequence)) * 100
    return at_content


# Count nucleotides
def count_nucleotides(sequence):
    a_count = sequence.count("A")
    t_count = sequence.count("T")
    g_count = sequence.count("G")
    c_count = sequence.count("C")

    return a_count, t_count, g_count, c_count


# Validate DNA sequence
def validate_sequence(sequence):
    valid_bases = set("ATGC")
    return set(sequence).issubset(valid_bases)


# Main analysis
def main():

    # Store sequence lengths
    lengths = []

    # Store sequence IDs
    sequence_ids = []

    # Store GC contents
    gc_contents = []

    # Store nucleotide counts
    a_counts = []
    t_counts = []
    g_counts = []
    c_counts = []

    # Count total sequences
    sequence_count = 0

    # Store invalid sequence IDs
    invalid_sequences = []


    # Create a CSV file for results
    with open("results/sequence_analysis.csv", "w", newline="") as file:

        writer = csv.writer(file)

        # Write column names
        writer.writerow([
            "Sequence ID",
            "Length",
            "GC Content",
            "AT Content",
            "A",
            "T",
            "G",
            "C"
        ])

        # Read sequences from FASTA file
        for record in SeqIO.parse("data/sequence.fasta", "fasta"):

            sequence = record.seq

            # Count the sequence
            sequence_count += 1

            # Validate sequence
            if not validate_sequence(sequence):
                invalid_sequences.append(record.id)

            # Calculate GC content
            gc_content = calculate_gc_content(sequence)

            # Calculate AT content
            at_content = calculate_at_content(sequence)

            # Count nucleotides
            a_count, t_count, g_count, c_count = count_nucleotides(sequence)

            # Store sequence data
            lengths.append(len(sequence))
            sequence_ids.append(record.id)
            gc_contents.append(gc_content)

            # Store nucleotide counts
            a_counts.append(a_count)
            t_counts.append(t_count)
            g_counts.append(g_count)
            c_counts.append(c_count)

            # Write results to CSV
            writer.writerow([
                record.id,
                len(sequence),
                round(gc_content, 2),
                round(at_content, 2),
                a_count,
                t_count,
                g_count,
                c_count
            ])


    # Sequence length statistics
    average_length = np.mean(lengths)
    minimum_length = np.min(lengths)
    maximum_length = np.max(lengths)
    std_length = np.std(lengths)


    # GC content statistics
    average_gc = np.mean(gc_contents)
    minimum_gc = np.min(gc_contents)
    maximum_gc = np.max(gc_contents)
    std_gc = np.std(gc_contents)


    # Find shortest and longest sequences
    shortest_id = sequence_ids[lengths.index(minimum_length)]
    longest_id = sequence_ids[lengths.index(maximum_length)]


    # Check for duplicate sequence IDs
    duplicate_ids = []

    for sequence_id in sequence_ids:

        if sequence_ids.count(sequence_id) > 1:

            if sequence_id not in duplicate_ids:
                duplicate_ids.append(sequence_id)


    # Display dataset summary
    print("Total Sequences:", sequence_count)

    print("Average Length:", round(average_length, 2))
    print("Minimum Length:", minimum_length)
    print("Maximum Length:", maximum_length)
    print("Standard Deviation of Length:", round(std_length, 2))

    print("Shortest Sequence:", shortest_id, "-", minimum_length)
    print("Longest Sequence:", longest_id, "-", maximum_length)

    print("Average GC Content:", round(average_gc, 2), "%")
    print("Minimum GC Content:", round(minimum_gc, 2), "%")
    print("Maximum GC Content:", round(maximum_gc, 2), "%")
    print("Standard Deviation of GC Content:", round(std_gc, 2), "%")


    if duplicate_ids:
        print("Duplicate IDs:", duplicate_ids)
    else:
        print("Duplicate IDs: None")


    if invalid_sequences:
        print("Invalid Sequences:", invalid_sequences)
    else:
        print("Invalid Sequences: None")


    # Save summary to a text file
    with open("results/summary.txt", "w") as file:

        file.write("FASTA SEQUENCE ANALYSIS SUMMARY\n")
        file.write("===============================\n\n")

        file.write(f"Total Sequences: {sequence_count}\n")
        file.write(f"Average Length: {average_length:.2f}\n")
        file.write(f"Minimum Length: {minimum_length}\n")
        file.write(f"Maximum Length: {maximum_length}\n")
        file.write(f"Standard Deviation of Length: {std_length:.2f}\n\n")

        file.write(f"Shortest Sequence: {shortest_id} - {minimum_length}\n")
        file.write(f"Longest Sequence: {longest_id} - {maximum_length}\n\n")

        file.write(f"Average GC Content: {average_gc:.2f}%\n")
        file.write(f"Minimum GC Content: {minimum_gc:.2f}%\n")
        file.write(f"Maximum GC Content: {maximum_gc:.2f}%\n")
        file.write(f"Standard Deviation of GC Content: {std_gc:.2f}%\n\n")

        if duplicate_ids:
            file.write(f"Duplicate IDs: {duplicate_ids}\n")
        else:
            file.write("Duplicate IDs: None\n")

        if invalid_sequences:
            file.write(f"Invalid Sequences: {invalid_sequences}\n")
        else:
            file.write("Invalid Sequences: None\n")


    # Create GC content graph
    plt.figure()
    plt.bar(sequence_ids, gc_contents)
    plt.xlabel("Sequence ID")
    plt.ylabel("GC Content (%)")
    plt.title("GC Content of DNA Sequences")
    plt.xticks(rotation=45)
    plt.tight_layout()
    plt.savefig("results/gc_content.png")
    plt.show()


    # Create sequence length graph
    plt.figure()
    plt.bar(sequence_ids, lengths)
    plt.xlabel("Sequence ID")
    plt.ylabel("Sequence Length (bp)")
    plt.title("Length of DNA Sequences")
    plt.xticks(rotation=45)
    plt.tight_layout()
    plt.savefig("results/sequence_lengths.png")
    plt.show()


    # Create nucleotide composition graph
    x = np.arange(len(sequence_ids))
    width = 0.2

    plt.figure()

    plt.bar(x - 1.5 * width, a_counts, width, label="A")
    plt.bar(x - 0.5 * width, t_counts, width, label="T")
    plt.bar(x + 0.5 * width, g_counts, width, label="G")
    plt.bar(x + 1.5 * width, c_counts, width, label="C")

    plt.xlabel("Sequence ID")
    plt.ylabel("Nucleotide Count")
    plt.title("Nucleotide Composition of DNA Sequences")
    plt.xticks(x, sequence_ids, rotation=45)
    plt.legend()
    plt.tight_layout()
    plt.savefig("results/nucleotide_composition.png")
    plt.show()


# Run the program
if __name__ == "__main__":
    main()