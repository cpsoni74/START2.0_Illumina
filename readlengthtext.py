import csv

# Define the input and output file paths
input_file = '/Users/chintansoni/Desktop/NGS/in mismatch filtered.txt'
output_file = '/Users/chintansoni/Desktop/NGS/inpout_barcodes.csv'

# Define the start and end positions for the 15-nucleotide sequence
start_pos = 26  # Base 27 (0-based index)
end_pos = 41    # Base 41 (0-based index, exclusive)

# Read the sequences from the input file
with open(input_file, 'r') as infile:
    sequences = infile.readlines()

# Process and trim each sequence
trimmed_sequences = []
i = 0
for sequence in sequences:
    print (i)
    i = i + 1
    sequence = sequence.strip()  # Remove any surrounding whitespace
    if len(sequence) >= end_pos:  # Ensure the sequence is long enough
        trimmed_sequence = sequence[start_pos:end_pos]
        trimmed_sequences.append([trimmed_sequence])  # Add as list to write as a single column

# Write the trimmed sequences to the output CSV file
with open(output_file, 'w', newline='') as outfile:
    csv_writer = csv.writer(outfile)
    # Write header (optional)
    #csv_writer.writerow(['Trimmed_Sequence'])
    # Write rows
    csv_writer.writerows(trimmed_sequences)

print(f"Trimmed sequences have been written to {output_file}")
