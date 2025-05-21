import pandas as pd

# Define the input and output file paths
input_file = '/Users/chintansoni/Desktop/NGS/inpout_barcodes.csv'
output_file = '/Users/chintansoni/Desktop/NGS/unique_sequences.csv'

# Read the sequences from the input file into a DataFrame
df = pd.read_csv(input_file, header=None, names=['Sequence'])

# Remove duplicate sequences
df_unique = df.drop_duplicates()

# Write the unique sequences to the output CSV file
df_unique.to_csv(output_file, index=False, header=False)

print(f"Unique sequences have been written to {output_file}")
