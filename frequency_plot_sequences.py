import csv
from Bio import motifs
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import matplotlib.pyplot as plt

# Initialize an empty list to store sequences
sequences = []
counts = []

# Open and read the CSV file
with open("/Users/chintansoni/Desktop/Top1% NST_mIF_noUaa_2by2_counts.csv", "r") as csvfile:
    reader = csv.DictReader(csvfile)
    for row in reader:
        count = int(row["count"])
        if count >= 1:
            sequence = row["aaRS"]
            # Convert count to an integer
            sequences.extend([sequence] * count)
            counts.extend([count])
    print (sequences, counts)


output_csv = "/Users/chintansoni/Desktop/Top1% NST_mIF_noUaa_2by2_With_Counts.csv"
with open(output_csv, "w", newline="") as outfile:
    writer = csv.writer(outfile)
    writer.writerow(["Sequence"])  # Write header
    for seq in zip(sequences):
        writer.writerow([seq])