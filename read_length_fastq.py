#Modifying fastq to only read a specific length

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

def read_specific_length(input_fastq, output_fastq, target_length, start_position):
    with open(input_fastq, "r") as input_handle, open(output_fastq, "w") as output_handle:
        for record in SeqIO.parse(input_handle, "fastq"):
            if len(record.seq) >= (start_position + target_length):
                truncated_seq = record.seq[start_position : start_position + target_length]
                truncated_qual = record.letter_annotations["phred_quality"][start_position : start_position + target_length]
                new_record = SeqRecord(truncated_seq, id=record.id, description=record.description)
                new_record.letter_annotations["phred_quality"] = truncated_qual
                SeqIO.write(new_record, output_handle, "fastq")

# Input FASTQ file path
input_fastq_file = "/Users/chintansoni/Desktop/NGS/967004-2_S63_L001_R1_001.fastq"

# Output FASTQ file path
output_fastq_file = "/Users/chintansoni/Desktop/NGS/in_967004-2_S63_L001_R1_001.fastq"

# Target length of the sequence to be read
target_length = 87  # Change this to the desired length

# Starting position to read from ....starts with 0
start_position = 4 # Change this to the desired starting position

# Call the function to read the specific length of the sequences from the starting position
read_specific_length(input_fastq_file, output_fastq_file, target_length, start_position)
