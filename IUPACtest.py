from Bio.Seq import Seq
from Bio.Data import IUPACData


sequence_data = 'AUGCGGG/AUCACAU'
sequence = sequence_data.replace('/ ','')
print (sequence)
sequence_holder = []
#sequence = IUPACData.unambiguous_rna_letters
sequence_holder.append((sequence))

print (type(sequence_holder))

print (str(sequence_holder))

