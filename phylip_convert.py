from Bio import SeqIO

records = SeqIO.parse("vp2_curated.fasta", "fasta")
count = SeqIO.write(records, "vp2_curated.phylip", "phylip")
print("Converted %i records" % count)