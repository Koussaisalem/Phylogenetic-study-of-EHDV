# fetch_ncbi_data.py

from Bio import Entrez
from Bio import SeqIO
import time

# Always tell NCBI who you are
Entrez.email = "your_email@example.com"

def fetch_sequences(term, db="nucleotide", retmax=1000):
    # Search for the term
    handle = Entrez.esearch(db=db, term=term, retmax=retmax)
    record = Entrez.read(handle)
    handle.close()

    # Fetch the sequences
    id_list = record["IdList"]
    handle = Entrez.efetch(db=db, id=id_list, rettype="gb", retmode="text")
    return handle

def main():
    # Fetch VP5 sequences
    vp5_handle = fetch_sequences("hemorrhagic disease virus VP5")
    with open("raw_vp5_data.gb", "w") as out_file:
        out_file.write(vp5_handle.read())
    
    time.sleep(3)  # Be nice to NCBI servers

    # Fetch VP2 sequences
    vp2_handle = fetch_sequences("hemorrhagic disease virus VP2")
    with open("raw_vp2_data.gb", "w") as out_file:
        out_file.write(vp2_handle.read())

if __name__ == "__main__":
    main()