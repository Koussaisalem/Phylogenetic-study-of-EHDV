import pandas as pd

def dataframe_to_fasta(csv_file, fasta_file):
    # Read the CSV file into a DataFrame
    df = pd.read_csv(csv_file)
    
    with open(fasta_file, 'w') as fasta_output:
        for _, row in df.iterrows():
            accession = row['accession']
            collection_date = row['collection_date']
            sequence = row['sequence']
            source = row['source']
            country = row['country']
            
            # Write in FASTA format
            header = f'>{accession}|{collection_date}|{source} |{country} '
            fasta_output.write(f'{header}\n')
            
            # Write sequence in lines of 60 characters
            for i in range(0, len(sequence), 60):
                fasta_output.write(f'{sequence[i:i+60]}\n')

# Example usage
csv_file = 'vp5_merged_dataset.csv'
fasta_file = 'VP5.fasta'
dataframe_to_fasta(csv_file, fasta_file)

