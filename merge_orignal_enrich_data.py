from Bio import Entrez
import pandas as pd
from tqdm import tqdm
import time

def save_as_fasta(df, filename):
    with open(filename, 'w') as f:
        for _, row in df.iterrows():
            f.write(f">{row['accession']}|{row['country']}|{row['collection_date']}|{row['source']}\n")
            f.write(f"{row['sequence']}\n")
    print(f"Saved {len(df)} sequences to {filename}")

Entrez.email = "p9lite.rk@gmail.com"  # Replace with your email address

def fetch_sequence_info(accession):
    try:
        handle = Entrez.efetch(db="nucleotide", id=accession, rettype="gb", retmode="xml")
        record = Entrez.read(handle)[0]
        
        # Extract country and collection date
        country = ""
        collection_date = ""
        for feature in record['GBSeq_feature-table']:
            if feature['GBFeature_key'] == 'source':
                for qualifier in feature['GBFeature_quals']:
                    if qualifier['GBQualifier_name'] == 'geo_loc_name':
                        country = qualifier['GBQualifier_value']
                    elif qualifier['GBQualifier_name'] == 'collection_date':
                        collection_date = qualifier['GBQualifier_value']
        
        return country, collection_date
    except Exception as e:
        print(f"Error fetching info for {accession}: {str(e)}")
        return "", ""

def update_dataset(input_file, output_file):
    df = pd.read_csv(input_file)
    
    # Add new columns if they don't exist
    if 'country' not in df.columns:
        df['country'] = ""
    if 'collection_date' not in df.columns:
        df['collection_date'] = ""
    
    for index, row in tqdm(df.iterrows(), total=df.shape[0], desc=f"Updating {input_file}"):
        if not row['country'] or not row['collection_date']:
            country, collection_date = fetch_sequence_info(row['accession'])
            df.at[index, 'country'] = country
            df.at[index, 'collection_date'] = collection_date
        
        # Add a small delay to avoid overwhelming the NCBI server
        time.sleep(0.5)
    
    df.to_csv(output_file, index=False)
    return df

def merge_datasets(original_file, enriched_file, output_csv, output_fasta ):
    original_df = pd.read_csv(original_file)
    enriched_df = pd.read_csv(enriched_file)
    
    # Add a column to distinguish original vs BLAST sequences
    original_df['source'] = 'original'
    enriched_df['source'] = 'BLAST'
    
    # Merge datasets
    merged_df = pd.concat([original_df, enriched_df], ignore_index=True)
    
    # Remove duplicates based on accession number, keeping the original if duplicate
    merged_df = merged_df.sort_values('source').drop_duplicates('accession', keep='first')
    
    merged_df.to_csv(output_csv, index=False)

     # Save as FASTA
    save_as_fasta(merged_df, output_fasta)
    
    return merged_df

def main():
    for segment in ['vp2', 'vp5']:
        print(f"Processing {segment} dataset...")
        
        # Update original dataset
        original_file = f"{segment}_dataset.csv"
        updated_original_file = f"{segment}_dataset_updated.csv"
        update_dataset(original_file, updated_original_file)
        
        # Update enriched dataset
        enriched_file = f"{segment}_enriched_dataset.csv"
        updated_enriched_file = f"{segment}_enriched_dataset_updated.csv"
        update_dataset(enriched_file, updated_enriched_file)
        
        # Merge datasets
        merged_file = f"{segment}_merged_dataset.csv"
        merged_fasta_file = f"{segment}_merged_dataset.fasta"
        merge_datasets(updated_original_file, updated_enriched_file, merged_file, merged_fasta_file)
        
        print(f"Completed processing {segment} dataset.")
        print(f"Merged CSV file saved as {merged_file}")
        print(f"Merged FASTA file saved as {merged_fasta_file}")
        
if __name__ == "__main__":
    main()