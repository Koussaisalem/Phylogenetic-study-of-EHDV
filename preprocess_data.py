# preprocess_data.py

# preprocess_data.py

from __future__ import print_function  # For Python 2 compatibility
from Bio import SeqIO
from Bio.Blast import NCBIWWW, NCBIXML
import xml.etree.ElementTree as ET
import pandas as pd
from Bio import Entrez
import re
import time
from tqdm import tqdm
import os

def is_mediterranean(country):
    mediterranean_countries = [
        "Spain", "France", "Monaco", "Italy", "Slovenia", "Croatia", "Bosnia and Herzegovina",
        "Montenegro", "Albania", "Greece", "Turkey", "Syria", "Lebanon", "Israel", 
        "Palestine", "Egypt", "Libya", "Tunisia", "Algeria", "Morocco", "Cyprus", "Malta"
    ]
    return any(country.lower() in c.lower() for c in mediterranean_countries)

def extract_info(record):
    country = ""
    collection_date = ""
    
    for feature in record.features:
        if feature.type == "source":
            if "country" in feature.qualifiers:
                country = feature.qualifiers["country"][0]
            elif "geo_loc_name" in feature.qualifiers:
                country = feature.qualifiers["geo_loc_name"][0]
            if "collection_date" in feature.qualifiers:
                collection_date = feature.qualifiers["collection_date"][0]
    
    return {
        "accession": record.id,
        "country": country,
        "collection_date": collection_date,
        "sequence": str(record.seq)
    }

def process_genbank_file(file_path, segment):
    data = []
    sequences = []

    if not os.path.exists(file_path):
        print("Error: File '{}' not found.".format(file_path))
        return None

    for record in SeqIO.parse(file_path, "genbank"):
        info = extract_info(record)
        if info["country"] and is_mediterranean(info["country"]):
            data.append(info)
            sequences.append(">{0}\n{1}".format(info['accession'], info['sequence']))

    if not data:
        print("No Mediterranean sequences found in the {} file.".format(segment))
        return None

    # Create DataFrame
    df = pd.DataFrame(data)
    df.to_csv("{0}_dataset.csv".format(segment), index=False)

    # Write FASTA file
    with open("{0}_sequences.fasta".format(segment), "w") as f:
        f.write("\n".join(sequences))

    print("Processed {0} records for {1}".format(len(data), segment))
    return df

Entrez.email = "p9lite.rk@gmail.com"  # Replace with your email address


def blast_search(accession, segment, max_wait_time=900):  # 15 minutes max wait time
    print("Performing BLAST search for {} (accession: {})...".format(segment, accession))
    try:
        # Step 1: Run the BLAST search
        handle = Entrez.epost("nuccore", id=accession)
        search_results = Entrez.read(handle)
        webenv = search_results["WebEnv"]
        query_key = search_results["QueryKey"]

        blast_handle = Entrez.elink(dbfrom="nuccore", db="nuccore", query_key=query_key, webenv=webenv, cmd="neighbor_score")
        blast_results = Entrez.read(blast_handle)

        # Check if blast_results is empty
        if not blast_results or not blast_results[0].get('LinkSetDb'):
            print("No BLAST results found for accession: {}".format(accession))
            return []

        # Step 2: Process the BLAST results
        blast_data = []
        for result in blast_results[0]['LinkSetDb'][0]['Link']:
            hit_id = result['Id']
            
            # Fetch the sequence information
            handle = Entrez.efetch(db="nucleotide", id=hit_id, rettype="gb", retmode="xml")
            records = Entrez.read(handle)
            
            if not records:
                print("No record found for hit_id: {}".format(hit_id))
                continue

            record = records[0]
            
            accession = record['GBSeq_primary-accession']
            description = record['GBSeq_definition']
            sequence = record['GBSeq_sequence']
            
            # Extract country information
            country = ""
            for feature in record['GBSeq_feature-table']:
                if feature['GBFeature_key'] == 'source':
                    for qualifier in feature['GBFeature_quals']:
                        if qualifier['GBQualifier_name'] == 'country':
                            country = qualifier['GBQualifier_value']
                            break
                    if country:
                        break
            
            blast_data.append({
                "accession": accession,
                "description": description,
                "country": country,
                "sequence": sequence
            })

        return blast_data

    except Exception as e:
        print("Error during BLAST search: {}".format(str(e)))
        return []

def enrich_dataset(df, segment):
    if df is None or df.empty:
        print("No data to enrich for {} dataset.".format(segment))
        return None

    enriched_data = set()  # Use a set to avoid redundancy
    original_accessions = set(df['accession'])
    
    for idx, row in tqdm(df.iterrows(), total=df.shape[0], desc="Enriching {} dataset".format(segment)):
        blast_results = blast_search(row['accession'], segment)
        if not blast_results:
            print("No BLAST results found for accession: {}".format(row['accession']))
        else:
            print("Found {} BLAST results for accession: {}".format(len(blast_results), row['accession']))
            
        for result in blast_results:
            print("Checking BLAST result: {}".format(result['accession']))
            if result['accession'] in original_accessions:
                print("  Excluded: Already in original dataset")
            #elif not result['country']:
             #   print("  Excluded: No country information")
            else:
                print("  Added to enriched dataset")
                enriched_data.add((
                    result['accession'],
                    result['country'],
                    "",  # BLAST results don't provide collection date
                    result['sequence']
                ))
        
        # Save partial results after each BLAST search
        if enriched_data:
            partial_df = pd.DataFrame(list(enriched_data), columns=['accession', 'country', 'collection_date', 'sequence'])
            partial_df.to_csv("{}_enriched_dataset_partial.csv".format(segment), index=False)
            
            with open("{}_enriched_sequences_partial.fasta".format(segment), "w") as f:
                for _, partial_row in partial_df.iterrows():
                    f.write(">{}\n{}\n".format(partial_row['accession'], partial_row['sequence']))
            
            print("Partial results saved. Current enriched dataset for {} has {} sequences".format(segment, len(partial_df)))
        else:
            print("No new sequences added to the enriched dataset.")
    
    if not enriched_data:
        print("No enriched data found for {} dataset.".format(segment))
        return None

    enriched_df = pd.DataFrame(list(enriched_data), columns=['accession', 'country', 'collection_date', 'sequence'])
    enriched_df.to_csv("{}_enriched_dataset.csv".format(segment), index=False)
    
    with open("{}_enriched_sequences.fasta".format(segment), "w") as f:
        for _, row in enriched_df.iterrows():
            f.write(">{}\n{}\n".format(row['accession'], row['sequence']))
    
    print("Final enriched dataset for {} has {} sequences".format(segment, len(enriched_df)))
    return enriched_df

def main():
    for segment in ["vp5", "vp2"]:
        file_path = "raw_{}_data.gb".format(segment)
        df = process_genbank_file(file_path, segment)
        if df is not None and not df.empty:
            enriched_df = enrich_dataset(df, segment)
        else:
            print("Skipping enrichment for {} due to lack of initial data.".format(segment))

if __name__ == "__main__":
    main()