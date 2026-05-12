# -*- coding: utf-8 -*-
"""
Created on Mon Mar  4 11:46:33 2024

@author: selcuk.1
"""

def extract_ali_values(file_path):
    extracted_data = {}
    current_key = None
    count = 0
    
    with open(file_path, 'r') as file:
        for line in file:
            line = line.strip()  # Strip the line first
            # Check if the line is a key
            if line.startswith(">>"):
                # Extract the key, stripping newlines and leading/trailing spaces
                current_key = line.split()[1]
                count = 0
            elif current_key:
                if "[No individual domains that satisfy reporting thresholds (although complete target did)]" in line:
                    current_key = None
                    continue
                count += 1
                if count < 3:
                    continue
                
                # Split the line into components
                line_parts = line.split()
                
                # Check if line_parts has enough elements to avoid IndexError
                if len(line_parts) <= 13:
                    continue
                
                ali_from = int(line_parts[12])
                ali_to = int(line_parts[13])
                
                # Store the values in the dictionary using the current key
                if current_key not in extracted_data:
                    extracted_data[current_key] = [ali_from, ali_to]
                else:
                    if ali_from < extracted_data[current_key][0]:
                        extracted_data[current_key][0] = ali_from
                    if ali_to > extracted_data[current_key][1]:
                        extracted_data[current_key][1] = ali_to
    
    print(len(extracted_data))
    return extracted_data

from Bio import SeqIO

def extract_regions_and_write_fasta(fasta_file_path, regions_dict):
    """
    Extracts regions from sequences in a FASTA file based on a dictionary of regions and writes them to a new FASTA file.
    
    Args:
    - fasta_file_path: Path to the input FASTA file.
    - output_file_path: Path where the output FASTA file will be written.
    - regions_dict: Dictionary with keys matching sequence identifiers and values being lists with [alifrom, ali_to] regions.
    """
    count=0
    count2=0
    output_file_path=fasta_file_path.replace(".fasta","_hmmregions.fasta")
    with open(output_file_path, 'w') as output_handle:
        for record in SeqIO.parse(fasta_file_path, "fasta"):
            # Construct the key as it appears in the regions_dict
            key = record.id
            if key in regions_dict:
                count2+=1
                # Adjust alifrom by decreasing it by 1 for inclusive boundary
                alifrom = regions_dict[key][0] - 1
                # Keep ali_to as it is
                ali_to = regions_dict[key][1]
                # Extract the region
                extracted_region = record.seq[alifrom:ali_to]
                if len(extracted_region)>19:
                    # Write the extracted region to the new FASTA file
                    output_handle.write(f">{record.description}\n{extracted_region}\n")
                    count+=1
    print(count2,"match header")
    print(count,"sequences are extracted!")
                  


import argparse

def main(hmm_file_path, fasta_file_path):
    filter_dict=extract_ali_values(hmm_file_path)
    extract_regions_and_write_fasta(fasta_file_path, filter_dict)


if __name__ == "__main__":
    # Initialize the parser
    parser = argparse.ArgumentParser(description="Process a text file and a FASTA file to extract specific regions from sequences based on HMM alignment.")

    # Adding arguments
    parser.add_argument("hmm_file_path", help="Path to the input hmm file.")
    parser.add_argument("fasta_file_path", help="Path to the input FASTA file.")

    # Parse the arguments
    args = parser.parse_args()

    # Execute the main functionality with the provided arguments
    main(args.hmm_file_path, args.fasta_file_path)