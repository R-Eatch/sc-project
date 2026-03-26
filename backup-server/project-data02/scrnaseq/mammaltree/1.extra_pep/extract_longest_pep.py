#!/usr/bin/env python3

from Bio import SeqIO
import re # Import regular expressions module for parsing
import os
from collections import defaultdict

# --- Configuration ---
# Your list of target gene symbols
gene_list = ['Cavin3', 'Dsg3', 'Lalba', 'Pip', 'Plet1']

# *** MODIFY THIS PATH ***
# Path to your input mouse peptide FASTA file (Ensembl format)
input_fasta_file = "../data/Mus_musculus.GRCm39.pep.all.fa" 

# Name for the output FASTA file containing the longest peptides
output_fasta_file = "longest_queries.fa"
# ---------------------

# --- Input Validation ---
if not os.path.exists(input_fasta_file):
    print(f"Error: Input FASTA file not found at '{input_fasta_file}'")
    print("Please modify the 'input_fasta_file' variable in the script.")
    exit(1)

# --- Data Structures ---
# Dictionary to store sequences grouped by gene symbol
# Structure: {gene_symbol: [(header, sequence), (header, sequence), ...]}
gene_sequences = defaultdict(list) 

# --- Main Processing ---
print(f"Reading sequences from: {input_fasta_file}")
print(f"Looking for genes: {', '.join(gene_list)}")

try:
    # Iterate through each sequence record in the FASTA file
    for record in SeqIO.parse(input_fasta_file, "fasta"):
        header = record.description # Get the full header line content
        sequence = str(record.seq) # Get the sequence as a string

        # Use regex to find the gene symbol in the header
        # This looks for 'gene_symbol:' followed by non-space characters
        match = re.search(r'gene_symbol:(\S+)', header) 
        
        if match:
            gene_symbol = match.group(1) # Extract the matched gene symbol
            
            # Check if this gene is in our target list
            if gene_symbol in gene_list:
                # Store the header and sequence for this gene
                gene_sequences[gene_symbol].append((header, sequence))

except FileNotFoundError:
    # This check is technically redundant due to the os.path.exists check earlier,
    # but good practice in case the file disappears between check and open.
    print(f"Error: Input FASTA file not found at '{input_fasta_file}'")
    exit(1)
except Exception as e:
    print(f"An error occurred while reading the FASTA file: {e}")
    exit(1)

print("Finished reading file. Identifying longest peptides...")

# --- Find Longest Peptide for Each Gene ---
longest_peptides_to_write = [] # List to store the final (header, sequence) tuples
genes_found = set()

for gene_symbol in gene_list:
    if gene_symbol in gene_sequences:
        sequences_for_gene = gene_sequences[gene_symbol]
        
        if sequences_for_gene:
            # Find the tuple (header, sequence) with the maximum sequence length
            longest_entry = max(sequences_for_gene, key=lambda item: len(item[1]))
            longest_peptides_to_write.append(longest_entry)
            genes_found.add(gene_symbol)
            print(f"  - Found {len(sequences_for_gene)} sequence(s) for {gene_symbol}. Longest length: {len(longest_entry[1])}")
        else:
            # This case should technically not happen with defaultdict if the symbol was found
            print(f"Warning: Gene {gene_symbol} was processed but no sequences were stored.")
    else:
        print(f"Warning: Gene '{gene_symbol}' from the list was not found in the input file's gene symbols.")

# --- Write Output FASTA File ---
print(f"\nWriting {len(longest_peptides_to_write)} longest peptide sequences to: {output_fasta_file}")

try:
    with open(output_fasta_file, "w") as outfile:
        for header, sequence in longest_peptides_to_write:
            outfile.write(f">{header}\n") # Write the header line (starts with '>')
            # Write the sequence. You can optionally wrap lines, but for BLAST queries, often one line is fine.
            outfile.write(f"{sequence}\n") 

except IOError as e:
    print(f"Error writing output file '{output_fasta_file}': {e}")
    exit(1)

print("\nScript finished successfully.")
print(f"Output file '{output_fasta_file}' created.")
missing_genes = set(gene_list) - genes_found
if missing_genes:
    print(f"Note: The following genes from your list were not found in the input file: {', '.join(missing_genes)}")
