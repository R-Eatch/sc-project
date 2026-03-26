#!/usr/bin/env python3

from Bio import SeqIO
import re # Import regular expressions module for parsing
import os
from collections import defaultdict

# --- Configuration ---
# Your list of target gene symbols for MOUSE (case-sensitive based on FASTA headers)
gene_list_mouse = ['Cavin3', 'Dsg3', 'Lalba', 'Pip', 'Plet1','Top2a'] # Mouse specific list
gene_list_mouse = ['Afm',
 'Art3',
 'Cavin3',
 'Cd300lg',
 'Cyp2s1',
 'Dgat2l6',
 'Dsg3',
 'Icam2',
 'Klk11',
 'Klk6',
 'Klk9',
 'Krt1',
 'Krt32',
 'Krt77',
 'Krt78',
 'Lalba',
 'Oas2',
 'Pip',
 'Plet1',
 'Plscr4',
 'S100a3',
 'Serpinb7',
 'Trim10',
 'Trim26',
 'Zfp263',
 'Zfp444',
 'Zkscan5',
 'Zkscan6',
 'Zkscan8',
 'Zscan20',
'Csn1s1','Top2a','Csn2','Ltf' ## three test gene
]
# Your list of target gene symbols for HUMAN (case-sensitive based on FASTA headers)
gene_list_human = ['CAVIN3', 'DSG3', 'LALBA', 'PIP', 'PLET1','TOP2A'] # Human specific list (Note potential case differences)

# Combine all target genes into a single set for searching across both files
# Use union operator | to combine sets, automatically handles duplicates
all_target_genes = set(gene_list_mouse) | set(gene_list_human)
if not all_target_genes:
    print("Warning: Both mouse and human gene lists are empty. No genes to search for.")
    # Exit early if nothing to do
    # exit(0) # Or let it run and report nothing found

# *** MODIFY THESE PATHS ***
# Path to your input MOUSE peptide FASTA file (Ensembl format)
input_fasta_file_mouse = "../data/Mus_musculus.GRCm39.pep.all.fa"
# Path to your input HUMAN peptide FASTA file (Ensembl format)
input_fasta_file_human = "../data/Homo_sapiens.GRCh38.pep.all.fa"

# Name for the output FASTA file containing the LONGEST peptide for each target gene found
output_fasta_file = "longest_peptides_human_mouse.fa" # Changed filename to reflect content
# ---------------------

# --- Input Validation ---
input_files_to_check = {
    "Mouse": input_fasta_file_mouse,
    "Human": input_fasta_file_human
}
missing_files = False
valid_files_processed = [] # Keep track of files that exist and will be processed
for species, filepath in input_files_to_check.items():
    if not filepath: # Check if path is empty or None
         print(f"Warning: Input FASTA file path for {species} is not set. Skipping this file.")
         # Mark as None so processing step skips it safely
         input_files_to_check[species] = None
    elif not os.path.exists(filepath):
        print(f"Error: Input FASTA file for {species} not found at '{filepath}'")
        print(f"Please modify the 'input_fasta_file_{species.lower()}' variable in the script.")
        missing_files = True
    else:
        valid_files_processed.append(species) # File exists

if missing_files:
    exit(1)

if not valid_files_processed:
    print("Error: No valid input FASTA files were specified or found. Exiting.")
    exit(1)

if not all_target_genes:
    print("No target genes specified in either list. Exiting.")
    exit(0)


# --- Data Structures ---
# Dictionary to store ALL sequences grouped by gene symbol from BOTH files
# Structure: {gene_symbol: [(header, sequence), (header, sequence), ...]}
gene_sequences = defaultdict(list)
# Set to keep track of which target gene symbols were found at least once across all files
genes_found_in_files = set()

# --- Helper Function for Processing ---
def process_fasta_collect_all(fasta_path, species_name, target_genes_set, gene_sequences_dict, found_set):
    """
    Reads a FASTA file, finds sequences for ANY gene in the target_genes_set,
    and stores them grouped by gene symbol in gene_sequences_dict.

    Args:
        fasta_path (str): Path to the input FASTA file.
        species_name (str): Name of the species (e.g., "Mouse", "Human") for logging.
        target_genes_set (set): The COMBINED set of gene symbols to search for.
        gene_sequences_dict (defaultdict): Dictionary to append found (header, sequence) tuples to, keyed by gene symbol.
        found_set (set): Set to add the symbol of any found gene to.
    """
    print(f"\n--- Processing {species_name} sequences from: {fasta_path} ---")
    print(f"Searching for any of these genes: {', '.join(sorted(list(target_genes_set)))}") # Log the full set being searched

    sequences_found_in_this_file = 0
    try:
        # Iterate through each sequence record in the FASTA file
        for record in SeqIO.parse(fasta_path, "fasta"):
            header = record.description # Get the full header line content
            sequence = str(record.seq) # Get the sequence as a string

            # Use regex to find the gene symbol in the header
            # Assumes Ensembl format 'gene_symbol:SYMBOL'
            match = re.search(r'gene_symbol:(\S+)', header)

            if match:
                gene_symbol = match.group(1) # Extract the matched gene symbol

                # Check if this gene is in our combined target list
                if gene_symbol in target_genes_set:
                    # Store the header and sequence for this gene symbol
                    gene_sequences_dict[gene_symbol].append((header, sequence))
                    found_set.add(gene_symbol) # Mark this gene symbol as found (at least once)
                    sequences_found_in_this_file += 1

    except FileNotFoundError:
        # Redundant check, but good practice
        print(f"Error: Input FASTA file for {species_name} not found during processing at '{fasta_path}'")
        return -1 # Indicate an error occurred
    except Exception as e:
        print(f"An error occurred while reading the {species_name} FASTA file: {e}")
        return -1 # Indicate an error occurred

    print(f"Finished reading {species_name} file. Found {sequences_found_in_this_file} sequences matching target gene symbols in this file.")
    return sequences_found_in_this_file

# --- Main Processing ---

# Process Mouse File (if path is valid)
if input_files_to_check.get("Mouse"):
    process_fasta_collect_all(
        fasta_path=input_fasta_file_mouse,
        species_name="Mouse",
        target_genes_set=all_target_genes, # Use the combined set
        gene_sequences_dict=gene_sequences,
        found_set=genes_found_in_files
    )

# Process Human File (if path is valid)
if input_files_to_check.get("Human"):
     process_fasta_collect_all(
        fasta_path=input_fasta_file_human,
        species_name="Human",
        target_genes_set=all_target_genes, # Use the combined set
        gene_sequences_dict=gene_sequences,
        found_set=genes_found_in_files
    )

print("\n--- Identifying Longest Peptides ---")

if not gene_sequences:
    print("No sequences found for any of the target genes in the provided files.")
    # Report missing genes based on the initial combined list
    print(f"Target genes searched for: {', '.join(sorted(list(all_target_genes)))}")
    exit(0)

# --- Find Longest Peptide for Each Gene Found ---
longest_peptides_to_write = [] # List to store the final (header, sequence) tuples

# Iterate through the genes we actually found sequences for
for gene_symbol in sorted(list(genes_found_in_files)): # Iterate based on what was found
    sequences_for_gene = gene_sequences[gene_symbol]

    if sequences_for_gene:
        # Find the tuple (header, sequence) with the maximum sequence length
        longest_entry = max(sequences_for_gene, key=lambda item: len(item[1]))
        longest_peptides_to_write.append(longest_entry)
        print(f"  - Gene '{gene_symbol}': Found {len(sequences_for_gene)} total sequence(s) across files. Selected longest (Length: {len(longest_entry[1])}).")
    else:
        # This case shouldn't happen if gene_symbol is from genes_found_in_files
        print(f"Warning: Gene {gene_symbol} was marked as found but no sequences were stored.")


# --- Report on Missing Genes ---
missing_genes = all_target_genes - genes_found_in_files
if missing_genes:
    print(f"\nWarning: The following target genes were NOT found in any processed input file: {', '.join(sorted(list(missing_genes)))}")
else:
    print("\nAll target genes were found in the processed input files.")


# --- Write Output FASTA File ---
if not longest_peptides_to_write:
     print("\nNo longest peptides identified to write (this shouldn't happen if genes were found).")
     exit(0)

print(f"\nWriting {len(longest_peptides_to_write)} longest peptide sequences to: {output_fasta_file}")

try:
    with open(output_fasta_file, "w") as outfile:
        for header, sequence in longest_peptides_to_write:
            outfile.write(f">{header}\n") # Write the header line (starts with '>')
            outfile.write(f"{sequence}\n") # Write the sequence

except IOError as e:
    print(f"Error writing output file '{output_fasta_file}': {e}")
    exit(1)

print("\nScript finished successfully.")
print(f"Output file '{output_fasta_file}' created, containing the longest peptide found for each targeted gene symbol across all input files.")
