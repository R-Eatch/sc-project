#!/usr/bin/env python3

import os
import subprocess # To run external commands like cat, makeblastdb, blastp
import re
from Bio import SeqIO
import sys # To exit on error
import glob # To find files matching a pattern

# --- Configuration ---
# Directory where this script is located and where results will be stored
output_dir = "./"

# Directory containing the input peptide FASTA files (downloaded species + mouse)
# *** MODIFY THIS PATH ***
data_dir = "../data/" 

# Name of the mouse peptide file (to exclude from the database)
# *** MODIFY THIS PATH if different ***
mouse_pep_file = "Mus_musculus.GRCm39.pep.all.fa"

# Path to the query FASTA file (longest peptides extracted in step 1)
# *** MODIFY THIS PATH ***
query_file = "../1.extra_pep/longest_queries.fa"

# Base name for the combined database file and BLAST database
db_base_name = "combined_species"
combined_fasta_db_file = f"{db_base_name}.pep.fa"
blast_db_name = f"{db_base_name}_db" # Name used by blastp's -db option

# BLASTp Parameters
# *** ADJUST AS NEEDED ***
num_threads = 20  # Number of CPU cores to use for BLASTp
evalue_threshold = "1e-4" # E-value cutoff
max_target_seqs = "100" # Max number of hits to retrieve (we can filter top 10 later)
# Output format 6: tabular, easy to parse. Customize columns as needed.
# qseqid: Query Seq-id, sseqid: Subject Seq-id, pident: Percentage of identical matches
# length: Alignment length, mismatch: Number of mismatches, gapopen: Number of gap openings
# qstart: Start of alignment in query, qend: End of alignment in query
# sstart: Start of alignment in subject, send: End of alignment in subject
# evalue: Expect value, bitscore: Bit score, stitle: Subject Title (useful!)
outfmt_string = "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore stitle"
# ---------------------


# --- Function to run shell commands ---
def run_command(command_list, cwd=None, check=True):
    """Runs a shell command using subprocess and checks for errors."""
    print(f"Running command: {' '.join(command_list)}")
    try:
        process = subprocess.run(command_list, cwd=cwd, check=check, 
                                 capture_output=True, text=True)
        print(process.stdout)
        if process.stderr:
            print("Stderr:", process.stderr, file=sys.stderr) # Show stderr messages
        return process
    except FileNotFoundError as e:
        print(f"Error: Command not found - {e}. Is BLAST+ installed and in PATH?")
        sys.exit(1)
    except subprocess.CalledProcessError as e:
        print(f"Error executing command: {' '.join(command_list)}")
        print(f"Return code: {e.returncode}")
        print(f"Output:\n{e.output}")
        print(f"Stderr:\n{e.stderr}")
        sys.exit(1)
    except Exception as e:
        print(f"An unexpected error occurred: {e}")
        sys.exit(1)

# --- Main Workflow ---
# 1. Create output directory
print(f"Creating output directory: {output_dir}")
os.makedirs(output_dir, exist_ok=True)

# 2. Combine non-mouse FASTA files
print("\nStep 2.1: Combining non-mouse FASTA files...")
# Find all .fa files in data_dir, excluding the mouse file
all_files = glob.glob(os.path.join(data_dir, "*.fa")) + \
            glob.glob(os.path.join(data_dir, "*.fasta")) # Add *.fasta just in case
files_to_combine = [f for f in all_files if os.path.basename(f) != mouse_pep_file]

if not files_to_combine:
    print(f"Error: No FASTA files found in '{data_dir}' to combine (excluding '{mouse_pep_file}').")
    sys.exit(1)

print(f"Found {len(files_to_combine)} files to combine:")
for f in files_to_combine:
    print(f"  - {os.path.basename(f)}")

# Use 'cat' command to combine files
# Note: Using shell=True here for simplicity with wildcard/piping, but be cautious with user input
# A safer way is to pass the list directly to cat if not too long.
# We will build the command list carefully.
cat_command = ["cat"] + files_to_combine
combined_fasta_path = os.path.join(output_dir, combined_fasta_db_file)

try:
    with open(combined_fasta_path, "w") as outfile:
        print(f"Running command: cat {' '.join(files_to_combine)} > {combined_fasta_path}")
        process = subprocess.run(cat_command, check=True, stdout=outfile, text=True, stderr=subprocess.PIPE)
        if process.stderr:
             print("Stderr from cat:", process.stderr, file=sys.stderr)
    print(f"Successfully combined files into: {combined_fasta_path}")
except Exception as e:
     print(f"Error combining files: {e}")
     sys.exit(1)


# 3. Build BLAST database using makeblastdb
print("\nStep 2.2: Building BLAST database...")
makeblastdb_command = [
    "makeblastdb",
    "-in", combined_fasta_path,
    "-dbtype", "prot",
    "-out", os.path.join(output_dir, blast_db_name), # Specify output path for db files
    "-title", "Combined Species Peptide Database",
    "-parse_seqids" # Recommended for easier sequence retrieval later
]
run_command(makeblastdb_command)
print(f"BLAST database '{blast_db_name}' created in '{output_dir}'.")

# 4. Perform BLASTp search for each query sequence
print("\nStep 2.3: Running BLASTp for each query sequence...")

# Check if query file exists
if not os.path.exists(query_file):
    print(f"Error: Query file not found at '{query_file}'")
    sys.exit(1)

try:
    # Iterate through each sequence in the multi-FASTA query file
    for record in SeqIO.parse(query_file, "fasta"):
        query_id = record.id 
        query_description = record.description
        # Try to extract gene symbol for a cleaner filename
        gene_symbol_match = re.search(r'gene_symbol:(\S+)', query_description)
        if gene_symbol_match:
            file_prefix = gene_symbol_match.group(1)
        else:
            # Fallback to using the first part of the sequence ID
            file_prefix = query_id.split('.')[0] 
            print(f"Warning: Could not parse 'gene_symbol:' from header for {query_id}. Using '{file_prefix}' as file prefix.")
        
        output_blast_file = os.path.join(output_dir, f"{file_prefix}.blastp.tsv")
        
        # Create a temporary file for the single query sequence
        # BLASTp can take a multi-fasta, but running individually makes output management easier
        temp_query_file = os.path.join(output_dir, f"temp_query_{file_prefix}.fa")
        SeqIO.write(record, temp_query_file, "fasta")
        
        print(f"  Running BLASTp for query: {file_prefix} ({query_id})")
        
        blastp_command = [
            "blastp",
            "-query", temp_query_file,
            "-db", os.path.join(output_dir, blast_db_name), # Path to the database name
            "-out", output_blast_file,
            "-num_threads", str(num_threads),
            "-evalue", evalue_threshold,
            "-max_target_seqs", max_target_seqs,
            "-outfmt", outfmt_string
        ]
        
        run_command(blastp_command, check=False) # Don't exit immediately if BLAST finds no hits
        
        # Clean up temporary query file
        os.remove(temp_query_file)
        
        print(f"  BLASTp results saved to: {output_blast_file}")

except FileNotFoundError:
    print(f"Error: Query FASTA file not found at '{query_file}'")
    sys.exit(1)
except Exception as e:
    print(f"An error occurred during BLASTp processing: {e}")
    sys.exit(1)

print("\nBLASTp searches completed.")
print(f"All results are saved in the '{output_dir}' directory.")
