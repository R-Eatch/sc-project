#!/usr/bin/env python3

import os
import subprocess # To run external commands like cat, makeblastdb, blastp
import sys # To exit on error
import glob # To find files matching a pattern
import re # For parsing gene symbols and potentially species IDs
from collections import defaultdict
from Bio import SeqIO

# --- Configuration ---

# Directory where this script is located and where results will be stored
main_output_dir = "./blast_results_reciprocal" # Main directory for all results

# Directory containing ALL input peptide FASTA files (downloaded species + mouse + human)
# *** MODIFY THIS PATH ***
data_dir = "../data/"

# Names of the peptide files for Mouse and Human (used for exclusion)
# *** MODIFY THESE PATHS if different ***
mouse_pep_file = "Mus_musculus.GRCm39.pep.all.fa"
human_pep_file = "Homo_sapiens.GRCh38.pep.all.fa"

# Path to the SOURCE query FASTA file (output from step 1, containing ALL extracted HUMAN+MOUSE peptides)
# *** MODIFY THIS PATH ***
source_query_file = "../1.extra_pep/longest_peptides_human_mouse.fa"

# BLASTp Parameters
# *** ADJUST AS NEEDED ***
num_threads = 20
evalue_threshold = "1e-5"
max_target_seqs = "100" # Max hits per query sequence within the gene's file
# Output format 6: tabular. Columns defined below.
outfmt_string = "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore stitle"
# Define the 0-based index of relevant columns based on outfmt_string
qseqid_column_index = 0   # Query Seq-id
sseqid_column_index = 1   # Subject Seq-id
bitsore_column_index = 11 # Bit score
# stitle_column_index = 12 # Subject Title (defined in outfmt, but not used for filtering)
# ---------------------

# --- Function Definitions (run_command, sort_blast_results, filter_blast_results_self_hits) ---
# These functions remain the same as in the original script.
# Include the full code for these functions here...
# --- Function to run shell commands ---
def run_command(command_list, cwd=None, check=True):
    """Runs a shell command using subprocess and checks for errors."""
    print(f"Running command: {' '.join(command_list)}")
    try:
        process = subprocess.run(command_list, cwd=cwd, check=check,
                                 capture_output=True, text=True, encoding='utf-8')
        # Limit printing long stdout unless necessary (like makeblastdb info)
        if len(process.stdout) < 2000 or "BLAST Database creation" in process.stdout:
             print("Stdout:", process.stdout)
        else: print("[stdout is long, truncated]")
        # Always print stderr unless it's just expected makeblastdb info
        if process.stderr and not ("BLAST Database creation" in process.stderr and "INFO" in process.stderr):
                 print("Stderr:", process.stderr, file=sys.stderr)
        return process
    except FileNotFoundError as e:
        print(f"Error: Command not found - {e}. Is BLAST+ installed and in your system's PATH?")
        sys.exit(1)
    except subprocess.CalledProcessError as e:
        print(f"Error executing command: {' '.join(command_list)}")
        print(f"Return code: {e.returncode}")
        print(f"Output (first 1000 chars):\n{e.output[:1000] if e.output else 'None'}")
        print(f"Stderr (first 1000 chars):\n{e.stderr[:1000] if e.stderr else 'None'}")
        sys.exit(1)
    except Exception as e:
        print(f"An unexpected error occurred while running command: {e}")
        sys.exit(1)

# --- Function to sort BLAST results ---
def sort_blast_results(input_file, output_file, column_index):
    """Reads BLAST TSV output, sorts by a specific column (numeric), and writes sorted output."""
    print(f"Sorting BLAST results from '{input_file}' by column {column_index+1} (Bitscore)...")
    try:
        with open(input_file, 'r', encoding='utf-8') as infile:
            lines = infile.readlines()
        header_lines = [line for line in lines if line.startswith('#')] # Keep potential headers
        data_lines = [line for line in lines if not line.startswith('#') and line.strip()] # Ignore empty lines
        if not data_lines:
            print(f"No data lines found in BLAST output '{input_file}' for sorting.")
            with open(output_file, 'w', encoding='utf-8') as outfile: outfile.writelines(header_lines)
            print(f"Empty sorted file (with headers, if any) written to: {output_file}")
            return True
        parsed_data = []
        line_offset = len(header_lines)
        for i, line in enumerate(data_lines):
            fields = line.strip().split('\t')
            if len(fields) > column_index:
                try:
                    score = float(fields[column_index])
                    parsed_data.append((score, line))
                except ValueError:
                    print(f"Warning (Sort): Could not parse score in column {column_index+1} on line {i+line_offset+1}. Skipping: {line.strip()}", file=sys.stderr)
            else:
                 print(f"Warning (Sort): Line {i+line_offset+1} has fewer columns than expected ({len(fields)}). Skipping: {line.strip()}", file=sys.stderr)
        # Sort in descending order (highest score first)
        parsed_data.sort(key=lambda x: x[0], reverse=True)
        with open(output_file, 'w', encoding='utf-8') as outfile:
            outfile.writelines(header_lines)
            for score, line in parsed_data: outfile.write(line) # line already includes newline
        print(f"Sorted BLAST results written to: {output_file}")
        return True
    except FileNotFoundError:
        print(f"Error: Cannot sort. Input BLAST file not found: '{input_file}'")
        return False
    except Exception as e:
        print(f"An error occurred during BLAST result sorting for '{input_file}': {e}")
        return False

# --- Function to filter BLAST results by removing self-hits (qseqid == sseqid) ---
def filter_blast_results_self_hits(input_file, output_file, qid_idx, sid_idx):
    """
    Reads a sorted BLAST TSV file and filters out hits where the query sequence ID
    is identical to the subject sequence ID.
    """
    print(f"Filtering BLAST results from '{input_file}' to remove self-hits (qseqid == sseqid)...")
    lines_written = 0; lines_read = 0; lines_filtered_out = 0
    try:
        with open(input_file, 'r', encoding='utf-8') as infile, \
             open(output_file, 'w', encoding='utf-8') as outfile:
            for i, line in enumerate(infile):
                lines_read += 1
                if line.startswith('#'): # Keep header lines
                    outfile.write(line)
                    lines_written +=1
                    continue
                if not line.strip(): # Skip empty lines
                     continue
                fields = line.strip().split('\t')
                if len(fields) > max(qid_idx, sid_idx):
                    query_id = fields[qid_idx]
                    subject_id = fields[sid_idx]
                    # Extract only the sequence ID part before any potential '|' or space
                    # Adjust this logic if your IDs have different structures
                    query_id_base = query_id.split('|')[0].split(' ')[0]
                    subject_id_base = subject_id.split('|')[0].split(' ')[0]
                    if query_id_base != subject_id_base:
                        outfile.write(line)
                        lines_written += 1
                    else:
                        lines_filtered_out += 1
                else:
                    print(f"Warning (Filter): Line {i+1} has fewer columns ({len(fields)}) than needed for qseqid/sseqid. Keeping line: {line.strip()}", file=sys.stderr)
                    outfile.write(line) # Keep lines that can't be parsed properly
                    lines_written += 1
        print(f"Filtering complete. Read {lines_read} lines, kept {lines_written} lines.")
        if lines_filtered_out > 0: print(f"Filtered out {lines_filtered_out} self-hits where qseqid base == sseqid base.")
        print(f"Filtered BLAST results written to: {output_file}")
        return True
    except FileNotFoundError:
        print(f"Error: Cannot filter. Input sorted BLAST file not found: '{input_file}'")
        return False
    except Exception as e:
        print(f"An error occurred during BLAST result filtering for '{input_file}': {e}")
        return False
# --- End of Function Definitions ---


# --- Main Workflow ---

# 1. Validate Inputs & Create Output Directory
print(f"Creating main output directory: {main_output_dir}")
os.makedirs(main_output_dir, exist_ok=True)

# Check required files
required_files = {
    "Source Query File": source_query_file,
    "Mouse Peptide File": os.path.join(data_dir, mouse_pep_file),
    "Human Peptide File": os.path.join(data_dir, human_pep_file)
}
abort = False
for name, path in required_files.items():
    if not os.path.exists(path):
        print(f"Error: {name} not found at '{path}'")
        abort = True
    elif os.path.getsize(path) == 0:
        print(f"Error: {name} at '{path}' is empty.")
        abort = True
if abort: sys.exit(1)

# Check if data_dir exists
if not os.path.isdir(data_dir):
    print(f"Error: Data directory '{data_dir}' not found.")
    sys.exit(1)


# 2. Prepare Query Sequences (Split by Species and Gene)
print("\nStep 1: Reading source query file and splitting by species and gene symbol...")
sequences_by_species_gene = {
    "Mouse": defaultdict(list),
    "Human": defaultdict(list)
}
source_records_count = 0
mouse_records_count = 0
human_records_count = 0
other_records_count = 0
parsing_warnings = 0

try:
    for record in SeqIO.parse(source_query_file, "fasta"):
        source_records_count += 1
        header = record.description
        seq_id = record.id # Use ID for species check (more reliable for Ensembl)

        # Determine Species (using Ensembl peptide ID prefix convention)
        origin_species = None
        if seq_id.startswith("ENSMUSP"): # Mouse peptide ID
            origin_species = "Mouse"
            mouse_records_count += 1
        elif seq_id.startswith("ENSP"): # Human peptide ID
            origin_species = "Human"
            human_records_count += 1
        else:
            other_records_count += 1
            # Optional: Try parsing species from header if needed, or just warn/skip
            # print(f"Warning: Could not determine species for query '{seq_id}' based on ID prefix. Skipping.")
            # parsing_warnings += 1
            # continue # Skip if species cannot be determined reliably

        # Extract Gene Symbol
        gene_symbol = None
        match = re.search(r'gene_symbol:(\S+)', header)
        if match:
            gene_symbol = match.group(1)
        else:
            print(f"Warning: Could not parse 'gene_symbol:' from header for query '{seq_id}' (Species: {origin_species}). Skipping.")
            parsing_warnings += 1
            continue # Skip if gene symbol cannot be determined

        # Store the sequence if species and gene symbol were found
        if origin_species and gene_symbol:
             sequences_by_species_gene[origin_species][gene_symbol].append(record)

except Exception as e:
    print(f"An error occurred while reading or parsing the source query file '{source_query_file}': {e}")
    sys.exit(1)

print(f"Parsed {source_records_count} sequences from source query file.")
print(f"  - Identified {mouse_records_count} Mouse sequences across {len(sequences_by_species_gene['Mouse'])} gene symbols.")
print(f"  - Identified {human_records_count} Human sequences across {len(sequences_by_species_gene['Human'])} gene symbols.")
if other_records_count > 0: print(f"  - Found {other_records_count} sequences with unrecognized species ID prefixes.")
if parsing_warnings > 0: print(f"  - Encountered {parsing_warnings} parsing warnings (see above).")


# 3. Loop Through Each Species Run (Mouse vs Non-Mouse, Human vs Non-Human)
species_runs = [
    {"query_species": "Mouse", "exclude_file": mouse_pep_file, "db_name": "non_mouse"},
    {"query_species": "Human", "exclude_file": human_pep_file, "db_name": "non_human"}
]

for run_config in species_runs:
    query_species = run_config["query_species"]
    exclude_filename = run_config["exclude_file"]
    db_suffix = run_config["db_name"]

    print(f"\n{'='*20} Starting Run: {query_species} Queries vs. {db_suffix} Database {'='*20}")

    # Define paths specific to this run
    run_output_dir = os.path.join(main_output_dir, f"{query_species.lower()}_vs_{db_suffix}")
    print(f"Creating output directory for this run: {run_output_dir}")
    os.makedirs(run_output_dir, exist_ok=True)

    combined_fasta_db_file = os.path.join(run_output_dir, f"{db_suffix}_species.pep.fa")
    blast_db_name = os.path.join(run_output_dir, f"{db_suffix}_species_db") # Full path for db

    # --- Database Creation Steps (for this run) ---
    print(f"\nStep 2 ({query_species} run): Combining FASTA files for {db_suffix} database...")
    all_files = glob.glob(os.path.join(data_dir, "*.fa")) + glob.glob(os.path.join(data_dir, "*.fasta"))
    exclude_file_path = os.path.abspath(os.path.join(data_dir, exclude_filename))
    files_to_combine = [f for f in all_files if os.path.abspath(f) != exclude_file_path]

    if not files_to_combine:
        print(f"Error: No FASTA files found in '{data_dir}' to combine (excluding '{exclude_filename}'). Skipping this run.")
        continue # Skip to the next species run

    print(f"Found {len(files_to_combine)} files to combine for the {db_suffix} database:")
    # for f in files_to_combine: print(f"  - {os.path.basename(f)}") # Can be verbose

    cat_command = ["cat"] + files_to_combine
    try:
        with open(combined_fasta_db_file, "w", encoding='utf-8') as outfile:
            # Use run_command for consistency, though direct subprocess might be slightly simpler here
            process = subprocess.run(cat_command, check=True, stdout=outfile, text=True, encoding='utf-8', stderr=subprocess.PIPE)
            if process.stderr: print("Stderr from cat:", process.stderr, file=sys.stderr)
        print(f"Successfully combined files into: {combined_fasta_db_file}")
    except Exception as e:
        print(f"Error combining files for {db_suffix} database: {e}. Skipping this run.")
        continue

    print(f"\nStep 3 ({query_species} run): Building BLAST database '{os.path.basename(blast_db_name)}'...")
    makeblastdb_command = [
        "makeblastdb", "-in", combined_fasta_db_file, "-dbtype", "prot",
        "-out", blast_db_name, # Use full path
        "-title", f"Combined {db_suffix} Peptide Database",
        "-parse_seqids"
    ]
    run_command(makeblastdb_command) # cwd defaults to current dir, output paths are absolute
    print(f"BLAST database '{os.path.basename(blast_db_name)}' created in '{run_output_dir}'.")


    # --- Per-Gene BLAST Steps (for this run) ---
    print(f"\nStep 4 ({query_species} run): Preparing query files and running BLASTp...")
    gene_symbols_to_process = sequences_by_species_gene[query_species]
    if not gene_symbols_to_process:
         print(f"No query sequences identified for {query_species}. Skipping BLAST steps for this run.")
         continue

    print(f"Found {len(gene_symbols_to_process)} gene symbols with {query_species} sequences to process.")
    temp_query_files = {} # Store paths for cleanup
    blast_success_count = 0
    sort_success_count = 0
    filter_success_count = 0
    genes_processed_count = 0

    for gene_symbol, sequences in sorted(gene_symbols_to_process.items()):
        if not sequences: continue # Should not happen with defaultdict(list)

        genes_processed_count += 1
        safe_gene_symbol = re.sub(r'[\\/*?:"<>|]', '_', gene_symbol) # Sanitize for filename
        print(f"\n  --- Processing Gene: {gene_symbol} ({query_species} query) ---")

        # Prepare temporary query file for this gene
        temp_file_path = os.path.join(run_output_dir, f"temp_query_{safe_gene_symbol}_{query_species.lower()}.fa")
        try:
            SeqIO.write(sequences, temp_file_path, "fasta")
            temp_query_files[gene_symbol] = temp_file_path # Store for cleanup
            # print(f"    Created temp query file: {os.path.basename(temp_file_path)}")
        except Exception as e:
            print(f"    Error writing temporary query file for gene '{gene_symbol}': {e}. Skipping this gene.")
            continue

        # Define output file paths for this gene's BLAST run
        raw_output_file = os.path.join(run_output_dir, f"{safe_gene_symbol}.blastp.raw.tsv")
        sorted_output_file = os.path.join(run_output_dir, f"{safe_gene_symbol}.blastp.sorted.tsv")
        filtered_output_file = os.path.join(run_output_dir, f"{safe_gene_symbol}.blastp.filtered_noself.tsv")

        # Run BLASTp
        print(f"    Running BLASTp against {db_suffix} database...")
        blastp_command = [
            "blastp", "-query", temp_file_path, "-db", blast_db_name, # Use full path to db
            "-out", raw_output_file, "-num_threads", str(num_threads),
            "-evalue", evalue_threshold, "-max_target_seqs", max_target_seqs,
            "-outfmt", outfmt_string
        ]
        # Run command, don't exit script on failure (check=False), just report and continue
        process = run_command(blastp_command, check=False)

        if process.returncode != 0:
             print(f"    BLASTp for gene '{gene_symbol}' finished with errors (exit code {process.returncode}). Check messages. Raw output (if any): {os.path.basename(raw_output_file)}")
             continue # Skip sorting/filtering for this gene

        print(f"    BLASTp completed. Raw results: {os.path.basename(raw_output_file)}")
        blast_success_count += 1

        # Sort Results
        if sort_blast_results(raw_output_file, sorted_output_file, bitsore_column_index):
            sort_success_count += 1
            # Filter Results (Self-Hits - still useful for multi-isoform queries)
            if filter_blast_results_self_hits(sorted_output_file, filtered_output_file, qseqid_column_index, sseqid_column_index):
                filter_success_count += 1
                print(f"    Sorting and Filtering successful. Final results: {os.path.basename(filtered_output_file)}")
            else:
                print(f"    Filtering (self-hits) failed for '{gene_symbol}'. Check logs.")
        else:
            print(f"    Sorting failed for '{gene_symbol}'. Check logs. Skipping filtering.")

    # --- Cleanup Temporary Query Files (for this run) ---
    print(f"\nStep 5 ({query_species} run): Cleaning up temporary query files...")
    for gene_symbol, temp_file_path in temp_query_files.items():
        try:
            if os.path.exists(temp_file_path):
                 os.remove(temp_file_path)
                 # print(f"  - Removed: {os.path.basename(temp_file_path)}")
            # else: # Already removed or failed to create
                 # print(f"  - Already removed or not created: {os.path.basename(temp_file_path)}")
                 pass
        except OSError as e:
            print(f"Warning: Could not remove temporary file '{temp_file_path}': {e}")

    # --- Summary for this run ---
    print(f"\n--- Summary for {query_species} Queries vs. {db_suffix} Database ---")
    print(f"Attempted to process {genes_processed_count} gene symbols.")
    print(f"Successfully completed BLASTp for {blast_success_count} gene symbol(s).")
    if blast_success_count > 0:
        print(f"Successfully sorted results for {sort_success_count} gene symbol(s).")
        if sort_success_count > 0:
             print(f"Successfully filtered results (removing self-hits) for {filter_success_count} gene symbol(s).")
    print(f"Results for this run are in: {run_output_dir}")


# --- Final Overall Summary ---
print(f"\n{'='*20} Overall Workflow Completed {'='*20}")
print(f"Processed runs for: {', '.join(run['query_species'] for run in species_runs)}")
print(f"All results are organized within the main directory: {main_output_dir}")
