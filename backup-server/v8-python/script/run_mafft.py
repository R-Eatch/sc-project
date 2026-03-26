#!/usr/bin/env python3

import os
import subprocess
import glob
import csv
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
import sys
import shutil

# --- Configuration ---
# *** NEW: Global flag to control selection logic ***
SELECT_BEST_PER_SPECIES = True  # Set to True to enable best-hit-per-species, False for original top N unique IDs

# Directory containing BLASTp results (.tsv files) and the combined FASTA DB
blast_dir = "../2.blastp/"
# Directory containing the original query sequences
query_dir = "../1.extra_pep/"
# Output directory for MAFFT alignments
mafft_dir = "./"  # Output to current directory

# Input file names
combined_db_fasta = os.path.join(blast_dir, "combined_species.pep.fa")
query_fasta = os.path.join(query_dir, "longest_queries.fa")

# BLAST result file pattern
blast_result_pattern = os.path.join(blast_dir, "*.blastp.tsv")

# Number of top hits OR top species representatives to select
top_n_hits = 10 # Max number of sequences (either unique IDs or species representatives)

# MAFFT Parameters
mafft_threads = 20
mafft_options = ["--auto"]

# Column indices in the BLAST TSV output
QSEQID_COL = 0
SSEQID_COL = 1
BITSCORE_COL = 11

# Dictionary for header modification
species_to_prefix_map = {
    'Chicken': 'ENSGALP',
    'Opossum': 'ENSMODP',
    'Human': 'ENSP',
    'Mouse': 'ENSMUSP',
    'Lizard': 'ENSPMRP',
    'Frog': 'ENSXETP'
}
prefix_to_species_map = {v: k for k, v in species_to_prefix_map.items()}
# ---------------------

# --- Helper Function to Get Species Name from ID ---
def get_species_from_id(sequence_id, prefix_map):
    """Tries to determine species name based on known prefixes."""
    for prefix, name in prefix_map.items():
        if sequence_id.startswith(prefix):
            return name
    return None # Return None if no known prefix matches

# --- Function to run shell commands ---
# (Keep the run_command function as it was)
def run_command(command_list, check=True, capture_output=False, text=True, stdout_file=None):
    """Runs a shell command using subprocess."""
    print(f"Running command: {' '.join(command_list)}")
    try:
        process = subprocess.run(command_list, check=check,
                                 capture_output=capture_output, text=text,
                                 stdout=stdout_file, stderr=subprocess.PIPE) # Capture stderr
        if capture_output and process.stdout:
            print(process.stdout)
        if process.stderr:
            print("Stderr:", process.stderr, file=sys.stderr)
        return process
    except FileNotFoundError as e:
        print(f"Error: Command not found - {e}. Is MAFFT installed and in PATH?")
        sys.exit(1)
    except subprocess.CalledProcessError as e:
        print(f"Error executing command: {' '.join(command_list)}")
        print(f"Return code: {e.returncode}")
        if capture_output: print(f"Output:\n{e.output}")
        print(f"Stderr:\n{e.stderr}")
        sys.exit(1)
    except Exception as e:
        print(f"An unexpected error occurred: {e}")
        sys.exit(1)

# --- Function to modify headers ---
# (Keep the modify_fasta_headers function as it was)
def modify_fasta_headers(input_fasta, output_fasta, prefix_map):
    """Reads a FASTA file, modifies headers based on prefix map, writes to new file."""
    print(f"  Modifying headers for: {os.path.basename(input_fasta)}")
    modified_records = []
    records_processed = 0
    headers_modified_count = 0
    
    try:
        for record in SeqIO.parse(input_fasta, "fasta"):
            records_processed += 1
            # Use the part before '|' if header was already modified, else full ID
            original_id = record.id.split('|')[0] 
            matched_prefix = None
            
            for prefix in prefix_map.keys():
                if original_id.startswith(prefix):
                    matched_prefix = prefix
                    break 

            if matched_prefix:
                species_name = prefix_map[matched_prefix]
                # Only add species tag if it's not already there
                if f"|{species_name}" not in record.id:
                    new_id = f"{original_id}|{species_name}"
                    
                    if record.description.startswith(original_id):
                        desc_start_index = len(original_id)
                        if desc_start_index < len(record.description) and record.description[desc_start_index] == ' ':
                             desc_start_index += 1
                        description_text = record.description[desc_start_index:]
                    else: 
                        description_text = record.description 

                    new_record = SeqRecord(record.seq, id=new_id, description=description_text, name='') 
                    modified_records.append(new_record)
                    headers_modified_count += 1
                else:
                     # Already has the correct tag, keep as is
                     modified_records.append(record)
            else:
                modified_records.append(record)

        with open(output_fasta, "w") as outfile:
            SeqIO.write(modified_records, outfile, "fasta")
        
        print(f"  Finished modifying headers. {headers_modified_count}/{records_processed} headers updated/verified.")
        return True 

    except Exception as e:
        print(f"Error modifying headers for {input_fasta}: {e}")
        return False

# --- Input Validation ---
# (Keep input validation as it was)
if shutil.which("mafft") is None: exit("Error: MAFFT not found.")
if not os.path.isdir(blast_dir): exit(f"Error: BLAST directory not found: {blast_dir}")
if not os.path.isdir(query_dir): exit(f"Error: Query directory not found: {query_dir}")
if not os.path.exists(combined_db_fasta): exit(f"Error: Combined DB FASTA not found: {combined_db_fasta}")
if not os.path.exists(query_fasta): exit(f"Error: Query FASTA not found: {query_fasta}")

# --- Main Workflow ---
print(f"Ensuring output directory exists: {os.path.abspath(mafft_dir)}")
os.makedirs(mafft_dir, exist_ok=True)

print(f"\nIndexing database FASTA: {combined_db_fasta} ...")
try: db_index = SeqIO.index(combined_db_fasta, "fasta")
except Exception as e: exit(f"Error indexing database FASTA: {e}")
print(f"Indexed {len(db_index)} sequences.")

print(f"Indexing query FASTA: {query_fasta} ...")
try: query_index = SeqIO.index(query_fasta, "fasta")
except Exception as e: exit(f"Error indexing query FASTA: {e}")
print(f"Indexed {len(query_index)} query sequences.")

blast_result_files = glob.glob(blast_result_pattern)
if not blast_result_files: exit(f"Error: No BLAST result files found: {blast_result_pattern}")

print(f"\nFound {len(blast_result_files)} BLAST result files to process.")
print(f"Selection mode: {'Best Hit Per Species' if SELECT_BEST_PER_SPECIES else 'Top N Unique Subject IDs'}") # Inform user

for blast_file in blast_result_files:
    base_name = os.path.basename(blast_file).split('.')[0]
    print(f"\nProcessing: {base_name} (File: {blast_file})")

    # --- Read ALL BLAST results ---
    hits = []
    query_id_for_this_file = None
    try:
        with open(blast_file, 'r', newline='') as tsvfile:
            reader = csv.reader(tsvfile, delimiter='\t')
            for row_num, row in enumerate(reader): # Add row number for better warnings
                if not row or row[0].startswith("#"): continue
                try:
                    if len(row) > max(QSEQID_COL, SSEQID_COL, BITSCORE_COL):
                        qseqid = row[QSEQID_COL]
                        sseqid = row[SSEQID_COL]
                        bitscore = float(row[BITSCORE_COL])
                        hits.append({'qseqid': qseqid, 'sseqid': sseqid, 'bitscore': bitscore})
                        if query_id_for_this_file is None: query_id_for_this_file = qseqid
                    else: print(f"Warning: Skipping malformed row {row_num+1} in {blast_file}: {row}", file=sys.stderr)
                except (ValueError, IndexError) as e: print(f"Warning: Error parsing row {row_num+1} in {blast_file}: {row} - {e}", file=sys.stderr)
    except FileNotFoundError: print(f"Error: BLAST result file not found: {blast_file}"); continue
    except Exception as e: print(f"Error reading BLAST file {blast_file}: {e}"); continue

    if not hits: print(f"  No valid hits found. Skipping."); continue
    if query_id_for_this_file is None: print(f"  Could not determine query ID. Skipping."); continue

    # Sort ALL hits by bitscore
    hits.sort(key=lambda x: x['bitscore'], reverse=True)

    # --- Apply Selection Logic ---
    final_subject_ids = []
    if SELECT_BEST_PER_SPECIES:
        print(f"  Applying 'Best Hit Per Species' logic (max {top_n_hits} species)...")
        species_already_selected = set()
        for hit in hits:
            sseqid = hit['sseqid']
            species_name = get_species_from_id(sseqid, prefix_to_species_map)

            if species_name is None:
                # Optional: Decide how to handle hits whose species prefix isn't known
                # print(f"    Warning: Unknown species prefix for {sseqid}. Skipping.")
                continue # Skip hits from unknown species based on prefix map

            if species_name not in species_already_selected:
                final_subject_ids.append(sseqid)
                species_already_selected.add(species_name)
                # Stop if we have reached the desired number of species representatives
                #if len(final_subject_ids) >= top_n_hits:
                #    break
        #print(f"  Selected {len(final_subject_ids)} hits representing unique species.")

    else: # Original logic: Top N unique subject IDs
        print(f"  Applying 'Top N Unique Subject IDs' logic (max {top_n_hits} IDs)...")
        seen_subjects = set()
        for hit in hits:
            sseqid = hit['sseqid']
            if sseqid not in seen_subjects:
                final_subject_ids.append(sseqid)
                seen_subjects.add(sseqid)
                if len(final_subject_ids) >= top_n_hits:
                    break
        print(f"  Selected {len(final_subject_ids)} top unique subject IDs.")

    # --- Proceed with selected IDs ---
    if not final_subject_ids:
        print(f"  No subject IDs selected after applying filter. Skipping alignment.")
        continue

    print(f"  Final subject IDs selected for alignment: {len(final_subject_ids)}")

    # --- Retrieve sequences ---
    sequences_to_align = []
    try:
        query_record = query_index[query_id_for_this_file]
        sequences_to_align.append(query_record)
        print(f"  Retrieved query sequence: {query_id_for_this_file}")
    except KeyError:
        print(f"Error: Query ID '{query_id_for_this_file}' not found. Skipping {base_name}.")
        continue

    missing_subjects = []
    retrieved_count = 0
    for sseqid in final_subject_ids: # Use the final selected list
        try:
            subject_record = db_index[sseqid]
            sequences_to_align.append(subject_record)
            retrieved_count += 1
        except KeyError:
            print(f"Warning: Subject ID '{sseqid}' not found in DB index. Excluded.")
            missing_subjects.append(sseqid)
    print(f"  Retrieved {retrieved_count} subject sequences.")

    if len(sequences_to_align) < 2:
         print(f"  Only {len(sequences_to_align)} sequence(s) total (need >= 2). Skipping MAFFT.")
         continue

    # --- Prepare for MAFFT ---
    temp_fasta_file = os.path.join(mafft_dir, f"{base_name}_temp_input.fa")
    alignment_output_file = os.path.join(mafft_dir, f"{base_name}.aln.fa")
    temp_modified_header_file = alignment_output_file + ".mod_tmp"

    # --- Main processing block with finally for cleanup ---
    header_mod_success = False # Initialize outside try/finally
    try:
        print(f"  Preparing temporary input FASTA: {temp_fasta_file}")
        try:
            with open(temp_fasta_file, "w") as temp_out:
                SeqIO.write(sequences_to_align, temp_out, "fasta")
        except IOError as e:
            print(f"Error writing temporary FASTA file {temp_fasta_file}: {e}")
            raise # Propagate error to outer handler

        # --- Run MAFFT ---
        print(f"  Running MAFFT (Output: {alignment_output_file})...")
        mafft_command = ["mafft", f"--thread", str(mafft_threads)] + mafft_options + [temp_fasta_file]

        mafft_success = False
        try:
            with open(alignment_output_file, "w") as outfile:
                run_command(mafft_command, stdout_file=outfile)
            print(f"  MAFFT alignment raw output saved to: {alignment_output_file}")
            mafft_success = True
        except Exception as e:
            print(f"  Error occurred during MAFFT execution for {base_name}: {e}")
            if os.path.exists(alignment_output_file):
                try: os.remove(alignment_output_file)
                except OSError: pass
            # Decide if MAFFT error should stop processing this gene; currently it continues to finally

        # --- Modify Headers (if MAFFT succeeded) ---
        if mafft_success and os.path.exists(alignment_output_file):
            try:
                if modify_fasta_headers(alignment_output_file, temp_modified_header_file, prefix_to_species_map):
                    os.replace(temp_modified_header_file, alignment_output_file)
                    print(f"  Successfully replaced alignment file with modified headers: {alignment_output_file}")
                    header_mod_success = True # Mark success only after replace
                else:
                    print(f"  Header modification function failed. Keeping original MAFFT output.")
                    if os.path.exists(temp_modified_header_file):
                         try: os.remove(temp_modified_header_file)
                         except OSError: pass
            except Exception as e: # Catch errors during modify/replace
                 print(f"Error during header modification or replacement: {e}")
                 if os.path.exists(temp_modified_header_file):
                      try: os.remove(temp_modified_header_file)
                      except OSError: pass

    # --- Outer exception handling ---
    except IOError: # Handled specific error from writing temp file
        print(f"  Skipping {base_name} due to error writing input file.")
    except Exception as e_outer: # Catch any other unexpected errors in the main try block
        print(f"  An unexpected error occurred processing {base_name}: {e_outer}")

    # --- Finally block for cleanup ---
    finally:
        print(f"  Cleaning up temporary files for {base_name}...")
        # Cleanup input file
        if os.path.exists(temp_fasta_file):
            try: os.remove(temp_fasta_file)
            except OSError as e: print(f"Warning: Could not remove {temp_fasta_file}: {e}")
        # Cleanup modification temp file only if modification/replacement didn't fully succeed
        if not header_mod_success and os.path.exists(temp_modified_header_file):
            try: os.remove(temp_modified_header_file)
            except OSError as e: print(f"Warning: Could not remove {temp_modified_header_file}: {e}")

# End of loop for blast_file

print("\nMAFFT alignment and header modification process completed.")
print(f"Aligned files are located in: {os.path.abspath(mafft_dir)}")
