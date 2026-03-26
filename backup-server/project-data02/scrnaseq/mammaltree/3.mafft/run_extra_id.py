#!/usr/bin/env python3

import os
import sys
import re # For filename sanitization
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord # Although we copy the whole record

# --- Configuration ---

# 1. Define the sequences to extract per group
#    Key: Name for the output FASTA file (without extension)
#    Value: List of exact sequence IDs to extract
sequences_to_extract = {
    'Cavin3_Analysis': ['ENSGALP00010014540.1', 'ENSMUSP00000044979.3', 'ENSP00000305675.4',
                        'ENSP00000307292.3', 'ENSP00000432047.1', 'ENSPMRP00000005185.1',
                        'ENSXETP00000107494.1'],
    'Dsg3_Analysis': ['ENSGALP00010005873.1', 'ENSMUSP00000064718.7', 'ENSMUSP00000157289.2',
                      'ENSP00000257189.4', 'ENSP00000261590.8', 'ENSP00000311859.4',
                      'ENSP00000352785.4', 'ENSP00000519121.1', 'ENSP00000519123.1',
                      'ENSP00000519125.1', 'ENSP00000519126.1', 'ENSXETP00000074721.2'],
    'Lalba_Analysis': ['ENSMUSP00000023726.4', 'ENSP00000301046.2', 'ENSP00000449780.1',
                       'ENSPMRP00000034493.1'],
    'Pip_Analysis': ['ENSMUSP00000074818.8', 'ENSMUSP00000144995.2', 'ENSP00000291009.3',
                     'ENSXETP00000037866.5', 'ENSXETP00000117478.1'],
    'Plet1_Analysis': ['ENSMUSP00000110118.2', 'ENSMUSP00000139422.2', 'ENSP00000341412.2']
}

# 2. List of Paths to the source peptide FASTA files
#    The script will search for IDs in ALL files listed here.
#    *** Updated based on your request ***
source_fasta_files = [
    "/data02/sunxuebo/project/scrnaseq/mammaltree/1.extra_pep/all_matching_peptides_human_mouse_separate_lists.fa",
    "/data02/sunxuebo/project/scrnaseq/mammaltree/2.blastp/blast_results_per_gene_noself/combined_species.pep.fa"
]

# 3. Output Directory for the extracted FASTA files
#    *** MODIFY AS NEEDED ***
output_dir = "./extracted_sequences_by_group_combined_source"

# --------------------- End of Configuration ---------------------

# --- Main Script ---

print("--- Starting Sequence Extraction by Group from Multiple Sources ---")

# 1. Validate Inputs
print("\n--- 1. Validating Inputs ---")
valid_source_files = []
for f_path in source_fasta_files:
    if not os.path.exists(f_path):
        print(f"Error: Source FASTA file not found: {f_path}. Cannot proceed.")
        # Exit if a source file is missing, as all IDs might not be found
        sys.exit(1)
    else:
        valid_source_files.append(f_path)

if not valid_source_files:
     print("Error: No source FASTA files provided or found.")
     sys.exit(1)

if not sequences_to_extract:
    print("Error: The 'sequences_to_extract' dictionary is empty. Nothing to do.")
    sys.exit(1)

# Create output directory
try:
    os.makedirs(output_dir, exist_ok=True)
    print(f"Output directory: {os.path.abspath(output_dir)}")
except OSError as e:
    print(f"Error: Could not create output directory '{output_dir}': {e}")
    sys.exit(1)

# 2. Index the Source FASTA Files
print("\n--- 2. Indexing Source Files ---")
source_indices = [] # Store index objects in a list
total_indexed_sequences = 0
for f_path in valid_source_files:
    print(f"Indexing: {os.path.basename(f_path)} ...")
    try:
        index = SeqIO.index(f_path, "fasta")
        count = len(index)
        if count == 0:
            print(f"  Warning: Source FASTA file '{f_path}' appears to be empty.")
            # Continue indexing other files, but the missing IDs might come from here
        else:
            print(f"  Indexed {count} sequences.")
            source_indices.append(index)
            total_indexed_sequences += count
    except Exception as e:
        print(f"  Error indexing source file {f_path}: {e}")
        print("  Cannot continue without indexing all source files.")
        sys.exit(1)

if not source_indices:
    print("Error: Failed to index any source FASTA files.")
    sys.exit(1)
print(f"Total sequences indexed across {len(source_indices)} file(s): {total_indexed_sequences}")

# 3. Extract Sequences for Each Group
print("\n--- 3. Extracting Sequences per Group ---")
groups_processed = 0
groups_with_missing_ids = 0

for group_name, id_list in sequences_to_extract.items():
    print(f"\n>>> Processing Group: {group_name}")
    sequences_for_this_group = []
    missing_ids_in_group = []
    found_ids_this_group = set() # Track IDs found to avoid duplicates if present in multiple files

    if not id_list:
        print("    Warning: ID list for this group is empty. Skipping.")
        continue

    # Sanitize group name for filename
    safe_group_name = re.sub(r'[\\/*?:"<>| ]', '_', group_name)
    output_filepath = os.path.join(output_dir, f"{safe_group_name}.fa")

    # Iterate through requested IDs for the current group
    for seq_id in id_list:
        record_found_in_any_file = None
        # Search for the ID in all indexed source files
        for index in source_indices:
            if seq_id in index:
                try:
                    # Retrieve the full SeqRecord
                    record_found_in_any_file = index[seq_id]
                    # If we found it, we can stop searching other files for this ID
                    break
                except Exception as e:
                    # Log if retrieval fails even if ID is listed in index
                    print(f"    Warning: Error retrieving '{seq_id}' from an index where it was listed: {e}")
                    # Continue searching other indices just in case
            # else: ID not in this specific index

        # Process the found record or note the missing ID
        if record_found_in_any_file:
            # Check if we already added this ID (e.g., if it was in multiple source files)
            if seq_id not in found_ids_this_group:
                 sequences_for_this_group.append(record_found_in_any_file)
                 found_ids_this_group.add(seq_id)
            # else: Already added this ID from another file, skip duplicate
        else:
            # ID was not found in *any* of the source indices
            missing_ids_in_group.append(seq_id)

    # Report missing IDs for this group
    if missing_ids_in_group:
        groups_with_missing_ids += 1
        print(f"    Warning: Could not find the following IDs for group '{group_name}' in any source file: {', '.join(missing_ids_in_group)}")

    # Write the FASTA file for the group if sequences were found
    if sequences_for_this_group:
        try:
            with open(output_filepath, "w") as outfile:
                SeqIO.write(sequences_for_this_group, outfile, "fasta")
            print(f"    Successfully extracted {len(sequences_for_this_group)} sequences to: {os.path.basename(output_filepath)}")
            groups_processed += 1
        except IOError as e:
            print(f"    Error writing output file '{output_filepath}': {e}")
    else:
        print(f"    No sequences were found or retrieved for group '{group_name}'. No output file created.")


# 4. Final Summary
print("\n--- 4. Extraction Summary ---")
total_groups = len(sequences_to_extract)
print(f"Attempted to process {total_groups} groups.")
print(f"Successfully created output files for {groups_processed} groups.")
if groups_with_missing_ids > 0:
     print(f"Warning: {groups_with_missing_ids} group(s) had one or more missing sequence IDs (see warnings above).")
if groups_processed == 0 and total_groups > 0:
    print("No output files were generated. Check source file content and requested IDs.")

print(f"\nExtraction process finished.")
print(f"Output FASTA files are located in: {os.path.abspath(output_dir)}")

# Close index file handles if using 'r' mode (not strictly necessary with default in-memory index)
# for index in source_indices:
#     if hasattr(index, 'close'):
#          index.close()
