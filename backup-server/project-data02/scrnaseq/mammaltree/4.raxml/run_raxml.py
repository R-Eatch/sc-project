#!/usr/bin/env python3

import os
import subprocess
import glob
import sys
import shutil # For checking executable existence

# --- Configuration ---
# Directory containing MAFFT alignment files (.aln.fa)
mafft_dir = "../3.mafft/"
# Output directory for RAxML v8 results
raxml_v8_output_dir = "./" 

# RAxML v8 Parameters
# *** ADJUST AS NEEDED ***

# !!! VERY IMPORTANT: Check the actual executable name installed by Conda !!!
# Common names: raxmlHPC-PTHREADS, raxmlHPC-SSE3, raxmlHPC-AVX, raxmlHPC, raxml
raxml_executable = "raxmlHPC-PTHREADS" 

num_threads = 20  # Number of CPU cores for RAxML (-T option)
bootstrap_reps = 500 # Number of bootstrap replicates (-# or -N option)
# Standard RAxML often requires two seeds for the -f a option:
random_seed_bootstrap = 12345 # Seed for bootstrapping (-x option)
random_seed_parsimony = 54321 # Seed for parsimony starting trees (-p option)

# Substitution model for standard RAxML (protein + Gamma + LG + empirical Freqs)
# Common codes: PROTGAMMALGF, PROTGAMMAWAGF, PROTGAMMAJTTF etc.
model = "PROTGAMMALGF" 

# Input alignment file pattern
alignment_pattern = os.path.join(mafft_dir, "*.aln.fa")
# ---------------------

# --- Function to run shell commands ---
def run_command(command_list, cwd=None, check=True):
    """Runs a shell command using subprocess and checks for errors."""
    print(f"Running command: {' '.join(command_list)}")
    try:
        # Standard RAxML can be very verbose, capture might fill memory for long runs
        # Let it print directly to console unless debugging is needed.
        process = subprocess.run(command_list, cwd=cwd, check=check, 
                                 text=True, stderr=subprocess.PIPE) # Capture stderr
        print(f"Command finished successfully.")
        if process.stderr:
             # RAxML v8 often prints informational messages also to stderr, 
             # so check return code first. Print stderr if error or for info.
             print("Stderr output:\n", process.stderr) 
        return process
    except FileNotFoundError as e:
        print(f"Error: Command '{command_list[0]}' not found - {e}.")
        print(f"Please verify the 'raxml_executable' variable is correct ('{raxml_executable}')")
        print("and that standard RAxML is installed in your Conda environment.")
        sys.exit(1)
    except subprocess.CalledProcessError as e:
        print(f"Error executing command: {' '.join(command_list)}")
        print(f"Return code: {e.returncode}")
        # Output is not captured by default here, check stderr
        print(f"Stderr:\n{e.stderr}") 
        sys.exit(1)
    except Exception as e:
        print(f"An unexpected error occurred: {e}")
        sys.exit(1)

# --- Input Validation ---
# Check if RAxML executable exists (basic check)
if shutil.which(raxml_executable) is None:
    print(f"Error: RAxML executable '{raxml_executable}' not found in PATH.")
    print("Please ensure standard RAxML is installed and the 'raxml_executable' variable is correct.")
    sys.exit(1)
else:
    print(f"Found RAxML executable: {shutil.which(raxml_executable)}")

# Check if MAFFT directory exists
if not os.path.isdir(mafft_dir):
    print(f"Error: MAFFT alignment directory not found: {mafft_dir}")
    sys.exit(1)

# --- Main Workflow ---
# 1. Create output directory
print(f"\nCreating output directory: {raxml_v8_output_dir}")
os.makedirs(raxml_v8_output_dir, exist_ok=True)

# 2. Find alignment files
alignment_files = glob.glob(alignment_pattern)
if not alignment_files:
    print(f"Error: No alignment files found matching '{alignment_pattern}' in {mafft_dir}")
    sys.exit(1)

print(f"Found {len(alignment_files)} alignment files to process using {raxml_executable}.")

# 3. Run RAxML v8 for each alignment
for aln_file in alignment_files:
    base_name = os.path.basename(aln_file).replace(".aln.fa", "")
    print(f"\nProcessing alignment: {base_name} (File: {aln_file})")

    # Define the output name/label (-n option) for RAxML files
    output_label = base_name 
    
    # Define the working directory where RAxML will write files (-w option)
    # Using absolute path is safer
    working_dir = os.path.abspath(raxml_v8_output_dir)

    # Check if final output tree already exists to potentially skip
    # The key output from '-f a' is RAxML_bipartitions.<output_label>
    final_tree_file = os.path.join(working_dir, f"RAxML_bipartitions.{output_label}")
    if os.path.exists(final_tree_file):
        print(f"  Output file {final_tree_file} already exists. Skipping RAxML run for {base_name}.")
        continue
# --- Check number of sequences in alignment ---
    try:
        num_sequences = 0
        with open(aln_file, 'r') as handle:
            # Simple way: count '>' characters
            for line in handle:
                if line.startswith('>'):
                    num_sequences += 1
        
        # More robust way using Biopython (if you prefer)
        # from Bio import AlignIO
        # alignment = AlignIO.read(aln_file, "fasta")
        # num_sequences = len(alignment)

        if num_sequences < 4: # Or use 4 if you want slightly more meaningful trees
            print(f"  Skipping RAxML for {base_name}: Found only {num_sequences} sequences (requires > 3).")
            continue # Skip to the next alignment file in the loop
        else:
            print(f"  Found {num_sequences} sequences in {aln_file}. Proceeding with RAxML.")

    except Exception as e:
        print(f"  Warning: Could not count sequences in {aln_file}. Error: {e}. Skipping RAxML.")
        continue

    # Construct the RAxML v8 command using the -f a workflow
    raxml_cmd = [
        raxml_executable,
        "-f", "a",                         # Rapid bootstrap + ML search
        "-m", model,                      # Substitution model (e.g., PROTGAMMALGF)
        "-s", os.path.abspath(aln_file),  # Input alignment file (use absolute path)
        "-w", working_dir,                # FULL path to output directory
        "-n", output_label,               # Output file label
        "-T", str(num_threads),           # Number of threads
        "-x", str(random_seed_bootstrap), # Random seed for rapid bootstrap
        "-p", str(random_seed_parsimony), # Random seed for parsimony inferences
        "-#", str(bootstrap_reps)         # Number of bootstrap replicates (use '#', sometimes '-N')
    ]

    # Run the command
    # No need to specify cwd if we use absolute paths for -s and -w
    run_command(raxml_cmd) 
    print(f"  RAxML v8 analysis complete for {base_name}.")
    print(f"  Key output file (best tree with BS support): {final_tree_file}")
    print(f"  Other files (like RAxML_info.{output_label}, RAxML_bestTree.{output_label}) are in: {working_dir}")


print("\nRAxML v8 tree building process completed.")
print(f"Output files are located in: {os.path.abspath(raxml_v8_output_dir)}")
print("The primary result files for visualization are typically the ones ending in 'RAxML_bipartitions.<label>'.")
print("Remember: These trees are also UNROOTED. You will root them using FigTree.")
