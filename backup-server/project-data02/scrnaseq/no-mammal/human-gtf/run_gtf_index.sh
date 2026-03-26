#!/bin/bash

# --- Configuration ---
# Input GTF file provided by the user
INPUT_GTF="Homo_sapiens.GRCh38.113.gtf"

# Check if input file exists
if [ ! -f "${INPUT_GTF}" ]; then
  echo "Error: Input GTF file '${INPUT_GTF}' not found!"
  exit 1
fi

# Define output file names based on the input name
SORTED_GTF="${INPUT_GTF%.gtf}.sorted.gtf"
COMPRESSED_GTF="${SORTED_GTF}.gz"
INDEX_FILE="${COMPRESSED_GTF}.tbi"

echo "--- Starting GTF Processing ---"
echo "Input file: ${INPUT_GTF}"
echo "Sorted output (temporary): ${SORTED_GTF}"
echo "Compressed output: ${COMPRESSED_GTF}"
echo "Index output: ${INDEX_FILE}"
echo "---------------------------------"

# --- Step 1: Sort the GTF file ---
# Keeps header lines (#) at the top, sorts data lines by chromosome (col 1, version sort)
# and then by start coordinate (col 4, numeric sort).
echo "[Step 1/3] Sorting GTF file..."
(grep "^#" "${INPUT_GTF}"; grep -v "^#" "${INPUT_GTF}" | sort -k1,1V -k4,4n) > "${SORTED_GTF}"

# Basic check if sorting worked (file exists and is not empty)
if [ ! -s "${SORTED_GTF}" ]; then
    echo "Error: Sorting failed or produced an empty file."
    exit 1
fi
echo "Sorting complete."
echo "---------------------------------"

# --- Step 2: Compress the sorted file with bgzip ---
# bgzip is required for tabix indexing. Standard gzip will not work.
echo "[Step 2/3] Compressing sorted GTF with bgzip..."
bgzip "${SORTED_GTF}"

# Check if compression worked
if [ ! -f "${COMPRESSED_GTF}" ]; then
    echo "Error: bgzip compression failed."
    # Clean up the potentially large sorted file if compression failed
    rm -f "${SORTED_GTF}"
    exit 1
fi
# bgzip automatically removes the input file (SORTED_GTF) upon successful compression
echo "Compression complete. Created ${COMPRESSED_GTF}"
echo "---------------------------------"

# --- Step 3: Index the compressed file with tabix ---
echo "[Step 3/3] Indexing compressed GTF with tabix..."
tabix -p gff "${COMPRESSED_GTF}"

# Check if indexing worked
if [ ! -f "${INDEX_FILE}" ]; then
    echo "Error: tabix indexing failed."
    exit 1
fi
echo "Indexing complete. Created ${INDEX_FILE}"
echo "---------------------------------"

echo "--- Processing Finished Successfully! ---"
echo "You can now load the file:"
echo "  ${COMPRESSED_GTF}"
echo "into IGV. Ensure the index file:"
echo "  ${INDEX_FILE}"
echo "is in the same directory."
echo "---------------------------------"

exit 0
