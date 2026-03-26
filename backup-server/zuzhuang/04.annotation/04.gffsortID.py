def read_and_group_records_by_id(file_name):
    """
    Reads a GFF file and groups records by their gene ID.
    """
    records_by_id = {}
    with open(file_name, 'r') as file:
        for line in file:
            if line.startswith('#') or not line.strip():
                continue
            parts = line.strip().split('\t')
            attributes = parts[8]

            # Extract ID from attributes
            id_attr = attributes.split(';')[0]
            id_value = id_attr.split('=')[1]

            if id_value not in records_by_id:
                records_by_id[id_value] = []

            records_by_id[id_value].append(line.strip())

    return records_by_id

def sort_records_by_id_numerically(records_by_id):
    """
    Sorts records by their ID numerically and combines them.
    """
    sorted_combined_records = []

    # Extract numeric part of IDs and sort
    sorted_ids = sorted(records_by_id.keys(), key=lambda x: int(x.replace('sugarglider', '')))

    for id_value in sorted_ids:
        sorted_combined_records.extend(records_by_id[id_value])

    return sorted_combined_records

def write_sorted_gff(sorted_records, output_file_name):
    """
    Writes the numerically sorted records to a new GFF file.
    """
    with open(output_file_name, 'w') as output_file:
        for record in sorted_records:
            output_file.write(record + '\n')

def main():
    input_gff_file = "C:\\Users\\eatch\\Desktop\\sugar_glider_new2.gff"  # Update to the actual input file path, which has been fixed by the first script
    output_gff_file = "C:\\Users\\eatch\\Desktop\\sugar_glider_new3.gff"  # Update to the actual output file path

    records_by_id = read_and_group_records_by_id(input_gff_file)
    sorted_combined_records = sort_records_by_id_numerically(records_by_id)
    write_sorted_gff(sorted_combined_records, output_gff_file)
    print(f"Numerically sorted GFF file has been saved to {output_gff_file}")

if __name__ == "__main__":
    main()
