def read_and_group_records(file_name):
    """
    Reads a GFF file and groups records by their gene ID.
    """
    records_by_gene = {}
    with open(file_name, 'r') as file:
        for line in file:
            if line.startswith('#') or not line.strip():
                continue
            parts = line.strip().split('\t')
            feature_type = parts[2]
            attributes = parts[8]

            # Extract gene ID from attributes
            gene_id = attributes.split(';')[0].split('=')[1]

            if gene_id not in records_by_gene:
                records_by_gene[gene_id] = {'gene': [], 'mRNA': [], 'CDS': [], 'UTR': []}

            # Assign records to the correct category within their gene
            if feature_type in ['three_prime_UTR', 'five_prime_UTR']:
                records_by_gene[gene_id]['UTR'].append(line.strip())
            elif feature_type == 'CDS':
                records_by_gene[gene_id]['CDS'].append(
                    (int(parts[3]), line.strip()))  # Include start position for sorting
            elif feature_type == 'mRNA':
                records_by_gene[gene_id]['mRNA'].append(line.strip())
            elif feature_type == 'gene':
                records_by_gene[gene_id]['gene'].append(line.strip())

    return records_by_gene


def sort_and_combine_records(records_by_gene):
    """
    Sorts CDS records within each gene based on their start positions and combines all records.
    """
    sorted_combined_records = []

    for gene_id, record_types in records_by_gene.items():
        # Sort CDS records by their start position
        sorted_cds = sorted(record_types['CDS'], key=lambda x: x[0])

        # Combine all records, ensuring CDS are between mRNA and UTR
        combined_records = record_types['gene'] + record_types['mRNA'] + [cds_record[1] for cds_record in sorted_cds] + \
                           record_types['UTR']
        sorted_combined_records.extend(combined_records)

    return sorted_combined_records


def write_sorted_gff(sorted_records, output_file_name):
    """
    Writes the sorted records to a new GFF file.
    """
    with open(output_file_name, 'w') as output_file:
        for record in sorted_records:
            output_file.write(record + '\n')


def main():
    input_gff_file = "C:\\Users\\eatch\\Desktop\\sugar_glider_new.gff"  # Update to the actual input file path, which has been fixed by the first script
    output_gff_file = "C:\\Users\\eatch\\Desktop\\sugar_glider_new2.gff"  # Update to the actual output file path

    records_by_gene = read_and_group_records(input_gff_file)
    sorted_combined_records = sort_and_combine_records(records_by_gene)
    write_sorted_gff(sorted_combined_records, output_gff_file)
    print(f"Sorted GFF file has been saved to {output_gff_file}")


if __name__ == "__main__":
    main()


