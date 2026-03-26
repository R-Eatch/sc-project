def read_gff(file_name):
    """
    Reads a GFF file and stores all lines.
    """
    records = []
    with open(file_name, 'r') as file:
        for line in file:
            if not line.startswith('#'):
                records.append(line.strip())
    return records

def fix_records(records):
    """
    Fixes the format of the records, ensuring:
    1. Removal of extra semicolons in gene rows.
    2. Correction of Parent format in UTR rows and removal of the colour attribute.
    3. Removal of 'utr_' and '_1_1' in UTR row IDs, while preserving the semicolon in Parent.
    """
    fixed_records = []
    for record in records:
        parts = record.split('\t')
        if len(parts) < 9:
            continue  # Skip incomplete records
        feature_type, attributes = parts[2], parts[8]

        # Handle extra semicolon in gene rows
        if feature_type == "gene":
            attributes = attributes.replace(";;", ";")

        # Prepare for attribute modifications
        attrs = {attr.split('=')[0]: attr.split('=')[1] for attr in attributes.split(';') if '=' in attr}

        # Special handling for UTR rows
        if "UTR" in feature_type:
            if 'ID' in attrs:
                attrs['ID'] = attrs['ID'].replace('utr_', '').replace('_1_1', '')
            if 'Parent' in attrs:
                parent_id = attrs['Parent'].split('_')[0]
                attrs['Parent'] = parent_id
            attrs.pop('colour', None)  # Remove colour attribute

        # Reconstruct attributes string and ensure it ends with a semicolon
        attributes = ';'.join(f"{k}={v}" for k, v in attrs.items()) + ';'
        parts[8] = attributes
        fixed_record = '\t'.join(parts)
        fixed_records.append(fixed_record)
    return fixed_records

def write_fixed_gff(fixed_records, output_file_name):
    """
    Writes the fixed records to a new GFF file.
    """
    with open(output_file_name, 'w') as output_file:
        for record in fixed_records:
            output_file.write(record + '\n')

def main():
    input_gff_file = "C:\\Users\\eatch\\Desktop\\sugar_glider.gff"  # Update to the actual input file path
    output_gff_file = "C:\\Users\\eatch\\Desktop\\sugar_glider_new.gff"  # Update to the actual output file path

    records = read_gff(input_gff_file)
    fixed_records = fix_records(records)  # Apply format corrections
    write_fixed_gff(fixed_records, output_gff_file)
    print(f"Processed GFF file has been saved to {output_gff_file}")

if __name__ == "__main__":
    main()
