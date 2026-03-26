def split_fasta(input_file):
    with open(input_file, 'r') as f:
        lines = f.readlines()

    seq = ''
    seq_name = ''
    for line in lines:
        line = line.strip()
        if line.startswith('>'):
            if seq:
                with open(f"{seq_name}.fasta", 'w') as out:
                    out.write(f">{seq_name}\n{seq}\n")
            seq_name = line[1:]  # Remove the '>'
            seq = ''
        else:
            seq += line
    if seq:
        with open(f"{seq_name}.fasta", 'w') as out:
            out.write(f">{seq_name}\n{seq}\n")

# 使用示例
split_fasta('lalba_exon.fa')
