import re

# 读取GFF文件
gff_file = 'filter.id.gff'
output_file = 'sugerglider.gff'

# 用于存储GFF行的列表
gff_content = {}

# 读取GFF文件并解析每一行
with open(gff_file, 'r') as file:
    for line in file:
        if line.startswith('#'):
            continue  # 跳过注释行
        fields = line.strip().split('\t')
        if len(fields) >= 9 and (fields[2] == 'CDS'):
            match = re.match('Parent=([^;]+)', fields[8])
            # print(line.strip())
            gene_id = match.group(1)
            if gene_id not in gff_content:
                gff_content[gene_id] = {}
                gff_content[gene_id]['content'] = ''
                gff_content[gene_id]['chrom'] = fields[0]
                gff_content[gene_id]['start'] = float('inf')
                gff_content[gene_id]['end'] = float('-inf')
                gff_content[gene_id]['strand'] = fields[6]
            gff_content[gene_id]['start'] = min(gff_content[gene_id]['start'], int(fields[3]))
            gff_content[gene_id]['end'] = max(gff_content[gene_id]['end'], int(fields[4]))
            gff_content[gene_id]['content'] += line

items = list(gff_content.items())

sorted_items = sorted(items, key=lambda item: (item[1]['chrom'], int(item[1]['start'])))
new_id_prefix = 'sugarglider'
new_id_counter = 0
with open(output_file, 'w') as output:
    for gene_id, attributes in sorted_items:
        new_id_counter += 1
        new_id = f"{new_id_prefix}{new_id_counter:06d}"
        chrom = attributes['chrom']
        start = attributes['start']
        end = attributes['end']
        strand = attributes['strand']
        content = re.split('\n', attributes['content'].strip())
        mRNA_line = f"{chrom}\tminiprot\tmRNA\t{start}\t{end}\t.\t{strand}\t.\tID={new_id};Parent={new_id};"
        output.write(f'{mRNA_line}\n')
        for line in content:
            fields = line.strip().split('\t')
            fields[8] = f'ID={new_id};Parent={new_id};'
            new_line = '\t'.join(fields)
            output.write(f'{new_line}\n')
    
