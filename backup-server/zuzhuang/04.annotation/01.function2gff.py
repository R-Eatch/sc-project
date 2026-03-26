def read_annotation(file_name):
    """
    读取注释文件，并存储每个ID的注释信息及统计信息。
    """
    annotation = {}
    stats = {'GN': 0, 'KEGG': 0, 'EC': 0, 'Unknown': 0}
    with open(file_name, 'r') as file:
        for line in file:
            parts = line.strip().split('\t')
            gene_id = parts[0]
            gene_name = None
            prediction_text = " ".join(parts[1:])  # 将剩余部分作为预测文本
            for part in parts:
                if 'GN=' in part:
                    gene_name = part.split('GN=')[1].split()[0]
                    stats['GN'] += 1
                    break
                elif part.startswith('K'):
                    if gene_name is None:
                        gene_name = part
                        stats['KEGG'] += 1
                elif part.startswith('EC:'):
                    if gene_name is None:
                        gene_name = part
                        stats['EC'] += 1
            if gene_name is None:
                gene_name = '-'
                stats['Unknown'] += 1
            annotation[gene_id] = (gene_name, prediction_text)
    return annotation, stats


def integrate_gff(gff_file_name, annotation, output_file_name):
    """
    将注释信息整合到GFF文件中，并仅在gene行添加预测文本。
    """
    with open(gff_file_name, 'r') as gff_file, open(output_file_name, 'w') as output_file:
        for line in gff_file:
            if line.startswith('#') or line.strip() == '':
                output_file.write(line)  # 复制注释和空行
                continue
            parts = line.strip().split('\t')
            feature_type = parts[2]  # 获取特征类型（gene, mRNA, CDS等）
            attributes = parts[-1]
            gene_id = attributes.split(';')[0].split('=')[1]
            if gene_id in annotation:
                gene_name, prediction_text = annotation[gene_id]
                # 仅在gene行添加基因名称和预测文本
                if feature_type == 'gene':
                    attributes += f';Name={gene_name};prediction_text={prediction_text}'
            parts[-1] = attributes
            output_file.write('\t'.join(parts) + '\n')


def main():
    annotation_file = "C:\\Users\\eatch\\Desktop\\test.pep.all.function"  # 替换为你的注释文件路径
    gff_file = "C:\\Users\\eatch\\Desktop\\test.noStopCodon.new.gff3"  # 替换为你的GFF文件路径
    output_gff_file = "C:\\Users\\eatch\\Desktop\\output.gff"  # 输出文件的路径

    annotation, stats = read_annotation(annotation_file)
    integrate_gff(gff_file, annotation, output_gff_file)
    print(f"Integration complete. Output file: {output_gff_file}")
    print(
        f"Statistics: GN names: {stats['GN']}, KEGG IDs: {stats['KEGG']}, EC numbers: {stats['EC']}, Unknown: {stats['Unknown']}")


if __name__ == "__main__":
    main()


