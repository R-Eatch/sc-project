import os
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

def translate_and_rename_cds(cds_dir, output_dir):
    for filename in os.listdir(cds_dir):
        if filename.endswith(".cds"):
            organ_name = filename.split("_")[0]
            input_path = os.path.join(cds_dir, filename)
            output_path = os.path.join(output_dir, f"{organ_name}_proteins.fasta")
            
            # 翻译CDS序列
            protein_records = []
            for record in SeqIO.parse(input_path, "fasta"):
                protein_seq = record.seq.translate(to_stop=True)
                protein_record = SeqRecord(protein_seq, id=record.id, description="Translated Protein")
                protein_records.append(protein_record)
            
            # 保存翻译后的蛋白序列
            SeqIO.write(protein_records, output_path, "fasta")

# 指定CDS目录和输出目录路径
cds_dir = "./cds"
output_dir = "./proteins"

# 创建输出目录
os.makedirs(output_dir, exist_ok=True)

translate_and_rename_cds(cds_dir, output_dir)


