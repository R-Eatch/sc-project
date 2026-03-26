from Bio import SeqIO
import os

def merge_and_rename_fasta_files(input_dir, output_file):
    """
    合并指定目录下的所有FASTA文件到一个单独的FASTA文件，并重命名序列ID。
    参数:
    input_dir : 包含FASTA文件的目录路径。
    output_file : 输出FASTA文件的路径。
    """
    fasta_files = [os.path.join(input_dir, f) for f in os.listdir(input_dir) if f.endswith('.fasta') or f.endswith('.fa')]
    all_sequences = []
    sequence_counter = 1  # 序列计数器，用于生成序列号

    for file in fasta_files:
        for record in SeqIO.parse(file, "fasta"):
            new_id = f"0.p{sequence_counter}"  # 生成新的序列ID
            record.id = new_id
            record.description = ""  # 清除原有的描述信息，只保留新的ID
            all_sequences.append(record)
            sequence_counter += 1

    SeqIO.write(all_sequences, output_file, "fasta")

input_dir = '../proteins'  # 更改为包含你的FASTA文件的目录
output_file = './new_merged.fasta'  # 指定你的输出文件路径

merge_and_rename_fasta_files(input_dir, output_file)

