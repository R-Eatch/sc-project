import os
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

def translate_dna_to_protein(dna_seq):
    """将DNA序列翻译成蛋白质序列。这里使用简化的密码子表。"""
    codon_table = {
        'ATA':'I', 'ATC':'I', 'ATT':'I', 'ATG':'M',
        'ACA':'T', 'ACC':'T', 'ACG':'T', 'ACT':'T',
        'AAC':'N', 'AAT':'N', 'AAA':'K', 'AAG':'K',
        'AGC':'S', 'AGT':'S', 'AGA':'R', 'AGG':'R',
        'CTA':'L', 'CTC':'L', 'CTG':'L', 'CTT':'L',
        'CCA':'P', 'CCC':'P', 'CCG':'P', 'CCT':'P',
        'CAC':'H', 'CAT':'H', 'CAA':'Q', 'CAG':'Q',
        'CGA':'R', 'CGC':'R', 'CGG':'R', 'CGT':'R',
        'GTA':'V', 'GTC':'V', 'GTG':'V', 'GTT':'V',
        'GCA':'A', 'GCC':'A', 'GCG':'A', 'GCT':'A',
        'GAC':'D', 'GAT':'D', 'GAA':'E', 'GAG':'E',
        'GGA':'G', 'GGC':'G', 'GGG':'G', 'GGT':'G',
        'TCA':'S', 'TCC':'S', 'TCG':'S', 'TCT':'S',
        'TTC':'F', 'TTT':'F', 'TTA':'L', 'TTG':'L',
        'TAC':'Y', 'TAT':'Y', 'TAA':'*', 'TAG':'*',
        'TGC':'C', 'TGT':'C', 'TGA':'*', 'TGG':'W',
    }
    protein_seq = ''
    for i in range(0, len(dna_seq), 3):
        codon = dna_seq[i:i+3]
        protein_seq += codon_table.get(codon, 'X')  # 使用'X'标记未知的密码子
    return protein_seq

def translate_and_rename_cds(cds_dir, output_dir):
    for filename in os.listdir(cds_dir):
        if filename.endswith(".cds"):
            organ_name = filename.split("_")[0]
            input_path = os.path.join(cds_dir, filename)
            output_path = os.path.join(output_dir, f"{organ_name}_proteins.fasta")
            
            # 翻译CDS序列
            protein_records = []
            for record in SeqIO.parse(input_path, "fasta"):
                protein_seq_str = translate_dna_to_protein(str(record.seq))
                protein_record = SeqRecord(Seq(protein_seq_str), id=record.id, description="Translated Protein")
                protein_records.append(protein_record)
            
            # 保存翻译后的蛋白序列
            SeqIO.write(protein_records, output_path, "fasta")

# 指定CDS目录和输出目录路径
cds_dir = "./cds"
output_dir = "./proteins"

# 创建输出目录
os.makedirs(output_dir, exist_ok=True)

translate_and_rename_cds(cds_dir, output_dir)

