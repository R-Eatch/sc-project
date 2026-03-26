from Bio import SeqIO
from Bio.Seq import Seq

def find_longest_orf(seq):
    longest_orf = ""
    for strand, nuc in [(+1, seq), (-1, seq.reverse_complement())]:
        for frame in range(3):
            length = 3 * ((len(seq) - frame) // 3)  # ORF长度必须是3的倍数
            for pro in nuc[frame:frame+length].translate(to_stop=True).split("*"):
                if len(pro) > len(longest_orf):
                    longest_orf = pro
    return longest_orf

input_file = 'TransLiG.fa'  # 请替换为您的文件路径
output_file = 'protein_midaiwu1.fa'

with open(output_file, 'w') as output_handle:
    for record in SeqIO.parse(input_file, "fasta"):
        longest_orf = find_longest_orf(record.seq)
        if longest_orf:  # 确保有ORF被找到
            record.seq = Seq(str(longest_orf))
            SeqIO.write(record, output_handle, "fasta")
