# -*- coding: utf-8 -*-
import re
import csv
import os
import itertools # Needed for zip_longest
from collections import defaultdict

# --- Configuration ---
INPUT_FILENAME = 'orthogroup_summary.txt' # 请替换为你的实际输入文件名
RAW_OUTPUT_FILENAME = 'filtered_raw_data.csv' # 第一个CSV：原始筛选行
FINAL_GENES_OUTPUT_FILENAME = 'mammalian_new_genes_final.csv' # 最终CSV：每个物种所有独立基因
STATS_OUTPUT_FILENAME = 'processing_stats.txt' # 统计结果TXT文件
TARGET_SPECIES = {'human', 'mouse', 'rat', 'pig'} # 目标物种集合

# --- Core Logic Function ---
def process_gene_data(input_file, raw_output_file, final_genes_output_file, stats_output_file, target_species_set):
    """
    Filters gene data lines, extracts all individual genes per target species,
    saves raw filtered lines and final split gene data, and writes statistics to a file.

    Args:
        input_file (str): Path to the input text file.
        raw_output_file (str): Path for the CSV file with raw filtered lines.
        final_genes_output_file (str): Path for the final CSV with all split genes per species.
        stats_output_file (str): Path for the text file to store statistics.
        target_species_set (set): A set containing the exact species names to filter for.
    """
    filtered_raw_lines = [] # List for CSV 1 (Raw lines)
    # Use defaultdict(list) to collect *all* individual genes for each target species across all lines
    all_species_genes = defaultdict(list)

    # Statistics Counters
    filtered_line_count = 0
    genes_eq_4_count = 0
    genes_gt_4_count = 0
    lines_missing_gene_info = 0
    total_genes_extracted_per_species = defaultdict(int) # Counts *after* potential splitting if needed (though splitting is handled implicitly by appending)

    # --- Regular Expressions ---
    all_species_pattern = re.compile(r"([a-zA-Z]+)\|([^,() ;]+)") # Find *ANY* species|gene
    target_gene_pattern = re.compile(r"(human|mouse|rat|pig)\|([^,() ;]+)") # Find *TARGET* species|gene
    info_pattern = re.compile(r"Genes:\s*(\d+);\s*Species:\s*(\d+)") # Extract "Genes: X; Species: Y;"

    print(f"开始处理文件: {input_file}")

    try:
        with open(input_file, 'r', encoding='utf-8') as infile:
            for line_num, line in enumerate(infile, 1):
                line = line.strip()
                if not line:
                    continue

                # --- Strict Filtering Logic ---
                all_matches_in_line = all_species_pattern.findall(line)
                all_found_species_in_line = {match[0] for match in all_matches_in_line}

                if all_found_species_in_line == target_species_set:
                    # Line passed the strict species filter
                    filtered_line_count += 1

                    # --- Store Raw Filtered Line (CSV 1) ---
                    filtered_raw_lines.append([line])

                    # --- Parse "Genes: X; Species: Y;" Info (For Stats) ---
                    info_match = info_pattern.search(line)
                    line_gene_count = -1
                    if info_match:
                        try:
                            line_gene_count = int(info_match.group(1))
                            if line_gene_count == 4:
                                genes_eq_4_count += 1
                            elif line_gene_count > 4:
                                genes_gt_4_count += 1
                        except (ValueError, IndexError):
                            # Silently ignore parse errors here, but count them below
                            lines_missing_gene_info += 1
                    else:
                        lines_missing_gene_info += 1

                    # --- Extract All Genes for Target Species (For Final CSV) ---
                    # This pattern finds the specific genes we want to put in the final CSV
                    target_matches_in_line = target_gene_pattern.findall(line)
                    genes_found_in_this_line = 0
                    for species, gene in target_matches_in_line:
                        # Directly append each found gene to the global list for that species
                        # No explicit splitting needed here, as findall already gives individual pairs
                        # e.g., (rat|GeneA1,rat|GeneA2) yields ('rat', 'GeneA1') and ('rat', 'GeneA2')
                        all_species_genes[species].append(gene)
                        total_genes_extracted_per_species[species] += 1 # Count every extracted gene
                        genes_found_in_this_line += 1

                    # Optional: Add a sanity check if the number of target genes found differs from Genes: count
                    if line_gene_count != -1 and genes_found_in_this_line != line_gene_count:
                         print(f"注意 (Line {line_num}): 'Genes:' 计数 ({line_gene_count}) 与实际提取的目标物种基因数 ({genes_found_in_this_line}) 不符。行: {line}")


    except FileNotFoundError:
        print(f"错误: 输入文件未找到于 {input_file}")
        return
    except Exception as e:
        print(f"在文件读取或处理过程中发生意外错误: {e}")
        return

    # --- Write Statistics to TXT file (Req 1) ---
    print(f"\n筛选完成。总共找到 {filtered_line_count} 行满足物种条件。")
    try:
        with open(stats_output_file, 'w', encoding='utf-8') as statsfile:
            statsfile.write("基因数据处理统计结果\n")
            statsfile.write("=" * 30 + "\n")
            statsfile.write(f"处理的输入文件: {input_file}\n")
            statsfile.write(f"筛选条件: 物种必须精确为 {', '.join(sorted(list(target_species_set)))}\n")
            statsfile.write("-" * 30 + "\n")
            statsfile.write(f"总共筛选出的行数: {filtered_line_count}\n")
            statsfile.write(f"  - 其中 'Genes: 4' 的行数: {genes_eq_4_count}\n")
            statsfile.write(f"  - 其中 'Genes: > 4' 的行数: {genes_gt_4_count}\n")
            unknown_gene_count = filtered_line_count - genes_eq_4_count - genes_gt_4_count
            if unknown_gene_count > 0:
                 statsfile.write(f"  - 其中 'Genes:' 数量未知或 < 4 的行数: {unknown_gene_count}\n")
            if lines_missing_gene_info > 0 and unknown_gene_count == lines_missing_gene_info :
                 statsfile.write(f"    (注: 上述 {lines_missing_gene_info} 行未能成功解析 'Genes: X; Species: Y;' 信息)\n")
            elif lines_missing_gene_info > 0:
                 statsfile.write(f"    (注: 有 {lines_missing_gene_info} 行未能成功解析 'Genes: X; Species: Y;' 信息)\n")

            statsfile.write("-" * 30 + "\n")
            statsfile.write(f"最终提取并写入 '{final_genes_output_file}' 的各物种基因总数:\n")
            header_for_stats = sorted(list(target_species_set))
            for species in header_for_stats:
                count = total_genes_extracted_per_species[species]
                statsfile.write(f"  - {species}: {count}\n")
        print(f"统计结果已保存到: {stats_output_file}")
    except IOError as e:
        print(f"错误: 写入统计文件 {stats_output_file} 时出错: {e}")
    except Exception as e:
        print(f"写入统计文件时发生意外错误: {e}")


    # --- Write Raw Filtered Data (CSV 1) ---
    if filtered_raw_lines:
        try:
            with open(raw_output_file, 'w', newline='', encoding='utf-8') as csvfile:
                writer = csv.writer(csvfile)
                writer.writerow(['Filtered_Line_Data'])
                writer.writerows(filtered_raw_lines)
            print(f"成功将原始筛选数据保存到: {raw_output_file}")
        except IOError as e:
            print(f"错误：写入原始筛选数据到 {raw_output_file} 时出错: {e}")
        except Exception as e:
            print(f"写入原始 CSV 时发生意外错误: {e}")
    else:
         print("没有原始数据被筛选出来。")


    # --- Write Final All Genes Data (CSV 4 - Req 2, 3, 4, 5) ---
    if any(all_species_genes.values()): # Check if any genes were collected
        header = sorted(list(target_species_set))
        # Prepare lists of genes in the correct order for zip_longest
        gene_lists_in_order = [all_species_genes[species] for species in header]

        try:
            with open(final_genes_output_file, 'w', newline='', encoding='utf-8') as csvfile:
                writer = csv.writer(csvfile)
                writer.writerow(header) # Write header row
                # Use zip_longest to write columns of potentially different lengths
                # '*' unpacks the list of lists into individual arguments for zip_longest
                # fillvalue='' ensures empty cells are written as empty strings
                writer.writerows(itertools.zip_longest(*gene_lists_in_order, fillvalue=''))
            print(f"成功将最终拆分的基因数据保存到: {final_genes_output_file}")
        except IOError as e:
            print(f"错误：写入最终基因数据到 {final_genes_output_file} 时出错: {e}")
        except Exception as e:
            print(f"写入最终基因 CSV 时发生意外错误: {e}")
    else:
        print(f"没有收集到任何目标物种的基因用于写入最终文件 ({final_genes_output_file})。")

    print("\n脚本执行完毕。")


# --- Script Execution ---
if __name__ == "__main__":
    # Create a dummy input file for testing if it doesn't exist
    if not os.path.exists(INPUT_FILENAME):
        print(f"警告: 输入文件 '{INPUT_FILENAME}' 未找到。正在创建一个用于测试的虚拟文件。")
        dummy_data = """
Node: n61; Genes: 4; Species: 4;        ((pig|KRT25,human|KRT25)n62,(rat|Krt25,mouse|Krt25)n63);
Node: n105; Genes: 5; Species: 5;       ((((mouse|Krt19,rat|Krt19)n108,human|KRT19)n107,chicken|ENSGALT00000090167)n106,frog|Frog020414);
Node: n33; Genes: 4; Species: 4;        (pig|KRT32,((rat|Krt32,mouse|Krt32)n35,human|KRT32)n34);
Node: n46; Genes: 6; Species: 6;        (((((rat|Krt12,mouse|Krt12)n50,human|KRT12)n49,pig|KRT12)n48,chicken|ENSGALT00000006107)n47,frog|Frog020416);
Node: n99; Genes: 4; Species: 4;        (human|GENEX,rat|GeneX,mouse|GeneX,pig|GENEX);
Node: n100; Genes: 3; Species: 3;       (human|GENEY,rat|GeneY,mouse|GeneY);
Node: n118; Genes: 8; Species: 8;       (frog|Frog031698,(lizard|IchBan12739,((human|OR6B1,(mouse|Olfr449,(rat|Olr811,rat|LOC103692138)n123)n122)n121,pig|ENSSSCT00000044414)n120)n119);
Node: n130; Genes: 4; Species: 4;       (human|GENEZ,rat|GeneZ,mouse|GeneZ,pig|GENEZ);
Node: n140; Genes: 5; Species: 4;       (human|GENEA,(rat|GeneA1,rat|GeneA2)n142,(mouse|GeneA,pig|GENEA)n143); # Valid species, Genes > 4, multiple rat genes
Node: n150; Genes: 4; Species: 4;       (human|GENEB,rat|GeneB,mouse|GeneB,pig|GENEB);
Node: n160; Species: 4; Genes: 4;        (human|GENEC,rat|GeneC,mouse|GeneC,pig|GENEC);
Node: n170; Info Missing;        (human|GENED,rat|GeneD,mouse|GeneD,pig|GENED);
Node: n180; Genes: 4; Species: 4;       (human|GENEE,rat|GeneE,mouse|GeneE,pig|GENEE); # Another standard case
        """
        try:
            with open(INPUT_FILENAME, 'w', encoding='utf-8') as f:
                f.write(dummy_data.strip())
        except IOError as e:
             print(f"错误：创建虚拟输入文件时出错: {e}")
             exit()

    # Run the main processing function with updated arguments
    process_gene_data(INPUT_FILENAME,
                      RAW_OUTPUT_FILENAME,
                      FINAL_GENES_OUTPUT_FILENAME, # Pass the final CSV name
                      STATS_OUTPUT_FILENAME,       # Pass the stats TXT name
                      TARGET_SPECIES)
