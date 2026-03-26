#!/usr/bin/env python3

import os
import subprocess
import glob
import csv
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
import sys
import shutil
import re # Import regular expressions module
from collections import defaultdict

# --- Configuration ---
# --- Configuration ---
# *** 选择模式 ***
# Options: "BEST_PER_SPECIES", "TOP_N_PER_SPECIES", "TOP_N_UNIQUE_IDS"
SELECTION_MODE = "TOP_N_PER_SPECIES"

# *** 模式对应的限制参数 ***
# 仅在 SELECTION_MODE = "TOP_N_UNIQUE_IDS" 时使用
TOP_N_UNIQUE_ID_LIMIT = 10
# 仅在 SELECTION_MODE = "TOP_N_PER_SPECIES" 时使用
TOP_HITS_PER_SPECIES_LIMIT = 3 # <--- 新增: 每个物种最多选几个 Hits

# *** 路径配置 (请确保相对于此脚本的运行位置正确) ***
blast_dir = "../2.blastp/"
query_dir = "../1.extra_pep/"
mafft_dir = "./"  # MAFFT 输出目录 (当前目录)

# 输入文件名 (基于上面的路径)
combined_db_fasta = os.path.join(blast_dir, "combined_species.pep.fa")
query_fasta = os.path.join(query_dir, "longest_queries.fa")

# BLAST 结果文件模式
blast_result_pattern = os.path.join(blast_dir, "*.blastp.tsv")

# *** 这个值只在 SELECT_BEST_PER_SPECIES = False 时生效 ***
top_n_hits = 10 # 作为 Top N 唯一 ID 的上限

# MAFFT 参数
mafft_threads = 20
mafft_options = ["--auto"] # MAFFT 比对选项

# BLAST TSV 列索引 (0-based)
QSEQID_COL = 0
SSEQID_COL = 1
BITSCORE_COL = 11

# *** 物种与前缀映射字典 (请确保准确无误) ***
species_to_prefix_map = {
    'Chicken': 'ENSGALP',
    'Opossum': 'ENSMODP',
    'Human': 'ENSP',
    'Mouse': 'ENSMUSP',
    'Lizard': 'ENSPMRP', # 确保这个前缀是正确的，并且比 'ENSP' 长
    'Frog': 'ENSXETP'
}
# 反向映射: 前缀 -> 物种名
prefix_to_species_map = {v: k for k, v in species_to_prefix_map.items()}

# *** 准备正则表达式用于物种识别 (处理前缀歧义) ***
try:
    # 1. 获取前缀并按长度降序排序
    sorted_prefixes = sorted(prefix_to_species_map.keys(), key=len, reverse=True)
    if not sorted_prefixes:
        raise ValueError("Prefix map is empty, cannot create regex.")
    # 2. 创建正则表达式字符串: ^(长前缀|...|短前缀)
    regex_pattern_str = r"^(" + "|".join(re.escape(p) for p in sorted_prefixes) + r")"
    # 3. 编译正则表达式以提高效率
    species_regex = re.compile(regex_pattern_str)
    print(f"Compiled species regex pattern: {species_regex.pattern}")
except Exception as e:
    print(f"Error preparing species regex: {e}")
    sys.exit(1)
# ---------------------

# --- 辅助函数: 使用正则表达式从 ID 获取物种名 ---
def get_species_from_id_regex(sequence_id, compiled_regex, prefix_map):
    """使用预编译的正则表达式确定物种名 (优先匹配长前缀)。"""
    match = compiled_regex.match(sequence_id)
    if match:
        matched_prefix = match.group(1) # 获取实际匹配的前缀
        return prefix_map.get(matched_prefix) # 使用匹配的前缀查找物种名
    return None # 没有已知前缀匹配

# --- 辅助函数: 运行外部命令 ---
def run_command(command_list, check=True, capture_output=False, text=True, stdout_file=None):
    """运行 shell 命令并处理基本错误。"""
    print(f"Running command: {' '.join(command_list)}")
    try:
        process = subprocess.run(command_list, check=check,
                                 capture_output=capture_output, text=text,
                                 stdout=stdout_file, stderr=subprocess.PIPE)
        # 打印标准错误（可能包含 MAFFT 的进度信息或警告）
        if process.stderr:
            print("Stderr:", process.stderr, file=sys.stderr)
        # 如果要求捕获输出，则打印标准输出
        if capture_output and process.stdout:
            print("Stdout:", process.stdout)
        return process
    except FileNotFoundError as e:
        print(f"Error: Command '{command_list[0]}' not found - {e}. Is the tool (e.g., MAFFT) installed and in PATH?")
        sys.exit(1) # 无法找到命令是致命错误
    except subprocess.CalledProcessError as e:
        # 命令执行失败 (非零退出码)
        print(f"Error executing command: {' '.join(command_list)}")
        print(f"Return code: {e.returncode}")
        if capture_output: print(f"Output (stdout):\n{e.stdout}")
        print(f"Output (stderr):\n{e.stderr}") # 打印详细错误信息
        # 根据需要决定是否退出，这里我们假设子进程错误应该被捕获并继续（例如 MAFFT 失败不停止整个脚本）
        # 但要抛出异常让调用者知道出错了
        raise e
    except Exception as e:
        print(f"An unexpected error occurred running command: {e}")
        raise e # 重新抛出异常

# --- 辅助函数: 使用正则表达式修改 FASTA 头部 ---
def modify_fasta_headers_regex(input_fasta, output_fasta, compiled_regex, prefix_map):
    """读取 FASTA, 使用 regex 修改头部 (优先长前缀), 写入新文件。"""
    print(f"  Modifying headers for: {os.path.basename(input_fasta)} using Regex")
    modified_records = []
    records_processed = 0
    headers_modified_count = 0

    try:
        for record in SeqIO.parse(input_fasta, "fasta"):
            records_processed += 1
            # 如果头部已包含'|', 取第一部分作为基础ID；否则用整个ID
            base_id = record.id.split('|')[0]
            species_name = None

            # 使用编译好的正则表达式匹配基础ID
            match = compiled_regex.match(base_id)
            if match:
                matched_prefix = match.group(1)
                species_name = prefix_map.get(matched_prefix)

            # 如果识别出物种名，并且原始ID中不包含该物种标签
            if species_name and f"|{species_name}" not in record.id:
                new_id = f"{base_id}|{species_name}" # 构建新ID

                # 仔细重构描述部分 (description)
                description_text = record.description
                if record.description.startswith(base_id): # 通常 description 以 ID 开头
                    desc_start_index = len(base_id)
                    # 跳过 ID 后可能的空格
                    if desc_start_index < len(record.description) and record.description[desc_start_index] == ' ':
                         desc_start_index += 1
                    description_text = record.description[desc_start_index:]

                # 创建带有新ID和保留描述的新记录
                new_record = SeqRecord(record.seq, id=new_id, description=description_text, name='') # name 通常为空
                modified_records.append(new_record)
                headers_modified_count += 1
            else:
                # 如果未识别出物种，或已包含正确的标签，则保留原始记录
                modified_records.append(record)

        # 将所有记录 (修改过的和未修改的) 写入输出文件
        with open(output_fasta, "w") as outfile:
            SeqIO.write(modified_records, outfile, "fasta")

        print(f"  Finished modifying headers. {headers_modified_count}/{records_processed} headers updated/verified.")
        return True # 表示成功

    except Exception as e:
        print(f"Error modifying headers for {input_fasta}: {e}")
        return False # 表示失败

# --- 输入验证 ---
print("--- Input Validation ---")
if shutil.which("mafft") is None:
    print("Error: 'mafft' command not found. Please ensure MAFFT is installed and in your PATH.")
    sys.exit(1)
if not os.path.isdir(blast_dir):
    print(f"Error: BLAST directory not found: {blast_dir}")
    sys.exit(1)
if not os.path.isdir(query_dir):
     print(f"Error: Query directory not found: {query_dir}")
     sys.exit(1)
if not os.path.exists(combined_db_fasta):
    print(f"Error: Combined database FASTA not found: {combined_db_fasta}")
    sys.exit(1)
if not os.path.exists(query_fasta):
    print(f"Error: Query FASTA file not found: {query_fasta}")
    sys.exit(1)
print("Input paths and MAFFT seem OK.")
print("--- End Validation ---")

# --- 主要工作流程 ---
print(f"\nEnsuring output directory exists: {os.path.abspath(mafft_dir)}")
os.makedirs(mafft_dir, exist_ok=True)

# 索引数据库和查询文件
print(f"\nIndexing database FASTA: {combined_db_fasta} ...")
try:
    db_index = SeqIO.index(combined_db_fasta, "fasta")
    print(f"Indexed {len(db_index)} sequences.")
except Exception as e:
    print(f"Error indexing database FASTA: {e}")
    sys.exit(1)

print(f"Indexing query FASTA: {query_fasta} ...")
try:
    query_index = SeqIO.index(query_fasta, "fasta")
    print(f"Indexed {len(query_index)} query sequences.")
except Exception as e:
    print(f"Error indexing query FASTA: {e}")
    sys.exit(1)

# 查找 BLAST 结果文件
blast_result_files = glob.glob(blast_result_pattern)
if not blast_result_files:
    print(f"Error: No BLAST result files found matching pattern: '{blast_result_pattern}'")
    sys.exit(1)

print(f"\nFound {len(blast_result_files)} BLAST result files to process.")
# 告知用户当前的筛选模式
#print(f"Sequence selection mode: {'Best Hit Per Species' if SELECT_BEST_PER_SPECIES else f'Top {top_n_hits} Unique Subject IDs'}")

# --- 循环处理每个 BLAST 文件 ---
for blast_file in blast_result_files:
    # 从文件名获取基因基本名
    base_name = os.path.basename(blast_file).split('.')[0]
    print(f"\n>>> Processing Gene: {base_name} (File: {blast_file})")

    # --- 读取并解析该基因的所有 BLAST Hits ---
    hits = []
    query_id_for_this_file = None
    line_num = 0
    print(f"  Reading BLAST hits...")
    try:
        with open(blast_file, 'r', newline='') as tsvfile:
            reader = csv.reader(tsvfile, delimiter='\t')
            for row in reader:
                line_num += 1
                if not row or row[0].startswith("#"): continue # 跳过空行/注释行
                try:
                    if len(row) > max(SSEQID_COL, BITSCORE_COL):
                        qseqid = row[QSEQID_COL]
                        sseqid = row[SSEQID_COL]
                        bitscore = float(row[BITSCORE_COL])
                        hits.append({'qseqid': qseqid, 'sseqid': sseqid, 'bitscore': bitscore, 'line': line_num})
                        if query_id_for_this_file is None: query_id_for_this_file = qseqid
                    else: print(f"    Warning: Skipping malformed row {line_num}: {row}", file=sys.stderr)
                except (ValueError, IndexError) as e: print(f"    Warning: Error parsing row {line_num}: {row} - {e}", file=sys.stderr); continue
    except FileNotFoundError: print(f"  Error: File not found during read: {blast_file}"); continue # 跳到下一个基因
    except Exception as e: print(f"  Error reading BLAST file {blast_file}: {e}"); continue # 跳到下一个基因

    if not hits: print(f"  No valid hits found in file. Skipping this gene."); continue
    if query_id_for_this_file is None: print(f"  Could not determine query ID from hits. Skipping this gene."); continue
    print(f"  Read {len(hits)} valid hits.")

    # --- 按 Bit Score 降序排序所有 Hits ---
    hits.sort(key=lambda x: x['bitscore'], reverse=True)

    # --- 应用筛选逻辑 (重构) ---
    final_subject_ids = []
    print(f"  Applying sequence selection logic ({SELECTION_MODE})...")

    if SELECTION_MODE == "TOP_N_PER_SPECIES":
        # 模式: 每个物种选 Top N (由 TOP_HITS_PER_SPECIES_LIMIT 定义)
        species_hit_counts = defaultdict(int) # 使用 defaultdict 简化计数
        for hit in hits:
            sseqid = hit['sseqid']
            species_name = get_species_from_id_regex(sseqid, species_regex, prefix_to_species_map)

            if species_name is None: continue # 跳过未知物种

            # 检查该物种已选数量是否小于限制
            if species_hit_counts[species_name] < TOP_HITS_PER_SPECIES_LIMIT:
                final_subject_ids.append(sseqid)
                species_hit_counts[species_name] += 1 # 增加计数
                # print(f"    Selected {sseqid} for {species_name} (Count: {species_hit_counts[species_name]})") # 可选调试

    elif SELECTION_MODE == "BEST_PER_SPECIES":
        # 模式: 每个物种选最佳 (Top 1)
        species_already_selected = set()
        for hit in hits:
            sseqid = hit['sseqid']
            species_name = get_species_from_id_regex(sseqid, species_regex, prefix_to_species_map)
            if species_name is None: continue
            if species_name not in species_already_selected:
                final_subject_ids.append(sseqid)
                species_already_selected.add(species_name)

    elif SELECTION_MODE == "TOP_N_UNIQUE_IDS":
        # 模式: Top N 个唯一 Subject ID (由 TOP_N_UNIQUE_ID_LIMIT 定义)
        seen_subjects = set()
        for hit in hits:
            sseqid = hit['sseqid']
            if sseqid not in seen_subjects:
                final_subject_ids.append(sseqid)
                seen_subjects.add(sseqid)
                if len(final_subject_ids) >= TOP_N_UNIQUE_ID_LIMIT:
                    print(f"      Reached limit of {TOP_N_UNIQUE_ID_LIMIT} unique IDs.")
                    break
    # --- 检查是否有选出的 ID ---
    if not final_subject_ids:
        print(f"  No subject IDs selected after applying filter. Skipping alignment for this gene.")
        continue
    print(f"  Selected {len(final_subject_ids)} subject IDs for alignment.")

    # --- 检索序列 (查询序列 + 选出的 Subject 序列) ---
    print(f"  Retrieving sequences...")
    sequences_to_align = []
    # 首先获取查询序列
    try:
        query_record = query_index[query_id_for_this_file]
        sequences_to_align.append(query_record)
        print(f"    Retrieved query: {query_id_for_this_file}")
    except KeyError:
        print(f"  Error: Query ID '{query_id_for_this_file}' not found in query index '{query_fasta}'. Skipping this gene.")
        continue

    # 获取选出的 Subject 序列
    retrieved_subject_count = 0
    missing_subjects = []
    for sseqid in final_subject_ids:
        try:
            subject_record = db_index[sseqid]
            sequences_to_align.append(subject_record)
            retrieved_subject_count += 1
        except KeyError:
            print(f"    Warning: Subject ID '{sseqid}' not found in the database index '{combined_db_fasta}'. It will be excluded.")
            missing_subjects.append(sseqid)
    print(f"    Retrieved {retrieved_subject_count} subject sequences.")

    # 检查是否有足够序列进行比对 (至少需要2条)
    if len(sequences_to_align) < 2:
         print(f"  Only {len(sequences_to_align)} sequence(s) retrieved in total (need >= 2). Skipping MAFFT for this gene.")
         continue

    # --- 准备临时文件和输出文件名 ---
    temp_fasta_file = os.path.join(mafft_dir, f"{base_name}_temp_input.fa")
    alignment_output_file = os.path.join(mafft_dir, f"{base_name}.aln.fa")
    temp_modified_header_file = alignment_output_file + ".mod_tmp" # 用于头部修改的临时文件

    # 定义将在 finally 中使用的变量，即使 try 失败
    mafft_success = False
    header_mod_success = False

    # --- 主要处理块 (写文件, 运行MAFFT, 修改头部) 带清理 ---
    try:
        # 1. 写入 MAFFT 的临时输入文件
        print(f"  Preparing temporary input FASTA: {os.path.basename(temp_fasta_file)}")
        try:
            with open(temp_fasta_file, "w") as temp_out:
                SeqIO.write(sequences_to_align, temp_out, "fasta")
        except IOError as e:
            print(f"  Error writing temporary FASTA file: {e}")
            raise # 抛出异常，跳到外层 except，然后执行 finally

        # 2. 运行 MAFFT
        print(f"  Running MAFFT (Output: {os.path.basename(alignment_output_file)})...")
        mafft_command = ["mafft", f"--thread", str(mafft_threads)] + mafft_options + [temp_fasta_file]
        try:
            with open(alignment_output_file, "w") as outfile:
                # 使用 run_command，它会处理子进程错误
                run_command(mafft_command, stdout_file=outfile, check=True) # check=True 确保非零退出码会抛异常
            print(f"  MAFFT alignment raw output saved.")
            mafft_success = True # 标记 MAFFT 成功
        except subprocess.CalledProcessError as e:
            # MAFFT 命令执行失败 (run_command 中已打印错误信息)
            print(f"  MAFFT execution failed.")
            # 清理可能产生的空或不完整的输出文件
            if os.path.exists(alignment_output_file):
                try: os.remove(alignment_output_file)
                except OSError: pass
            # 不再继续尝试修改头部，mafft_success 保持 False
        except Exception as e: # 捕获 run_command 可能抛出的其他错误
             print(f"  An unexpected error occurred during MAFFT execution: {e}")
             if os.path.exists(alignment_output_file):
                 try: os.remove(alignment_output_file)
                 except OSError: pass

        # 3. 修改头部 (仅当 MAFFT 成功时)
        if mafft_success and os.path.exists(alignment_output_file):
            print(f"  Attempting header modification...")
            try:
                # 使用修正后的 regex 函数进行修改
                if modify_fasta_headers_regex(alignment_output_file, temp_modified_header_file, species_regex, prefix_to_species_map):
                    # 替换原始文件
                    try:
                        os.replace(temp_modified_header_file, alignment_output_file)
                        print(f"  Successfully replaced alignment file with modified headers.")
                        header_mod_success = True # 标记头部修改成功
                    except OSError as e_rep:
                         print(f"  Error replacing alignment file with modified version: {e_rep}")
                         # 如果替换失败，清理临时的修改文件
                         if os.path.exists(temp_modified_header_file):
                              try: os.remove(temp_modified_header_file)
                              except OSError: pass
                else:
                    # 修改函数返回 False (内部已打印错误)
                    print(f"  Header modification function failed. Keeping original MAFFT output.")
                    # 清理失败的临时修改文件
                    if os.path.exists(temp_modified_header_file):
                         try: os.remove(temp_modified_header_file)
                         except OSError: pass
            except Exception as e_mod: # 捕获修改/替换过程中的其他意外错误
                 print(f"  An unexpected error occurred during header modification/replacement: {e_mod}")
                 if os.path.exists(temp_modified_header_file):
                      try: os.remove(temp_modified_header_file)
                      except OSError: pass

    # --- 外层错误处理 ---
    except IOError: # 处理写入 temp_fasta_file 的错误
        print(f"  Skipping further processing for {base_name} due to error writing input file.")
    except Exception as e_outer: # 处理这个基因主 try 块中的其他意外错误
        print(f"  An unexpected error occurred processing {base_name} before cleanup: {e_outer}")

    # --- 清理 ---
    finally:
        # 这个块总会执行，无论 try 中是否发生错误
        print(f"  Cleaning up temporary files for {base_name}...")
        # 清理 MAFFT 的临时输入文件
        if os.path.exists(temp_fasta_file):
            try:
                os.remove(temp_fasta_file)
                # print(f"    - Removed {os.path.basename(temp_fasta_file)}") # 可选
            except OSError as e:
                print(f"    Warning: Could not remove temporary input file {temp_fasta_file}: {e}")

        # 清理头部修改的临时文件，仅当修改/替换未成功时
        if not header_mod_success and os.path.exists(temp_modified_header_file):
            try:
                os.remove(temp_modified_header_file)
                # print(f"    - Removed intermediate header mod file {os.path.basename(temp_modified_header_file)}") # 可选
            except OSError as e:
                print(f"    Warning: Could not remove temp modification file {temp_modified_header_file}: {e}")
        print(f"  Cleanup finished for {base_name}.")

# --- 脚本结束 ---
print(f"\n--- Workflow finished ---")
print(f"Final alignment files (with modified headers where applicable) are located in: {os.path.abspath(mafft_dir)}")
