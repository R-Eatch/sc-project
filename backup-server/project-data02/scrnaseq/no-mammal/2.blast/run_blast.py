#!/usr/bin/env python
"""
脚本功能：
1. 使用步骤1生成的查询蛋白 FASTA 文件作为查询序列，
   对包含非哺乳动物蛋白序列的 FASTA 文件构建 BLAST 数据库，
   并调用 blastp 进行比对。
2. 比对参数设定较严格（例如：E-value 1e-10、最小覆盖率80%、最小百分比身份50%）。
3. 脚本新增 --threads 参数，使 blastp 使用多线程运算（默认4个线程）。
4. BLAST 输出为 tabular 格式，默认文件名 blast_results.tsv。
注意：运行此脚本前需确保 BLAST+ 软件已安装（blastp、makeblastdb 可用）。
"""

import argparse
import subprocess
import os

def parse_args():
    parser = argparse.ArgumentParser(
        description="使用 blastp 比对查询蛋白与非哺乳动物蛋白数据库"
    )
    parser.add_argument("-q", "--query", required=True,
                        help="步骤1生成的查询蛋白 FASTA 文件")
    parser.add_argument("-d", "--db", required=True,
                        help="包含非哺乳动物蛋白序列的 FASTA 文件，用于构建 BLAST 数据库")
    parser.add_argument("-o", "--out", default="blast_results.tsv",
                        help="输出 BLAST 结果文件 (default: blast_results.tsv)")
    parser.add_argument("--evalue", default="1e-10",
                        help="E-value 阈值 (default: 1e-10)")
    parser.add_argument("--perc_identity", default="50",
                        help="最小百分比身份 (default: 50)")
    parser.add_argument("--qcov_hsp_perc", default="80",
                        help="最小比对覆盖率 (default: 80)")
    parser.add_argument("--threads", default="4",
                        help="使用的线程数 (default: 4)")
    return parser.parse_args()

def make_blast_db(db_fasta):
    """
    检查 BLAST 数据库文件是否存在，如不存在则构建数据库
    需要的文件后缀：.pin, .phr, .psq
    """
    db_files = [db_fasta + ext for ext in [".pin", ".phr", ".psq"]]
    if not all(os.path.exists(f) for f in db_files):
        print("BLAST 数据库不存在，正在构建...")
        cmd = ["makeblastdb", "-in", db_fasta, "-dbtype", "prot"]
        subprocess.run(cmd, check=True)
    else:
        print("检测到 BLAST 数据库，跳过构建。")

def run_blast(query, db, out, evalue, perc_identity, qcov_hsp_perc, threads):
    """
    调用 blastp 进行比对，比对结果以 tabular 格式输出，字段：
    qseqid sseqid evalue pident qcovs
    使用 --threads 参数指定多线程数
    """
    cmd = [
        "blastp",
        "-query", query,
        "-db", db,
        "-out", out,
        "-evalue", evalue,
        "-perc_identity", perc_identity,
        "-qcov_hsp_perc", qcov_hsp_perc,
        "-num_threads", threads,
        "-outfmt", "6 qseqid sseqid evalue pident qcovs"
    ]
    print("执行 BLAST 命令：\n", " ".join(cmd))
    subprocess.run(cmd, check=True)
    print("BLAST 比对结果已写入：", out)

def main():
    args = parse_args()
    make_blast_db(args.db)
    run_blast(args.query, args.db, args.out, args.evalue, args.perc_identity, args.qcov_hsp_perc, args.threads)

if __name__ == "__main__":
    main()

