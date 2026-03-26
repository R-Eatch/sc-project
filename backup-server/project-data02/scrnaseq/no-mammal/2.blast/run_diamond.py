#!/usr/bin/env python
"""
脚本功能：
1. 使用输入的非哺乳动物蛋白 FASTA 文件构建 Diamond 数据库，
2. 利用 Diamond 的 blastp 模式，比对步骤1生成的查询蛋白序列（query_proteins.fasta）与该数据库，
   使用参数 '--ultra-sensitive'，并设置 E-value 阈值（默认 1e-5）过滤不可靠的比对结果，
3. 结果以 tabular 格式输出（字段：qseqid, sseqid, evalue, pident）。

示例文献参数：
  --ultra-sensitive
  并过滤掉 E-value > 1e-5 的比对结果
"""

import argparse
import subprocess
import os

def parse_args():
    parser = argparse.ArgumentParser(
        description="使用 Diamond 对查询蛋白与非哺乳动物蛋白数据库进行比对"
    )
    parser.add_argument("-q", "--query", required=True,
                        help="查询蛋白 FASTA 文件（如 extract_longest_protein.py 生成的文件）")
    parser.add_argument("-d", "--db", required=True,
                        help="非哺乳动物蛋白 FASTA 文件，用于构建 Diamond 数据库")
    parser.add_argument("-o", "--out", default="diamond_results.tsv",
                        help="输出 Diamond 比对结果文件 (default: diamond_results.tsv)")
    parser.add_argument("--evalue", default="1e-5",
                        help="E-value 阈值，过滤掉大于该值的比对 (default: 1e-5)")
    parser.add_argument("--threads", default="8",
                        help="使用的线程数 (default: 8)")
    return parser.parse_args()

def make_diamond_db(db_fasta, db_prefix):
    """
    检查 Diamond 数据库是否已构建，如果未构建则使用 diamond makedb 构建数据库。
    构建后的数据库前缀为 db_prefix
    """
    # Diamond 构建的数据库会生成一个 .dmnd 文件
    db_file = db_prefix + ".dmnd"
    if not os.path.exists(db_file):
        print("Diamond 数据库不存在，正在构建...")
        cmd = [
            "diamond", "makedb",
            "--in", db_fasta,
            "-d", db_prefix
        ]
        subprocess.run(cmd, check=True)
    else:
        print("检测到 Diamond 数据库，跳过构建。")
    return db_file

def run_diamond(query, db_prefix, out, evalue, threads):
    """
    使用 Diamond blastp 进行比对，使用 '--ultra-sensitive' 参数。
    输出格式：tabular，字段：qseqid, sseqid, evalue, pident
    """
    cmd = [
        "diamond", "blastp",
        "--query", query,
        "--db", db_prefix,
        "--out", out,
        "--ultra-sensitive",
        "--evalue", evalue,
        "--threads", threads,
        "--outfmt", "6",
	"--max-target-seqs","10"
    ]
    print("执行 Diamond 命令：\n", " ".join(cmd))
    subprocess.run(cmd, check=True)
    print("Diamond 比对结果已写入：", out)

def main():
    args = parse_args()
    # 使用非哺乳动物蛋白 FASTA 文件构建数据库，数据库前缀以 db 文件名为基础
    db_prefix = os.path.splitext(args.db)[0] + "_dmnd"
    make_diamond_db(args.db, db_prefix)
    run_diamond(args.query, db_prefix, args.out, args.evalue, args.threads)

if __name__ == "__main__":
    main()

