#!/bin/bash

# --- 配置 ---
# 默认的TSV文件路径 (如果命令行没有提供)
DEFAULT_TSV_FILE="blast_results.tsv"
# 默认提取的Top行数 (如果命令行没有提供)
DEFAULT_TOP_N=10

# --- 函数定义 ---
usage() {
  echo "用法: $0 <tsv_file> <top_n>"
  echo "   <tsv_file>: BLASTP输出的TSV文件路径。"
  echo "   <top_n>   : 要处理的文件前N行。"
  echo ""
  echo "如果没有提供参数，脚本将尝试使用默认值："
  echo "   TSV文件: $DEFAULT_TSV_FILE"
  echo "   Top N  : $DEFAULT_TOP_N"
  exit 1
}

# --- 参数处理 ---
TSV_FILE="${1:-$DEFAULT_TSV_FILE}" # 如果$1未设置或为空，则使用默认值
TOP_N="${2:-$DEFAULT_TOP_N}"      # 如果$2未设置或为空，则使用默认值

# --- 输入验证 ---
# 检查行数是否为正整数
if ! [[ "$TOP_N" =~ ^[1-9][0-9]*$ ]]; then
  echo "错误: <top_n> ('$TOP_N') 必须是一个正整数。" >&2
  usage
fi

# 检查文件是否存在且可读
if [ ! -f "$TSV_FILE" ]; then
  echo "错误: 文件 '$TSV_FILE' 不存在。" >&2
  usage
fi
if [ ! -r "$TSV_FILE" ]; then
  echo "错误: 文件 '$TSV_FILE' 不可读。" >&2
  usage
fi

# --- 核心逻辑 ---
# 1. 使用 head 读取文件的前 N 行
# 2. 使用 cut 提取前两列 (字段1和字段2)，分隔符为Tab (默认)
# 3. 使用 tr 将Tab替换为换行符，使每个ID占一行
# 4. 使用 sort -u 对所有ID进行排序并去重
# 5. 使用 awk 构建Python列表格式的字符串
#    - BEGIN { printf "[" }: 在开始时打印左方括号
#    - NR > 1 { printf ", " }: 对于第二行及之后的行，在打印ID前先打印逗号和空格
#    - { printf "\047%s\047", $0 }: 打印当前行(ID)，并用单引号包裹 (\047 是单引号的八进制表示)
#    - END { print "]" }: 在处理完所有行后打印右方括号
# 6. 将最终结果用 echo 输出 (虽然awk已经打印了，但为了明确是脚本的最终输出，可以保留echo)

unique_ids_list=$(head -n "$TOP_N" "$TSV_FILE" | cut -f 1,2 | tr '\t' '\n' | sort -u | awk '
BEGIN { printf "[" }
NR > 1 { printf ", " }
{ printf "\047%s\047", $0 }
END { print "]" }
')

# --- 输出结果 ---
echo "$unique_ids_list"

exit 0

