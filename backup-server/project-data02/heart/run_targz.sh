#!/usr/bin/env bash
set -euo pipefail

# 用法：
#   ./pack_result.sh [result目录路径] [输出tar.gz文件名]
# 例如：
#   ./pack_result.sh result
#   ./pack_result.sh /path/to/result my_pack.tar.gz

DIR=${1:-result}
OUT=${2:-"${DIR%/}_csv_png_pdf_$(date +%Y%m%d_%H%M%S).tar.gz"}

if [[ ! -d "$DIR" ]]; then
  echo "目录不存在：$DIR" >&2
  exit 1
fi

parent_dir="$(cd "$(dirname "$DIR")" && pwd)"
base_dir="$(basename "$DIR")"

# 临时文件保存以 NUL 分隔的文件清单
filelist="$(mktemp)"
trap 'rm -f "$filelist"' EXIT

(
  cd "$parent_dir"
  find "$base_dir" -type f \( -iname '*.csv' -o -iname '*.png' -o -iname '*.pdf' \) -print0 \
    > "$filelist"
)

# 注意 -C 必须放在 --files-from 之前
tar -C "$parent_dir" --null -czf "$OUT" --files-from "$filelist"

echo "完成：$OUT"

