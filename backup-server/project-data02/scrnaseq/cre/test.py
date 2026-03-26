import pandas as pd

def test_csv_separator(csv_file):
    try:
        # 尝试读取CSV文件，使用制表符作为分隔符
        df = pd.read_csv(csv_file, delimiter='\t', header=None)
        
        # 打印文件的前几行，检查列数
        print(f"文件 {csv_file} 使用制表符作为分隔符：")
        print(f"每行列数: {df.shape[1]}")
        
        # 查看前几行数据
        print(df.head())
        
    except Exception as e:
        print(f"无法读取文件，出现错误: {e}")
        print("请确保文件使用制表符（\\t）分隔列。")

# 测试函数
csv_file = "output/alignment.csv"  # 请替换成你的CSV文件路径
test_csv_separator(csv_file)

