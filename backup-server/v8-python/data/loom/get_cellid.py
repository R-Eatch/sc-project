#!/usr/bin/env python
import scanpy as sc

def main():
    # 定义 S-SG 类型的 loom 文件字典
    s_sg_dict = {
        'S-SG-1MTH-1': 'S-SG-1MTH-1.loom',
        'S-SG-1MTH-2': 'S-SG-1MTH-2.loom',
        'S-SG-8M-2': 'S-SG-8M-2.loom',
        'S-SG-8M-3': 'S-SG-8M-3.loom',
        'S-SG-GES14': 'S-SG-GES14.loom',
        'S-SG-LA-1': 'S-SG-LA-1.loom',
        'S-SG-LA-2': 'S-SG-LA-2.loom',
        'S-SG-P20': 'S-SG-P20.loom'
    }
    
    # 定义 M-SG 类型的 loom 文件字典
    m_sg_dict = {
        'M-SG-3WK-1': 'M-SG-3WK-1.loom',
        'M-SG-3WK-2': 'M-SG-3WK-2.loom',
        'M-SG-8WK-1': 'M-SG-8WK-1.loom',
        'M-SG-8WK-2': 'M-SG-8WK-2.loom',
        'M-SG-E16_5': 'M-SG-E16_5.loom',
        'M-SG-GES13_5': 'M-SG-GES13_5.loom',
        'M-SG-GES16_5': 'M-SG-GES16_5.loom',
        'M-SG-LA-1': 'M-SG-LA-1.loom',
        'M-SG-LA-2': 'M-SG-LA-2.loom',
        'M-SG-P1': 'M-SG-P1.loom'
    }
    
    # 构造最终结果字典
    result_dict = {}
    
    # 处理 S-SG 类型文件
    for sample_name, filename in s_sg_dict.items():
        # 读取 loom 文件得到 AnnData 对象
        ldata2 = sc.read_loom(filename)
        # 从 obs_names 第一个元素中以冒号分割后取第一个部分作为值
        value = ldata2.obs_names[0].split(":")[0]
        result_dict[sample_name] = value

    # 处理 M-SG 类型文件
    for sample_name, filename in m_sg_dict.items():
        ldata2 = sc.read_loom(filename)
        value = ldata2.obs_names[0].split(":")[0]
        result_dict[sample_name] = value

    # 打印最终构建的字典
    print(result_dict)

if __name__ == '__main__':
    main()

