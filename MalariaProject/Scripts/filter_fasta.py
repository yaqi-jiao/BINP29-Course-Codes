import sys

def filter_fasta(id_file, fasta_file, output_file):
    # 1. 读取需要移除的 ID
    try:
        with open(id_file, 'r') as f:
            # 使用 set 提高查找速度
            remove_ids = {line.strip() for line in f if line.strip()}
    except FileNotFoundError:
        print(f"错误: 找不到 ID 文件 '{id_file}'")
        return

    # 2. 流式处理 FASTA 文件
    try:
        with open(fasta_file, 'r') as f_in, open(output_file, 'w') as f_out:
            keep_sequence = True
            count_removed = 0
            
            for line in f_in:
                if line.startswith('>'):
                    # 只要 Header 包含任一 ID，就标记为不保留
                    if any(cid in line for cid in remove_ids):
                        keep_sequence = False
                        count_removed += 1
                    else:
                        keep_sequence = True
                    
                    if keep_sequence:
                        f_out.write(line)
                else:
                    # 如果当前序列块标记为保留，则写入
                    if keep_sequence:
                        f_out.write(line)
            
            print(f"处理完成！")
            print(f"已过滤掉 {count_removed} 条序列，结果保存至: {output_file}")
            
    except FileNotFoundError:
        print(f"错误: 找不到 FASTA 文件 '{fasta_file}'")

if __name__ == "__main__":
    # 检查命令行参数数量
    # sys.argv[0] 是脚本名，1-3 是你传入的参数
    if len(sys.argv) != 4:
        print("用法提示: python filter_fasta.py <ID列表文件> <输入FASTA> <输出FASTA>")
        print("例如: python ../../Scripts/filter_fasta.py scaffolds.txt ../Ht.fna filtered_Ht.fna")
    else:
        filter_fasta(sys.argv[1], sys.argv[2], sys.argv[3])
