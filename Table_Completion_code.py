import os
import sys
import subprocess

def install(package):
    subprocess.check_call([sys.executable, "-m", "pip", "install", package])

try:
    import pandas as pd
except ImportError:
    install("pandas")
    import pandas as pd

try:
    import openpyxl
except ImportError:
    install("openpyxl")
    import openpyxl

def main():
    # 定义输入和输出文件路径
    file_path = "/Users/shenyz/Downloads/2025-7-23_py/human_geckov2_library_b_09mar2015 (2).xlsx"
    output_file = "final_human_geckov2_library_b_09mar2015 (2).xlsx"

    # 工作目录和文件检查
    print(f"当前工作目录: {os.getcwd()}")
    if not os.path.exists(file_path):
        sys.exit(f"找不到输入文件：{file_path}")
    print(f"找到输入文件: {file_path}")

    # 读取 Excel
    df = pd.read_excel(file_path)

    # 标准化列名（去除前后隐藏空格）
    df.columns = df.columns.str.strip()


    # 必需列及默认填充
    required_columns = [
        'sgRNA_Seq', 'sgRNA_Type',
        'Gene_ID', 'Gene_Name', 'Transcript_ID', 'sgRNA_ID'  # 包含重复检查列
    ]

    # 检查缺失列
    missing = [col for col in required_columns if col not in df.columns]
    if missing:
        print(f"错误：缺少以下列：{', '.join(missing)}")
        sys.exit(1)

    # 处理 sgRNA_Type 列
    if 'sgRNA_Type' not in df.columns:
        user_input = input("未检测到列 'sgRNA_Type'。是否该列全部为 Target？(Y/N): ").strip().upper()
        if user_input == 'Y':
            df['sgRNA_Type'] = 'Target'
        else:
            sys.exit("请添加列 'sgRNA_Type'，或默认填充为 Target。")
    else:
        type_col = df['sgRNA_Type'].astype(str).str.strip()
        has_non_targeting = (type_col == 'Non-TargetingControl').any()
        has_target = (type_col == 'Target').any()

        if not has_non_targeting and not has_target:
            user_input = input("检测到列 'sgRNA_Type'，但无 'Non-TargetingControl' 或 'Target' 项。是否默认全部为 Target？(Y/N): ").strip().upper()
            if user_input == 'Y':
                df['sgRNA_Type'] = 'Target'
            else:
                sys.exit("请补全 'sgRNA_Type' 列，或默认填充为 Target。")
        elif has_non_targeting and not has_target:
            df.loc[type_col != 'Non-TargetingControl', 'sgRNA_Type'] = 'Target'
        elif has_target and not has_non_targeting:
            df.loc[type_col != 'Target', 'sgRNA_Type'] = 'Non-TargetingControl'
        else:
            # 若二者都存在但部分缺失，则不处理，保持原样
            pass

    # 补全 sgRNA_ID：对空白或 NaN 值顺延编号
    prefix = "SGR"
    mask_blank = df['sgRNA_ID'].isna() | (df['sgRNA_ID'].astype(str).str.strip() == '')
    # 提取已有的符合格式前缀 + 数字的 ID
    existing_ids = df.loc[~mask_blank, 'sgRNA_ID'].astype(str)
    nums = existing_ids.str.extract(rf"^{prefix}(\d+)$")[0].dropna().astype(int)
    next_idx = nums.max() + 1 if not nums.empty else 1
    # 填充缺失 ID
    for idx in df[mask_blank].index:
        df.at[idx, 'sgRNA_ID'] = f"{prefix}{next_idx:07d}"
        next_idx += 1

    # 删除重复 sgRNA_ID 行（保留首条）
    df = df.drop_duplicates(subset='sgRNA_ID', keep='first')

    # 对 Gene_ID, Gene_Name, Transcript_ID 任意缺失的单元格，填充 'Nan'
    gene_cols = ['Gene_ID', 'Gene_Name', 'Transcript_ID']
    for col in gene_cols:
        # 先把真正的 NaN 填上
        df[col] = df[col].fillna('Nan')
        # 再把空字符串也处理成 'Nan'
        df.loc[df[col].astype(str).str.strip() == '', col] = 'Nan'
    

    # 保存为新文件
    try:
        df.to_excel(output_file, index=False)
        print(f"处理完成，已保存为 {output_file}")
    except Exception as e:
        sys.exit(f"保存文件失败: {e}")


if __name__ == '__main__':
    main()
