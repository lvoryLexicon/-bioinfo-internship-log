# 导入并安装所需库（仅需安装一次）
import subprocess
import sys

def install(package):
    subprocess.check_call([sys.executable, "-m", "pip", "install", package])

try:
    import pandas as pd
except ImportError:
    install("pandas")
    import pandas as pd

# 读取Excel文件
file_path = 'library.xlsx'
df = pd.read_excel(file_path)

# 标准化列名（去除隐藏空格）
df.columns = df.columns.str.strip()

# 定义所需列并确保存在
required_columns = ['sgRNA_Seq', 'SgRNA_ID', 'sgRNA_Tpye', 'Gene_ID', 'Gene_Name', 'Transcript_ID']
for col in required_columns:
    if col not in df.columns:
        df[col] = pd.NA

# 补全 SgRNA_ID：对空白或 NaN 值顺延生成编号
mask_id_blank = df['SgRNA_ID'].isna() | (df['SgRNA_ID'].astype(str).str.strip() == '')
if mask_id_blank.any():
    prefix = "SGR"
    existing = df.loc[~mask_id_blank, 'SgRNA_ID'].astype(str)
    nums = existing.str.extract(rf"{prefix}(\d+)")[0].dropna().astype(int)
    next_idx = nums.max() + 1 if not nums.empty else 1
    for idx in df[mask_id_blank].index:
        df.at[idx, 'SgRNA_ID'] = f"{prefix}{next_idx:07d}"
        next_idx += 1

# 删除重复 SgRNA_ID 行（保留首条）
df = df.drop_duplicates(subset='SgRNA_ID', keep='first')

# 对第4、5、6列（Gene_ID, Gene_Name, Transcript_ID）同时空白或缺失的行，填充字符串 'Nan'
cols_456 = ['Gene_ID', 'Gene_Name', 'Transcript_ID']
# 构建一个布尔 DataFrame，标记每列是否为空白或缺失
blank_df = df[cols_456].apply(lambda col: col.isna() | (col.astype(str).str.strip() == ''))
# 仅当三列都为 True 时，才填充 'Nan'
mask_all_missing = blank_df.all(axis=1)
df.loc[mask_all_missing, cols_456] = 'Nan'

# 保存为新文件
output_file = 'final_library.xlsx'
df.to_excel(output_file, index=False)
print(f"处理完成，已保存为 {output_file}")
