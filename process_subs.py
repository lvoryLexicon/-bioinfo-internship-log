# process_subs.py
import os
import sys
import re
import pandas as pd

def norm(s):
    return re.sub(r'[^0-9a-z]', '', str(s).lower())

# 映射候选关键词（按优先级）
CANDIDATES = {
    'Gene_ID': ['geneid','gene_id','gene id','geneid'],
    'sgRNA_Seq': ['sgrna_seq','sgrna','sgrna sequence','sgrna_seq','sgrnaseq','sgRNA','sgRNA_seq'],
    'Gene_Name': ['#uniqueid','uniqueid','gene_name','genename','gene name'],
    'Transcript_ID': ['transcript_id','transcriptid','transcript id','chromosome','chr']
}

def map_columns(df):
    cols = list(df.columns)
    nmap = {norm(c): c for c in cols}
    out = {}
    for target, keys in CANDIDATES.items():
        for k in keys:
            for nc, orig in nmap.items():
                if k == nc or k in nc or nc in k:
                    out[target] = orig
                    break
            if target in out:
                break
    return out

def _normalize_seq(s):
    if pd.isna(s):
        return ''
    s = str(s).upper().strip()
    # keep only A C G T U N; convert U->T
    s = re.sub(r'[^ACGTUN]', '', s)
    s = s.replace('U', 'T')
    return s

# 核心处理函数（保持原有行为）
def process_df(df):
    df = df.copy()
    df.columns = df.columns.str.strip()
    required = ['sgRNA_Seq', 'sgRNA_Type', 'Gene_ID', 'Gene_Name', 'Transcript_ID', 'sgRNA_ID']
    missing = [c for c in required if c not in df.columns]

    if 'sgRNA_ID' in missing:
        df['sgRNA_ID'] = pd.NA
        missing = [c for c in missing if c != 'sgRNA_ID']

    if missing:
        if missing == ['sgRNA_Type']:
            df['sgRNA_Type'] = 'Target'
            missing = []
    if missing:
        raise ValueError(f"缺少必要列: {', '.join(missing)}")

    # 统一 sgRNA_Type：非 Non-TargetingControl 的视为 Target（简单、明确）
    df['sgRNA_Type'] = df['sgRNA_Type'].astype(str).str.strip()
    if (df['sgRNA_Type'] == 'Non-TargetingControl').any():
        df.loc[df['sgRNA_Type'] != 'Non-TargetingControl', 'sgRNA_Type'] = 'Target'
    else:
        df['sgRNA_Type'] = 'Target'

    # 补全 sgRNA_ID (格式 SGR + 7 位数字)
    prefix = "SGR"
    mask_blank = df['sgRNA_ID'].isna() | (df['sgRNA_ID'].astype(str).str.strip() == '')
    existing = df.loc[~mask_blank, 'sgRNA_ID'].astype(str)
    nums = existing.str.extract(rf"^{prefix}(\d+)$")[0].dropna().astype(int) if not existing.empty else pd.Series(dtype=int)
    next_idx = int(nums.max()) + 1 if not nums.empty else 1
    for idx in df[mask_blank].index:
        df.at[idx, 'sgRNA_ID'] = f"{prefix}{next_idx:07d}"
        next_idx += 1

    df = df.drop_duplicates(subset='sgRNA_ID', keep='first')

    for col in ['Gene_ID', 'Gene_Name', 'Transcript_ID']:
        if col in df.columns:
            df[col] = df[col].fillna('Nan')
            df.loc[df[col].astype(str).str.strip() == '', col] = 'Nan'
    return df

def find_seq_column(df):
    """
    在 NC 表中自动识别哪一列是 sgRNA-like 序列（仅 ACGTUN 等，长度典型在 >10）
    选择匹配率最高且满足阈值的列；若只有一列则直接返回它。
    """
    best_col = None
    best_score = 0.0
    rows = len(df)
    if rows == 0:
        return None
    for col in df.columns:
        vals = df[col].astype(str).fillna('').str.strip()
        if vals.eq('').all():
            continue
        match = vals.str.match(r'^[ACGTUNacgtun]+$')
        score = match.sum() / rows
        avg_len = vals[match].str.len().mean() if match.sum() > 0 else 0
        if score > best_score and avg_len >= 8:
            best_score = score
            best_col = col
    if best_col is None and len(df.columns) == 1:
        return df.columns[0]
    return best_col if best_score >= 0.4 else None

def main(path):
    print(f"开始处理文件：{os.path.basename(path)}")
    if not os.path.exists(path):
        print("错误：找不到输入文件。"); sys.exit(1)
    
    book = pd.ExcelFile(path)
    sheets = {s.lower(): s for s in book.sheet_names}

    base, ext = os.path.splitext(os.path.basename(path))
    out_base = f"final_{base}{ext}"

    # 情况1：没有 CRISPRa 或 NC，直接读主表并处理
    if 'crispra' not in sheets:
        print("未检测到'CRISPRa'子表，将处理主表。")
        df = pd.read_excel(path)
        print("正在处理数据...")
        df_processed = process_df(df)
        df_processed.to_excel(out_base, index=False)
        print(f"处理完成，结果已保存至：{out_base}")
        return

    # 如果存在 CRISPRa：读取并处理；若存在 NC 则把 NC 的序列行直接追加到 CRISPRa 表中
    sheet_name = sheets['crispra']
    print(f"找到 'CRISPRa' 子表，正在读取：'{sheet_name}'")
    df_crispr = pd.read_excel(path, sheet_name=sheet_name)
    
    # 自动识别并重命名列
    mapping = map_columns(df_crispr)
    rename_map = {}
    for target, orig in mapping.items():
        if orig in df_crispr.columns and target not in df_crispr.columns:
            rename_map[orig] = target
    if rename_map:
        df_crispr = df_crispr.rename(columns=rename_map)
        print("成功映射并重命名列。")
    
    # 如果没有 sgRNA_Seq 列，试图找 'sg' 之类（保底）
    if 'sgRNA_Seq' not in df_crispr.columns:
        candidates = [c for c in df_crispr.columns if 'sg' in norm(c)]
        if candidates:
            df_crispr = df_crispr.rename(columns={candidates[0]: 'sgRNA_Seq'})
            print(f"自动将 '{candidates[0]}' 识别为 sgRNA 序列列。")
        else:
            raise ValueError("CRISPRa 表中未识别到 sgRNA 序列列，请确认列名包含 'sgRNA' 或类似关键字。")

    # 提取 CRISPRa 的 sgRNA 序列集
    crispr_seqs = set(df_crispr['sgRNA_Seq'].dropna().astype(str).apply(_normalize_seq))

    # 若存在 NC 子表：提取 NC 的序列列，标准化后追加
    if 'nc' in sheets:
        print(f"找到 'NC' 子表，正在读取：'{sheets['nc']}'")
        df_nc = pd.read_excel(path, sheet_name=sheets['nc'])
        seq_col = find_seq_column(df_nc)
        seqs = []
        if seq_col is not None:
            seqs = df_nc[seq_col].astype(str).fillna('').apply(_normalize_seq)
            seqs = seqs[seqs.str.len() > 0].unique().tolist()
        else:
            if df_nc.shape[1] == 1:
                seqs = df_nc.iloc[:,0].astype(str).fillna('').apply(_normalize_seq)
                seqs = seqs[seqs.str.len() > 0].unique().tolist()
        
        # 检查重复序列
        nc_seqs = set(seqs)
        common_seqs = crispr_seqs.intersection(nc_seqs)
        if common_seqs:
            print(f"\n警告：在 'CRISPRa' 和 'NC' 表中发现 {len(common_seqs)} 个重复的 sgRNA 序列。")
            print("部分重复序列：", list(common_seqs)[:5], "\n")
        
        if seqs:
            print(f"正在从 'NC' 表追加 {len(seqs)} 条非靶向sgRNA。")
            cols = list(df_crispr.columns)
            if 'sgRNA_Type' not in cols:
                df_crispr['sgRNA_Type'] = pd.NA
                cols = list(df_crispr.columns)
            
            rows = []
            for s in seqs:
                row = {c: pd.NA for c in cols}
                row['sgRNA_Seq'] = s
                row['sgRNA_Type'] = 'Non-TargetingControl'
                rows.append(row)
            df_append = pd.DataFrame(rows, columns=cols)
            df_crispr = pd.concat([df_crispr, df_append], ignore_index=True)
        else:
            print("未在 'NC' 表中发现有效的 sgRNA 序列，跳过追加。")

    # 对现存的（非 NC 追加）行，如果没有显式标为 Non-TargetingControl，统一标为 Target
    if 'sgRNA_Type' not in df_crispr.columns:
        df_crispr['sgRNA_Type'] = 'Target'
    else:
        df_crispr['sgRNA_Type'] = df_crispr['sgRNA_Type'].fillna('')
        df_crispr.loc[df_crispr['sgRNA_Type'] != 'Non-TargetingControl', 'sgRNA_Type'] = 'Target'
    
    print("正在对所有数据进行最终处理...")
    processed = process_df(df_crispr)
    processed.to_excel(out_base, index=False)
    print(f"CRISPRa（含NC追加）处理完成，结果已保存至：{out_base}")

if __name__ == '__main__':
    if len(sys.argv) < 2:
        print("用法: python process_subs.py <输入文件.xlsx>")
        sys.exit(1)
    main(sys.argv[1])