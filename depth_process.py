#!/usr/bin/env python3
"""
plot_bigwig_coverage_tracks.py

从 bigWig(depth.bw) 读取每条染色体的覆盖度，不做 log10 变换并把每条染色体绘制成一条水平 track（带填充）。
输出图像 (PNG / PDF)。

依赖:
  pip install pyBigWig numpy matplotlib tqdm

使用:
  python plot_bigwig_coverage_tracks.py \
      --bw /Users/shenyz/Downloads/depth_process/depth.bw \
      --out coverage_tracks.png
"""

import argparse
import os
import math
import pyBigWig
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
from tqdm import tqdm

# ---- 辅助函数 ----
def canonical_chr_name(chrom):
    """把染色体名字标准化为 'chrN' 格式（如果输入已经是 '1' 则返回 'chr1'）"""
    if chrom.lower().startswith("chr"):
        return chrom
    else:
        return "chr" + chrom

def bin_array_mean(arr, n_bins):
    """
    将 1D 数组 arr 分到 n_bins 个 bin 上，每个 bin 返回平均值（忽略 nan）。
    返回长度为 n_bins 的 numpy 数组。
    """
    if len(arr) == 0:
        return np.zeros(n_bins, dtype=float)
    # 使用等长切分（最后一个 bin 可能稍长）
    # 将 arr 切成 n_bins 大致相等的片段并取每段平均
    split = np.array_split(arr, n_bins)
    binned = []
    for seg in split:
        if seg.size == 0:
            binned.append(0.0)
        else:
            # 忽略 nan
            seg = seg[~np.isnan(seg)]
            if seg.size == 0:
                binned.append(0.0)
            else:
                binned.append(float(np.mean(seg)))
    return np.array(binned, dtype=float)

def read_and_bin_chrom(bw, chrom, n_bins):
    """从 bigWig 对象读取某条染色体全部值并做 bin 平均，返回 binned ndarray"""
    length = bw.chroms()[chrom]
    # pyBigWig.values 会返回长度 = length 的列表，nan 表示没有数据
    vals = np.array(bw.values(chrom, 0, length), dtype=float)
    # 将 nan 替换为 0
    vals = np.nan_to_num(vals, nan=0.0)
    if len(vals) <= n_bins:
        # 如果序列长度小于 bin 数，直接 pad
        if len(vals) == 0:
            return np.zeros(n_bins, dtype=float)
        result = np.zeros(n_bins, dtype=float)
        result[:len(vals)] = vals
        return result
    else:
        return bin_array_mean(vals, n_bins)

# ---- 主绘图函数 ----
def plot_tracks_from_bigwig(bw_path,
                            out_path,
                            chrom_order=None,
                            include_chrs=None,
                            base_bins=500,  # 减少基准bin数，避免图像过大
                            figsize=(10, 8),  # 调整默认尺寸
                            dpi=100,  # 降低DPI
                            log_transform=False,
                            title="Distribution of read coverage"):
    """
    bw_path: bigWig 文件路径
    out_path: 输出图像路径
    chrom_order: list of chromosome names preferred order (like ['chr1', 'chr2', ...])
    include_chrs: set/list of chromosomes to include (if None 则使用 bw 中所有)
    base_bins: 基准bin数，实际bin数会根据染色体长度按比例调整
    """
    if not os.path.exists(bw_path):
        raise FileNotFoundError(f"bigWig file not found: {bw_path}")

    bw = pyBigWig.open(bw_path)
    if not bw:
        raise IOError("Failed to open bigWig file")

    # 从 bigWig 获取染色体及长度
    bw_chroms = bw.chroms()  # dict: {chrom: length}
    chroms_in_bw = list(bw_chroms.keys())

    # 标准染色体顺序
    default_order = [f"chr{i}" for i in range(1, 23)] + ["chrX", "chrY"]
    chroms_norm = {c: canonical_chr_name(c) for c in chroms_in_bw}
    
    # 生成最终使用的染色体顺序
    final_list = []
    if chrom_order:
        for c in chrom_order:
            cn = canonical_chr_name(c)
            if cn in bw_chroms:
                final_list.append(cn)
    else:
        for c in default_order:
            if c in bw_chroms:
                final_list.append(c)
        for c in chroms_in_bw:
            if c not in final_list:
                final_list.append(c)

    # 如果用户指定 include_chrs，进行过滤
    if include_chrs:
        inc = set(canonical_chr_name(x) for x in include_chrs)
        final_list = [c for c in final_list if c in inc]

    if not final_list:
        final_list = chroms_in_bw

    # 计算最长染色体的长度，用于确定bin数比例
    max_length = max(bw_chroms[chrom] for chrom in final_list)

    # 读取并处理每条染色体（按长度比例确定bin数）
    tracks = []
    print("Reading and binning chromosomes:")
    for chrom in tqdm(final_list):
        try:
            # 根据染色体长度比例确定bin数，但设置上限避免图像过大
            chrom_length = bw_chroms[chrom]
            n_bins = max(10, min(2000, int(base_bins * chrom_length / max_length)))  # 设置上限2000
            binned = read_and_bin_chrom(bw, chrom, n_bins=n_bins)
        except Exception as e:
            print(f"Warning: failed to read {chrom}: {e}")
            n_bins = max(10, min(2000, int(base_bins * 1000000 / max_length)))  # 假设长度1000000
            binned = np.zeros(n_bins, dtype=float)

        # ---------- 在这里加入 log10(x+1) 变换 ----------
        if log_transform:
            binned = np.log10(binned + 1.0)
        # ----------------------------------------------

        tracks.append((chrom, binned, chrom_length))

    bw.close()

    # ---- 绘图 ----
    n_tracks = len(tracks)
    track_height = 1.0
    total_height = n_tracks * track_height

    # 动态调整图像高度
    fig_h = min(20, max(6, n_tracks * 0.8))  # 根据染色体数量调整高度
    fig_w = 10  # 固定宽度
    
    fig = plt.figure(figsize=(fig_w, fig_h), dpi=dpi)
    ax = fig.add_subplot(1, 1, 1)

    # 计算最大覆盖度用于统一y轴范围
    max_coverage_all = max(np.nanmax(vals) for _, vals, _ in tracks) if tracks else 1.0
    
    # 为每条 track 计算 y 偏移（从上到下）
    max_bins = 0
    for i, (chrom, vals, chrom_length) in enumerate(tracks):
        y_offset = total_height - (i + 1) * track_height

        # 使用实际的bin数作为x轴
        x = np.arange(len(vals))
        max_bins = max(max_bins, len(vals))

        # ---------- 每条染色体按自身最大值归一化（与之前的实现一致） ----------
        # 这里使用每条染色体的最大值作为归一化基准，避免某条极值染色体把其他染色体压扁
        per_chrom_max = np.nanmax(vals) if np.nanmax(vals) > 0 else 1.0
        norm_vals = vals / per_chrom_max * (track_height * 0.9)
        # -------------------------------------------------------------------------

        # 绘制填充区域
        ax.fill_between(x, y_offset, y_offset + norm_vals, step='mid', alpha=0.9, color='navy')

        # 在左侧添加覆盖度标签（0和最大值）
        max_val = np.nanmax(vals) if np.nanmax(vals) > 0 else 0
        ...

        # 格式化数值：如果大于1就取整，否则保留1位小数
        if max_val >= 1:
            max_label = f"{max_val:.0f}"
        else:
            max_label = f"{max_val:.1f}"
            
        # 在track左侧添加0和最大值标签
        label_x_pos = -max_bins * 0.02  # 左侧位置
        ax.text(label_x_pos, y_offset + track_height*0.1, '0', 
                va='bottom', ha='right', fontsize=8, color='black')  # 改为黑色
        ax.text(label_x_pos, y_offset + track_height*0.8, max_label, 
                va='top', ha='right', fontsize=8, color='black')  # 改为黑色
        
        # 在右侧写染色体名字并用方框框起来
        chrom_label_x = max_bins + max_bins * 0.02  # 右侧位置
        # 添加方框
        box_width = max_bins * 0.07
        box_height = track_height * 0.4
        box_x = chrom_label_x - box_width * 0.5
        box_y = y_offset + (track_height - box_height) * 0.5
        rect = Rectangle((box_x, box_y), box_width, box_height, 
                        linewidth=1, edgecolor='black', facecolor='white')
        ax.add_patch(rect)
        # 添加染色体名字（黑色，居中对齐）
        ax.text(chrom_label_x, y_offset + track_height*0.5, chrom,
                va='center', ha='center', fontsize=9, color='black')

    # 在最左侧添加竖着的"Read Density"标签（黑色）
    ax.text(-max_bins*0.15, total_height/2, 'Read Density', 
            rotation=90, va='center', ha='center', fontsize=12, color='black')

    # 图形美化
    # 调整x轴位置，下移一点点（增加底部边距）
    fig.subplots_adjust(bottom=0.13)
    ax.set_xlim(-max_bins*0.0, max_bins + max_bins*0.15)  # 给右侧标签留更多空间
    ax.set_ylim(0, total_height)
    ax.set_yticks([])  # 隐藏左侧 y 轴刻度
    ax.set_xlabel("Genomic position (binned by chromosome length)", color='black')
    ax.set_title(title, color='black')
    
    # 设置坐标轴标签颜色为黑色
    ax.tick_params(axis='x', colors='black')
    
    # 移除所有脊柱（边框），只保留底部
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['left'].set_visible(False)
    ax.spines['bottom'].set_visible(True)
    ax.spines['bottom'].set_color('black')  # 底部轴线设为黑色

    # 保存文件
    fig.tight_layout()
    outdir = os.path.dirname(out_path)
    if outdir and not os.path.exists(outdir):
        os.makedirs(outdir)
    
    # 检查图像尺寸是否过大
    width_pixels = int(fig_w * dpi)
    height_pixels = int(fig_h * dpi)
    if width_pixels > 10000 or height_pixels > 10000:
        print(f"Warning: Image size ({width_pixels}x{height_pixels}) is large, reducing DPI...")
        dpi = min(dpi, 72)  # 进一步降低DPI
    
    fig.savefig(out_path, bbox_inches='tight', dpi=dpi)
    print(f"Saved figure to {out_path}")
    plt.close(fig)


# ---- CLI 支持 ----
def parse_args():
    p = argparse.ArgumentParser(description="Plot coverage tracks from bigWig")
    p.add_argument("--bw", required=True, help="Path to bigWig file (e.g. depth.bw)")
    p.add_argument("--out", required=True, help="Output image path (png/pdf)")
    p.add_argument("--bins", type=int, default=500, help="Base number of bins for longest chromosome (default: 500)")
    p.add_argument("--figwidth", type=float, default=10, help="Figure width (inches)")
    p.add_argument("--figheight", type=float, default=8, help="Figure height (inches)")
    p.add_argument("--dpi", type=int, default=100, help="Output DPI (default: 100)")
    p.add_argument("--log", dest="log", action="store_true", help="Enable log10(x+1) transform (default: disabled)")
    return p.parse_args()

if __name__ == "__main__":
    args = parse_args()
    chrom_order = [f"chr{i}" for i in range(1, 23)] + ["chrX", "chrY"]
    plot_tracks_from_bigwig(
        bw_path=args.bw,
        out_path=args.out,
        chrom_order=chrom_order,
        base_bins=args.bins,
        figsize=(args.figwidth, args.figheight),
        dpi=args.dpi,
        log_transform=args.log
    )
