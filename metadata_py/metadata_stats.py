"""
metadata_stats.py


用法示例：
    from metadata_stats import plot_metadata_stats
    stats, seg_stats = plot_metadata_stats(
        "/Users/shenyz/Downloads/metadata/metadata.csv",
        sep=None,            # None -> pandas 自动推断（csv 或 tsv 都可以）
        top_n=10,
        save_path="metadata_stats.png",
        show=True
    )

依赖：pandas, matplotlib
"""

from typing import Optional, Tuple
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt


def plot_metadata_stats(path: str,
                        sep: Optional[str] = None,
                        top_n: int = 10,
                        save_path: Optional[str] = None,
                        show: bool = True) -> Tuple[dict, pd.DataFrame]:
    """
    读取 metadata 文件，计算统计量并绘图（左侧文本统计，右侧分割方法饼图），
    背景为黑色、文本为白色。

    参数:
        path: metadata 文件路径（csv 或 tsv 都可以）。
        sep: 传给 pandas.read_csv 的 sep；若为 None 将自动推断（使用 engine='python'）。
        top_n: 只显示分组计数的前 top_n 项。
        save_path: 若给定，则把绘图保存为该文件（PNG/PDF 均可）。
        show: 是否在函数中调用 plt.show()。

    返回:
        stats: 一个字典，包含 total_cells, median/mean 等值。
        seg_stats: 一个 DataFrame，包含 segmentation_method, n, percent（按 n 降序，最多 top_n 行）。
    """

    # 读取数据（兼容 csv 或 tsv）
    try:
        meta = pd.read_csv(path, sep=sep, engine='python')
    except Exception as e:
        raise RuntimeError(f"无法读取文件 {path}: {e}")

    # 强制转换为数值并计算统计量
    for col in ["nCount_Xenium", "nFeature_Xenium"]:
        if col in meta.columns:
            meta[col] = pd.to_numeric(meta[col], errors='coerce')
        else:
            # 如果列缺失，创建全是 NaN 的列以避免后续 KeyError
            meta[col] = np.nan

    stats = {
        "total_cells": int(len(meta)),
        "median_nCount": float(np.nanmedian(meta['nCount_Xenium'])),
        "median_nFeature": float(np.nanmedian(meta['nFeature_Xenium'])),
        "mean_nCount": float(np.nanmean(meta['nCount_Xenium'])),
        "mean_nFeature": float(np.nanmean(meta['nFeature_Xenium']))
    }

    # segmentation_method 计数与百分比（包含 NaN 作为一个类别）
    seg_counts = (meta['segmentation_method']
                  .fillna('<NA>')
                  .value_counts(dropna=False)
                  .rename_axis('segmentation_method')
                  .reset_index(name='n'))

    seg_counts['percent'] = seg_counts['n'] / seg_counts['n'].sum() * 100
    seg_stats = seg_counts.sort_values('n', ascending=False).head(top_n).reset_index(drop=True)

    # 开始绘图（黑底，白字）
    fig, axes = plt.subplots(1, 2, figsize=(12, 6))
    fig.patch.set_facecolor('black')
    for ax in axes:
        ax.set_facecolor('black')

    # 左侧文本块
    left_ax = axes[0]
    left_ax.axis('off')

    lines = [
        f"Number of cells detected: {stats['total_cells']:,}",
        f"Median nCount: {stats['median_nCount']:.2f}",
        f"Median nFeature: {stats['median_nFeature']:.2f}",
        f"Average nCount: {stats['mean_nCount']:.2f}",
        f"Average nFeature: {stats['mean_nFeature']:.2f}"
    ]

    # 从上到下绘制文本，保持左对齐
    y_start = 0.9
    y_step = 0.16
    for i, line in enumerate(lines):
        left_ax.text(0.01, y_start - i * y_step, line, transform=left_ax.transAxes,
                     fontsize=12, color='white', va='top', ha='left')

    left_ax.set_title('Metadata Statistics', color='white', pad=10)

    # 右侧饼图
    right_ax = axes[1]
    sizes = seg_stats['n'].values
    labels = seg_stats['segmentation_method'].astype(str).values

    # 为只有百分比大于 1 的扇区显示百分比标签
    percents = seg_stats['percent'].values
    pct_labels = [f"{p:.1f}%" if p > 1 else "" for p in percents]

    # 绘制饼图（不直接使用 autopct，以便更灵活地控制显示）
    wedges, texts = right_ax.pie(sizes, startangle=90, labels=None)

    # 在扇区中心写入百分比（对于小于或等于 1% 的不写）
    for wedge, pct in zip(wedges, pct_labels):
        if pct:
            ang = (wedge.theta2 + wedge.theta1) / 2.0
            x = 0.6 * np.cos(np.deg2rad(ang))
            y = 0.6 * np.sin(np.deg2rad(ang))
            right_ax.text(x, y, pct, ha='center', va='center', color='white', fontsize=9)

    right_ax.set(aspect="equal")
    right_ax.set_title('Segmentation Method (top {})'.format(top_n), color='white')

    # 图例放在右侧
    legend = right_ax.legend(wedges, labels, title='Segmentation Method',
                         bbox_to_anchor=(1.3, 0.5), loc='center left', frameon=False)

    # 将图例文字设为白色
    plt.setp(legend.get_texts(), color='white')
    if legend.get_title():
        plt.setp(legend.get_title(), color='white')

    plt.tight_layout()

    # 保存文件（若指定）
    if save_path:
        try:
            fig.savefig(save_path, facecolor=fig.get_facecolor(), bbox_inches='tight', dpi=150)
        except Exception as e:
            raise RuntimeError(f"保存图片失败: {e}")

    if show:
        plt.show()
    else:
        plt.close(fig)

    return stats, seg_stats


if __name__ == '__main__':
    # 简单的命令行示例（便于直接运行）
    import argparse

    parser = argparse.ArgumentParser(description='Plot metadata statistics')
    parser.add_argument('--path', required=True, help='metadata file path (csv or tsv)')
    parser.add_argument('--sep', default=None, help='delimiter for file (default: auto)')
    parser.add_argument('--top_n', type=int, default=10, help='top N segmentation methods')
    parser.add_argument('--save', default=None, help='save plot to this path')
    args = parser.parse_args()

    plot_metadata_stats(args.path, sep=args.sep, top_n=args.top_n, save_path=args.save, show=True)
