#!/usr/bin/env python3
"""
mosdepth_plot.py

将 mosdepth 生成的 summary 与 global.dist 文件读取并绘图（左：每染色体 mean depth，右：depth vs percent(>=depth) 曲线）。

用法示例：
    python mosdepth_plot.py \
      --summary genome.mosdepth.summary.txt \
      --globaldist genome.mosdepth.global.dist.txt \
      --out-prefix mosdepth_plots \
      --xlim 200 \
      --sort-by mean \
      --dpi 150 \
      --no-show

依赖:
    pandas, matplotlib, numpy
    pip install pandas matplotlib numpy
"""
from __future__ import annotations
import argparse
import sys
from pathlib import Path
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker


def read_summary(summary_fn: Path) -> pd.DataFrame:
    """读取 summary 文件并标准化 chrom_id，筛选数字和 X/Y 染色体。"""
    try:
        summary = pd.read_csv(
            summary_fn, sep='\t', header=None,
            names=['chrom', 'length', 'bases', 'mean', 'min', 'max'],
            comment='#', dtype={'chrom': str}
        )
    except Exception as e:
        raise RuntimeError(f'无法读取 summary 文件 {summary_fn}: {e}')

    summary['mean'] = pd.to_numeric(summary['mean'], errors='coerce')

    # 从 chrom 列提取染色体编号（兼容 'chr1','1','ChrX' 等），只取数字或 X/Y
    pattern = r'(?i)^(?:chr)?(\d+|x|y)\b'
    summary['chrom_id'] = summary['chrom'].astype(str).str.extract(pattern)[0].str.upper()

    # 构造要保留的染色体：所有出现的数字染色体（按数字排序）加上出现的 X/Y
    nums = sorted({c for c in summary['chrom_id'].unique() if pd.notna(c) and c.isdigit()}, key=lambda x: int(x))
    keep = list(map(str, nums))
    for s in ('X', 'Y'):
        if s in summary['chrom_id'].values:
            keep.append(s)

    summary = summary[summary['chrom_id'].isin(keep)].copy()
    if summary.empty:
        raise RuntimeError('summary 文件中没有找到有效的染色体记录 (chr1..chr22/X/Y)。')

    return summary


def read_global_dist(dist_fn: Path) -> pd.DataFrame:
    """读取 global.dist，自动识别 2/3 列格式并返回 depth/percent DataFrame。"""
    try:
        raw = pd.read_csv(dist_fn, delim_whitespace=True, header=None, comment='#', dtype=str)
    except Exception as e:
        raise RuntimeError(f'无法读取 global.dist 文件 {dist_fn}: {e}')

    # 判定列数并命名
    if raw.shape[1] == 3:
        raw.columns = ['group', 'depth', 'percent']
    elif raw.shape[1] == 2:
        raw.columns = ['depth', 'percent']
    else:
        raw = raw.iloc[:, :3]
        raw.columns = ['group', 'depth', 'percent']

    # 选择合适的 group（优先 total）
    if 'group' in raw.columns and 'total' in raw['group'].values:
        dsel = raw[raw['group'] == 'total'][['depth', 'percent']].copy()
    elif 'group' in raw.columns and raw.shape[1] == 3:
        first_group = raw['group'].iloc[0]
        dsel = raw[raw['group'] == first_group][['depth', 'percent']].copy()
    else:
        dsel = raw[['depth', 'percent']].copy()

    # 类型转换 & 清理
    dsel['depth'] = pd.to_numeric(dsel['depth'], errors='coerce')
    dsel['percent'] = pd.to_numeric(dsel['percent'], errors='coerce')
    dsel = dsel.dropna(subset=['depth', 'percent'])

    if dsel.empty:
        raise RuntimeError('无法从 global.dist 中提取有效 depth/percent 数据。')

    # 如果 percent 看起来是 0-100，则转换为 0-1
    if dsel['percent'].max() > 1.0:
        dsel['percent'] = dsel['percent'] / 100.0

    return dsel


def detect_and_compute_coverage(dsel: pd.DataFrame) -> (pd.Index, pd.Series, str):
    """
    输入 per-depth 或累积的 depth/percent 数据（dsel），
    返回 depth index、percent(>=depth)（以 0-100 表示）和检测说明字符串。
    """
    s = dsel.groupby('depth', as_index=True)['percent'].sum()
    s = s.sort_index()  # 索引为深度，升序

    total = s.sum()
    is_decreasing = s.is_monotonic_decreasing
    is_increasing = s.is_monotonic_increasing

    # 判定逻辑（同你原脚本）
    if abs(total - 1.0) < 0.1:
        coverage = s.iloc[::-1].cumsum().iloc[::-1] * 100.0
        detected = "per-depth -> computed coverage as sum(percent for >=depth)"
    elif is_decreasing and s.iloc[0] >= 0.9:
        coverage = s * 100.0
        detected = "already cumulative (>= depth) -> use directly"
    elif is_increasing and s.iloc[-1] >= 0.9:
        pct_le = s  # cumulative <= depth
        pct_ge = 1.0 - pct_le.shift(fill_value=0.0)
        coverage = pct_ge * 100.0
        detected = "already cumulative (<= depth) -> converted to >= depth"
    else:
        coverage = s.iloc[::-1].cumsum().iloc[::-1] * 100.0
        detected = "ambiguous -> assumed per-depth (fallback)"

    return coverage.index.values, coverage.values, detected


def plot_mosdepth(left_df: pd.DataFrame,
                   depths,
                   coverage_vals,
                   out_prefix: Path,
                   xlim: float | None = 200,
                   sort_by: str = 'mean',
                   figsize=(16, 6),
                   dpi: int = 100,
                   show: bool = True):
    """绘图并保存。"""
    # 根据参数决定染色体排序显示
    if sort_by == 'mean':
        plot_df = left_df.sort_values(by='mean', ascending=False)
    elif sort_by == 'chrom':
        # 按自然染色体顺序（chrom_id 数字再 X Y）
        # 构造 sort key：数字优先按 int，X->100, Y->101
        def chrom_key(x):
            cid = str(x)
            if cid.isdigit():
                return int(cid)
            if cid == 'X':
                return 100
            if cid == 'Y':
                return 101
            return 999
        plot_df = left_df.copy()
        plot_df['sort_key'] = plot_df['chrom_id'].map(chrom_key)
        plot_df = plot_df.sort_values(by='sort_key').drop(columns=['sort_key'])
    else:
        plot_df = left_df.copy()

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=figsize)

    # 左：只显示已筛选的染色体（按 plot_df）
    ax1.bar(plot_df['chrom'], plot_df['mean'])
    ax1.set_xlabel('')
    ax1.set_ylabel('Mean depth')
    ax1.grid(axis='y', linestyle='--', linewidth=0.5)
    plt.setp(ax1.get_xticklabels(), rotation=45, ha='right', fontsize='small')

    # 右：depth vs percent(>= depth)
    ax2.plot(depths, coverage_vals, linestyle='-', marker=None, linewidth=2)
    ax2.set_xlabel('Depth')
    ax2.set_ylabel('Percent (%) ≥ depth')
    ax2.set_ylim(0, 100.5)

    if xlim is not None:
        ax2.set_xlim(0, float(xlim))

    # 不使用科学计数法
    fmt = mticker.ScalarFormatter(useOffset=False)
    fmt.set_scientific(False)
    ax2.yaxis.set_major_formatter(fmt)
    ax2.grid(linestyle='--', linewidth=0.5)

    plt.tight_layout()

    png_fn = out_prefix.with_suffix('.png')
    pdf_fn = out_prefix.with_suffix('.pdf')

    fig.savefig(str(png_fn), dpi=dpi, bbox_inches='tight')
    fig.savefig(str(pdf_fn), format='pdf', bbox_inches='tight')

    if show:
        plt.show()
    else:
        plt.close(fig)


def parse_args():
    p = argparse.ArgumentParser(description='Plot mosdepth summary & global.dist coverage curve.')
    p.add_argument('--summary', required=True, type=Path, help='mosdepth summary file (genome.mosdepth.summary.txt)')
    p.add_argument('--globaldist', required=True, type=Path, help='mosdepth global.dist file (genome.mosdepth.global.dist.txt)')
    p.add_argument('--out-prefix', default='mosdepth_plots', type=Path, help='输出文件前缀（不含扩展名）')
    p.add_argument('--xlim', default=200, type=float, help='右侧图像 x 轴上限（深度），设为 0 则自动')
    p.add_argument('--sort-by', choices=['mean', 'chrom', 'none'], default='mean', help='左图染色体排序方式')
    p.add_argument('--dpi', default=100, type=int, help='输出 PNG 的 DPI')
    p.add_argument('--no-show', action='store_true', help='不在屏幕显示图像，仅保存')
    return p.parse_args()


def main():
    args = parse_args()

    summary_fp = args.summary
    dist_fp = args.globaldist
    out_prefix = args.out_prefix

    if not summary_fp.exists():
        sys.exit(f'错误：summary 文件不存在: {summary_fp}')
    if not dist_fp.exists():
        sys.exit(f'错误：global.dist 文件不存在: {dist_fp}')

    # 读取文件
    summary = read_summary(summary_fp)

    # 如果 xlim 为 0，表示自动根据平均深度设置
    if args.xlim == 0:
        avg_depth = summary['mean'].mean()
        xlim = 200 if pd.isna(avg_depth) or avg_depth < 100 else avg_depth * 2
    else:
        xlim = args.xlim

    dsel = read_global_dist(dist_fp)
    depths, coverage_vals, detected = detect_and_compute_coverage(dsel)

    # debug 信息输出
    print(f"[Info] coverage detection: {detected}")
    print(f"[Info] saving to {out_prefix}.png and {out_prefix}.pdf (dpi={args.dpi})")

    plot_mosdepth(
        left_df=summary,
        depths=depths,
        coverage_vals=coverage_vals,
        out_prefix=out_prefix,
        xlim=xlim,
        sort_by=args.sort_by,
        figsize=(16, 6),
        dpi=args.dpi,
        show=not args.no_show
    )


if __name__ == '__main__':
    main()
