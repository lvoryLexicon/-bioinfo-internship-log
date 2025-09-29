# mosdepth_plot

一个用于绘制 mosdepth 输出的 summary 和全局覆盖度分布的命令行小工具。

## 文件说明

- `mosdepth_plot.py` — 主脚本。读取 `*.summary.txt`（每条染色体的深度汇总）和 `*.global.dist.txt`（覆盖度分布），生成两张图：左图显示每条染色体的平均深度；右图显示深度 ≥ 某值的百分比曲线。输出 `OUTPREFIX.png` 和 `OUTPREFIX.pdf`。
- `README.md` — 本文件。
- `requirements.txt` — Python 依赖。

------

## 环境依赖

- Python 3.8+
- 依赖库见 `requirements.txt`，安装命令：

```bash
python -m pip install -r requirements.txt
```

------

## 快速使用示例

```bash
python mosdepth_plot.py \
  --summary genome.mosdepth.summary.txt \
  --globaldist genome.mosdepth.global.dist.txt \
  --out-prefix mosdepth_plots \
  --xlim 200 \
  --sort-by mean \
  --dpi 150
```

参数说明

- `--summary`（必填）：mosdepth summary 文件，如 `genome.mosdepth.summary.txt`。
- `--globaldist`（必填）：mosdepth global.dist 文件，如 `genome.mosdepth.global.dist.txt`。
- `--out-prefix`：输出文件前缀（默认: `mosdepth_plots`），会生成 `.png` 和 `.pdf`。
- `--xlim`：右图 x 轴上限（深度），设置为 `0` 表示自动根据平均深度设置，默认 `200`。
- `--sort-by`：左图染色体排序方式，选项：`mean`（默认，按平均深度降序）、`chrom`（自然染色体顺序 1..22,X,Y）、`none`（保持原顺序）。
- `--dpi`：PNG 输出分辨率，默认 `100`。
- `--no-show`：不在屏幕显示图像，仅保存文件。

------

## 功能说明

- 脚本支持常见的 mosdepth 输出格式（2 列或 3 列 global.dist）。如果 `group` 列中包含 `total`，优先使用该组；否则使用第一组或两列数据。
- `percent` 值可接受 0..1 或 0..100 的范围，脚本会自动归一化。
- 脚本会检测 global.dist 是按每深度频率还是累积覆盖率，并转换为 "深度 ≥ x 的百分比" 曲线。
- 左图默认筛选染色体 `1..22` 和 `X/Y`，支持从 `chr1`、`1`、`ChrX` 等格式自动提取染色体 ID。

------

## 示例

- 保存图像，不弹出窗口：

```bash
python mosdepth_plot.py --summary genome.mosdepth.summary.txt --globaldist genome.mosdepth.global.dist.txt --no-show
```

- 使用自然染色体顺序：

```bash
python mosdepth_plot.py --summary genome.mosdepth.summary.txt --globaldist genome.mosdepth.global.dist.txt --sort-by chrom
```

------

## 故障排查

- 如果提示文件读取错误，请检查文件路径是否正确以及文件是否为空。
- 如果覆盖度曲线看起来异常，请不要使用 `--no-show`，检查打印的检测信息，脚本会输出类似 `per-depth -> computed coverage as ...` 的提示帮助诊断。

------

## 开发 / 集成

- 如果希望在 notebook 或其他 Python 代码中调用绘图函数，可以直接导入 `mosdepth_plot.py` 中的函数：

  ```python
  from mosdepth_plot import plot_mosdepth, read_summary, read_global_dist, detect_and_compute_coverage
  ```



------

## requirements.txt（推荐内容）

```
pandas>=1.5
matplotlib>=3.5
numpy>=1.24
```