# README.md

## Metadata Stats (Python)

该小工具读取 metadata（CSV/TSV），计算若干统计量，并生成黑底白字的可视化：左侧为文本统计信息，右侧为 `segmentation_method` 的饼图（仅显示前 `top_n` 项）。

------

### 文件

- `metadata_stats.py` — 包含主函数 `plot_metadata_stats(path, sep=None, top_n=10, save_path=None, show=True)`，同时支持命令行调用。

------

### 功能简介

- 读取 metadata 文件（CSV 或 TSV）
- 计算统计量：总细胞数、`nCount_Xenium` 与 `nFeature_Xenium` 的中位数与均值
- 统计 `segmentation_method` 的频次与百分比，输出前 `top_n` 项
- 绘制 2 列图：左为文本统计、右为饼图（黑底白字风格）
- 支持保存为图片文件（PNG/PDF）

------

### 环境要求

- Python 3.8+
- 安装依赖（见 `requirements.txt`）

建议使用 virtualenv 或 conda 环境：

```bash
python -m venv .venv
source .venv/bin/activate      # macOS / Linux
.venv\Scripts\activate        # Windows（PowerShell/CMD）
pip install -r requirements.txt
```

------

### 使用方法

#### 作为模块导入

```python
from metadata_stats import plot_metadata_stats

stats, seg_stats = plot_metadata_stats(
    "/path/to/metadata.csv",  # 支持相对/绝对路径
    sep=None,                  # None -> pandas 自动推断，若是 tsv 请传 "\t"
    top_n=10,
    save_path="out.png",     # 可选：保存输出
    show=True                 # 是否显示窗口
)
print(stats)
print(seg_stats.head())
```

#### 命令行使用

```bash
python metadata_stats.py --path /Users/you/metadata.csv --save metadata_plot.png
# 如果是 tsv：
python metadata_stats.py --path metadata.csv --sep '\t' --save metadata_plot.png
```

CLI 参数说明：

- `--path`：必需，metadata 文件路径（csv/tsv）。
- `--sep`：可选，分隔符（默认 `None`，使用 pandas 自动推断；若文件为 TSV 请传 `'\t'`）。
- `--top_n`：可选，显示 segmentation_method 的前 N 项，默认 10。
- `--save`：可选，保存绘图为该文件路径。

------

### 常见问题与排查

- **FileNotFoundError / 无法读取文件**：确认你在正确的工作目录或使用了正确的绝对路径（不要写成 `/metadata.csv` 除非文件确实位于根目录）。
- **分隔符问题**：如果数据是制表符分隔（TSV），传 `--sep '\t'`。如果第一行是表头但列名不对，打开文件检查列名是否为 `nCount_Xenium` / `nFeature_Xenium` / `segmentation_method`。
- **图像保存失败**：确认 `save_path` 的目录存在，且有写权限。

------

### 扩展建议（可选）

- 使用 `seaborn` 优化美观风格或改成条形图显示 `segmentation_method`，便于展示较小比例的类别。
- 输出更多统计量（quartiles、missing rate 等）。
  