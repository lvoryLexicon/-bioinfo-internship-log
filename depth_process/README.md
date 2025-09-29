#### depth_process.py

简洁说明：从 bigWig (depth.bw) 读取覆盖度并绘制每条染色体的水平 track（PNG/PDF）。

#### 要求
- Python 3.8+

#### 安装（建议）
```bash
python3 -m venv env
source env/bin/activate      # Windows: env\Scripts\activate
pip install --upgrade pip
pip install -r requirements.txt
```

#### 运行示例

```
python depth_process.py \
  --bw /path/to/depth.bw \
  --out /path/to/coverage_tracks.png \
  --bins 1000 --log 
```
