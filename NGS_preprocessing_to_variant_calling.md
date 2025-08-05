#### 登录服务器

```bash
ssh shenyz@222.29.188.221
```

#### 上传文件

```bash
scp /Users/shenyz/Downloads/delly_v1.3.3_linux_x86_64bit shenyz@222.29.188.221:/public/home/shenyz
```

#### 设置执行权限

```bash
chmod +x delly_v1.3.3_linux_x86_64bit
```

#### 执行程序

```bash
./delly_v1.3.3_linux_x86_64bit
```

#### 对比参考基因组和 .bam 文件

```bash
./delly_v1.3.3_linux_x86_64bit call -g fixed_human_virus.fa -o delly_output.bcf zzq_hv.sorted.bam
```

#### 使用多线程（可能不支持）

```bash
OMP_NUM_THREADS=12 ./delly_v1.3.3_linux_x86_64bit call -g fixed_human_virus.fa -o delly_output.bcf zzq_hv.sorted.bam
```

#### 检查插入 INS

```bash
OMP_NUM_THREADS=12 ./delly_v1.3.3_linux_x86_64bit call -g fixed_human_virus.fa -o delly_output.bcf map.sorted.bam
```

#### 只检测特定变异类型（如 INS）

```bash
OMP_NUM_THREADS=4 ./delly_v1.3.3_linux_x86_64bit call -t INS -g fixed_human_virus.fa -o delly_DEL.bcf map.sorted.bam
```

#### 重新比对（去掉 -ps 参数）

```bash
bwa mem -t 12 -v 1 fixed_human_virus.fa ../zzq_1.clean.fq.gz ../zzq_2.clean.fq.gz > map.sam
```

#### 重新排序

```bash
samtools sort -@ 12 -o map.sorted.bam map.sam
```

#### 查看文件大小

```bash
ll -rth
```

#### 提取 BND

```bash
bcftools view -i 'SVTYPE="BND"' -Ov -o delly_BND.vcf delly_ALL.vcf
```

#### 精确筛选 chr17 的 BND

```bash
bcftools view -i 'CHROM=="chr17" && SVTYPE=="BND"' -Ov -o delly_BND_chr17.vcf delly_ALL.vcf
```

#### 对数据进行清洗

```bash
fastp -i zzq2_1.fq.gz -I zzq2_2.fq.gz \
      -o zzq2_1.trim.fq.gz -O zzq2_2.trim.fq.gz \
      --thread 12 --detect_adapter_for_pe
```

#### 对数据进行比对

```bash
bwa mem -t 12 fixed_human_virus.fa \
        zzq2_1.trim.fq.gz zzq2_2.trim.fq.gz \
  | samtools view -b - \
  | samtools sort -@ 4 -o zzq2.sorted.bam
```

#### 建立索引

```bash
samtools index zzq2.sorted.bam
```

#### 添加 Bioconda & Conda-forge 源，再重新创建环境

```bash
conda config --add channels defaults
conda config --add channels bioconda
conda config --add channels conda-forge
conda config --set channel_priority strict
```

#### 创建新环境并安装 bcftools

```bash
conda create -n bcf_env bcftools=1.17 -y
```

#### 激活新环境

```bash
conda activate bcf_env
```

#### 转换 bcf 为 vcf

```bash
bcftools view all_sv.bcf -Ov -o delly_ALL2.vcf
```

#### 提取 BND

```bash
bcftools view -i 'SVTYPE="BND"' -Ov -o delly_BND2.vcf delly_ALL2.vcf
```

#### 精确筛选 chr17

```bash
bcftools view -i 'CHROM=="chr17" && SVTYPE=="BND"' -Ov -o delly_BND2_chr17.vcf delly_ALL2.vcf
```

#### 搜索是否有病毒

```bash
grep MN772835 delly_ALL2.vcf
```

#### 检查环境的命令

```bash
conda environments:
```

#### 激活之前的环境

```bash
conda activate bcftools_env
```

#### 查看版本号

```bash
bcftools 1.17
```

#### bgzip 压缩 VCF 文件

```bash
bgzip delly_BND_chr17.vcf
bgzip delly_BND2_chr17.vcf
```

#### 建立 tabix 索引

```bash
tabix -p vcf delly_BND_chr17.vcf.gz
tabix -p vcf delly_BND2_chr17.vcf.gz
```

#### 重新运行合并命令

```bash
bcftools concat -a -O v -o merged_BND_chr17.vcf delly_BND_chr17.vcf.gz delly_BND2_chr17.vcf.gz
```

#### 合并 BAM（注意一定得有索引文件）

```bash
samtools merge -@ 12 merged.bam zzq2.sorted.bam map.sorted.bam
```

#### 对 BAM 进行索引

```bash
samtools index merged.bam
```

#### 跑 delly

```bash
./delly_v1.3.3_linux_x86_64bit call \
  -g reference.fa \
  -o merged.sv.bcf \
  merged.bam
```

#### 筛选 BND

```bash
grep "SVTYPE=BND" merged.sv.vcf > merged.BND.vcf
```

#### 切换环境

```bash
conda env list
conda activate base
fastp
```

#### 比对（nohup 后台运行）

```bash
nohup bash -c "bwa mem -t 12 fixed_human_virus.fa zzq2_1.trim.fq.gz zzq2_2.trim.fq.gz | samtools view -b - | samtools sort -@ 4 -o zzq2.sorted.bam" &
```

#### 查看虚拟环境

```bash
conda info --envs
```

或

```bash
conda env list
```

#### 进入虚拟环境

```bash
conda activate myenv
```

#### 查看已安装包

```bash
conda list
```

#### 跑 delly

```bash
./delly_v1.3.3_linux_x86_64bit call \
  -g fixed_human_virus.fa \
  -o merged.sv.bcf \
  zzq2.sorted.bam
```

#### 切换环境

```bash
conda activate bcftools_env
```

#### 提取所有 BND 类型变异

```bash
bcftools view -i 'INFO/SVTYPE="BND"' merged.sv.vcf -Ov -o merged.BND.vcf
```

#### 查找 MN772835 染色体相关变异

```bash
grep MN772835 merged.sv.vcf > MN772835_related_variants.vcf
```

#### 重新对zzq_1和zzq_2进行清洗（7.24）

```bash
fastp -i zzq_1.fq.gz -I zzq_2.fq.gz \
      -o zzq_1.trim.fq.gz -O zzq_2.trim.fq.gz \
      --thread 12 --detect_adapter_for_pe

# 还可以加的参数（常用）：
      --length_required 30          # 设置最小 read 长度（过滤太短的）
      --qualified_quality_phred 20  # 设定 Q 值阈值（默认 Q20）
      --unqualified_percent_limit 40  # 超过这个比例的碱基低于质量阈值就丢弃
      --n_base_limit 5              # N 碱基不能超过几个

```

#### 对zzq_1和zzq_2进行比对和排序

```bash
bwa mem -t 12 fixed_human_virus.fa \
        zzq_1.trim.fq.gz zzq_2.trim.fq.gz \
  | samtools view -b - \
  | samtools sort -@ 4 -o zzq1.sorted.bam
```

#### 从 BAM 中提取目标 reads（包括 target 区域 + unmapped）
#### 提取比对到 chr17 和 MN772835.1 的 reads：

```bash
samtools view -b zzq2.sorted.bam chr17 MN772835.1 > target.chr.bam #注意病毒的名字要打全
```

#### 提取 unmapped 的 reads：

```bash
samtools view -b -f 4 zzq2.sorted.bam > unmapped.bam
```

#### 从这两个 BAM 文件中提取 read ID

```bash
samtools view target.chr.bam | awk '{print $1}' > target_reads.txt
samtools view unmapped.bam | awk '{print $1}' > unmapped_reads.txt
```

#### 合并并去重：

```bash
cat target_reads.txt unmapped_reads.txt | sort | uniq > final_reads.txt
```

#### 用 seqkit grep 从原始 fastq.gz 中提取这些 reads（失败）

```bash
seqkit grep -f final_reads.txt zzq_1.trim.fq.gz > subset_1.fq
seqkit grep -f final_reads.txt zzq_2.trim.fq.gz > subset_2.fq
```

#### 可能之前的samtools sort失败了，先只跑前半部分并保存 BAM 到中间文件

```bash
bwa mem -t 12 fixed_human_virus.fa zzq_1.trim.fq.gz zzq_2.trim.fq.gz \
  | samtools view -b - > zzq1.raw.bam
```
#### 然后再排序：

```bash
samtools sort -@ 4 -o zzq1.sorted.bam zzq1.raw.bam
···

## 匹配 ID 部分

```bash
seqkit grep -n -f final_reads.txt zzq_2.trim.fq.gz > subset_2.fq
```
#### 输入了新的命令，去掉了-n，因为名字不全
#### 输入`seqkit grep -h`可以查看命令

```bash
seqkit grep -f final_reads.txt zzq_1.trim.fq.gz > subset_1
```

#### 提取目标区域比对结果：

```bash
samtools view -b zzq1.sorted.bam chr17 MN772835.1 > target1.chr.bam
```
#### 提取未比对 reads

```bash
samtools view -b -f 4 zzq1.sorted.bam > unmapped1.bam
```

#### 提取 read ID：
```bash
samtools view target1.chr.bam | awk '{print $1}' > target1_reads.txt
samtools view unmapped1.bam | awk '{print $1}' > unmapped1_reads.txt
```


#### 合并并去重 read ID：
```bash
cat target1_reads.txt unmapped1_reads.txt | sort | uniq > final1_reads.txt
```

#### 用 seqkit grep 从原始 fastq.gz 中提取这些 reads

```bash
seqkit grep -f final1_reads.txt zzq_1.trim.fq.gz > subset1_1
seqkit grep -f final1_reads.txt zzq_2.trim.fq.gz > subset1_2
seqkit grep -f final_reads.txt zzq2_1.trim.fq.gz > subset2_1
seqkit grep -f final_reads.txt zzq2_2.trim.fq.gz > subset2_2
```

#### 用来查看一下压缩文件

```bash
gunzip -c zzq2_2.trim.fq.gz |head
```

```bash
head final1_reads.txt
tail final1_reads.txt
```

#### 查看文件生成的命令的历史

```bash
history | grep subset
```

#### 合并子集 FASTQ（R1 和 R2 各自合并）

```bash
cat subset1_1 subset2_1 > combined_R1.fq
cat subset1_2 subset2_2 > combined_R2.fq
```

#### 对fa文件restart一下，因为断点在334到1788，所以选1000

```bash
seqkit restart -i 1000 fixed_human_virus.fa > fixed_human_virus_restart1000.
```

#### Bwa建索引

```bash
bwa index fixed_human_virus_restart1000.fa
```

#### 比对 paired-end reads（R1 / R2）

```bash
bwa mem -t 12 -v 1 fixed_human_virus_restart1000.fa \
  combined_R1.fq combined_R2.fq > map_restart1000.sam
```

#### SAM 转 BAM
```bash
samtools view -@ 12 -Sb map_restart1000.sam > map_restart1000.bam
```

#### 排序
```bash
samtools sort -@ 12 -o map_restart1000.sorted.bam map_restart1000.bam
```

#### 建立索引
```bash
samtools index map_restart1000.sorted.bam
```

#### 跑delly

```bash
./delly_v1.3.3_linux_x86_64bit call \
  -g fixed_human_virus_restart1000.fa \
  -o merged_new.sv.bcf \
  map_restart1000.sorted.bam
```

#### 提取所有 BND 类型变异

```bash
bcftools view -i 'INFO/SVTYPE="BND"' merged_new.sv.vcf -Ov -o merged_new.BND.vcf
```

#### 查找 MN772835 染色体相关变异(失败)

```bash
grep MN772835 merged_new.sv.vcf > MN772835_related_variants_new.vcf
```

#### 对病毒restart一下

```bash
seqkit restart -i 1000 TTMV.fa > TTMV_restart1000.
```

#### 合并到人类基因组参考

```bash
cat GRCh38_fixed.fa TTMV_restart1000.fa > GRCh38_TTMV_restart1000.fa
```

#### 用 `bwa index` 对新参考建索引

```bash
bwa index GRCh38_TTMV_restart1000.fa
```

#### 比对双端 reads 到新的合并参考

```bash
bwa mem -t 12 GRCh38_TTMV_restart1000.fa combined_R1.fq combined_R2.fq > map_human_virus.sam
```

#### 转换为 BAM、排序、索引

```bash
samtools view -Sb map_human_virus.sam > map_human_virus.bam
samtools sort -o map_human_virus.sorted.bam map_human_virus.bam
samtools index map_human_virus.sorted.bam
```

#### 跑delly

```bash
./delly_v1.3.3_linux_x86_64bit call \
  -g GRCh38_TTMV_restart1000.fa \
  -o merged_human_virus.sv.bcf \
  map_human_virus.sorted.bam
```

#### 将 BCF 转换为 VCF

```bash
bcftools view merged_human_virus.sv.bcf > merged_human_virus.sv.vcf
```

#### 查找病毒

```bash
grep MN772835.1 merged_human_virus.sv.vcf > MN_related_variants.vcf
```

#### 一次性下载完

```bash
scp shenyz@222.29.188.221:/public/home/shenyz/genome/map_human_virus.sorted.bam ~/Downloads/for_igv_new/
scp shenyz@222.29.188.221:/public/home/shenyz/genome/map_human_virus.sorted.bam.bai ~/Downloads/for_igv_new/
scp shenyz@222.29.188.221:/public/home/shenyz/genome/GRCh38_TTMV_restart1000.fa ~/Downloads/for_igv_new/
scp shenyz@222.29.188.221:/public/home/shenyz/genome/GRCh38_TTMV_restart1000.fa.fai ~/Downloads/for_igv_new/
scp shenyz@222.29.188.221:/public/home/shenyz/genome/merged_human_virus.sv.vcf ~/Downloads/for_igv_new/
scp shenyz@222.29.188.221:/public/home/shenyz/genome/MN_related_variants.vcf ~/Downloads/for_igv_new/
```

#### 检查结构变异(失败)

```bash

./delly_v1.3.3_linux_x86_64bit call \
  -g GRCh38_TTMV_restart1000.fa \
  -o MN772835_DEL.bcf \
  -x human.hg38.excl.tsv \
  -r MN772835.1 \
  -t DEL \
  map_human_virus.sorted.bam

```

#### 提取 MN772835.1 的 BAM 区段

```bash
samtools view -h -b map_human_virus.sorted.bam MN772835.1 > map_MN772835.1.bam
samtools index map_MN772835.1.bam
```

#### 用 `delly` 检测结构变异

``` bash
./delly_v1.3.3_linux_x86_64bit call \
  -g GRCh38_TTMV_restart1000.fa \
  -o MN772835_sv.bcf \
  map_human_virus.sorted.bam
```

#### 转为vcf格式

```bash
bcftools view -r MN772835.1:1400-1700 MN772835_sv.vcf > MN_region.vcf
```

#### 用 bgzip 压缩 VCF 文件

```
bgzip MN772835_sv.vcf
```

这会生成一个压缩文件：`MN772835_sv.vcf.gz`

#### 用 tabix 建立索引

```
tabix -p vcf MN772835_sv.vcf.gz
```

这会生成索引文件：`MN772835_sv.vcf.gz.tbi`

#### 用 bcftools 提取指定区域

```
bcftools view -r MN772835.1:1400-1700 MN772835_sv.vcf.gz > MN_region.vcf
```





7.29

#### 建立新参考基因组索引

```bash
bwa index GRCh38_TTMV_inserted.fa
samtools faidx GRCh38_TTMV_inserted.fa
```



#### 全量数据重新比对

```bash
# 比对（使用原始trimmed数据）
bwa mem -t 32 GRCh38_TTMV_inserted.fa \
    zzq_1.trim.fq.gz zzq_2.trim.fq.gz > realigned.sam

# 转换排序索引
samtools view -@ 32 -Sb realigned.sam | samtools sort -@ 32 -o realigned.sorted.bam
samtools index realigned.sorted.bam
```



#### 挂在后台运行

```bash
nohup bash -c '
bwa index GRCh38_TTMV_inserted.fa \
&& samtools faidx GRCh38_TTMV_inserted.fa \
&& bwa mem -t 32 GRCh38_TTMV_inserted.fa zzq_1.trim.fq.gz zzq_2.trim.fq.gz > realigned.sam \
&& samtools view -@ 32 -Sb realigned.sam | samtools sort -@ 32 -o realigned.sorted.bam \
&& samtools index realigned.sorted.bam
' > pipeline.log 2>&1 &
```



#### 下载数据

```bash
scp shenyz@222.29.188.221:/public/home/shenyz/genome/GRCh38_TTMV_inserted.fa .
scp shenyz@222.29.188.221:/public/home/shenyz/genome/GRCh38_TTMV_inserted.fa.fai .
scp shenyz@222.29.188.221:/public/home/shenyz/genome/realigned.sorted.bam .
scp shenyz@222.29.188.221:/public/home/shenyz/genome/realigned.sorted.bam.bai .
```

#### 比对运行

```bash
nohup bash -c '
bwa index modified_genome.fa \
&& samtools faidx modified_genome.fa \
&& bwa mem -t 32 modified_genome.fa zzq_1.trim.fq.gz zzq_2.trim.fq.gz > realigned_new.sam \
&& samtools view -@ 32 -Sb realigned_new.sam | samtools sort -@ 32 -o realigned_new.sorted.bam \
&& samtools index realigned.sorted_new.bam
' > pipeline.log 2>&1 &
```

#### new比对

```bash
nohup bash -c '
bwa index modified_genome.fa \
&& samtools faidx modified_genome.fa \
&& bwa mem -t 12 modified_genome.fa zzq_1.trim.fq.gz zzq_2.trim.fq.gz > realigned_new.sam \
&& samtools view -@ 12 -Sb realigned_new.sam | samtools sort -@ 12 -o realigned_new.sorted.bam \
&& samtools index realigned_new.sorted.bam
' > pipeline.log 2>&1 &

```

#### 观察基因是不是对的上，找到insertion的原因

```bash
samtools faidx GRCh38.primary_assembly.genome.fa chr17:40332575-40332590

samtools faidx TTMV_restart1000.fa MN772835.1:1910-1930
```

#### 创建脚本文件

``` bash
nano pipeline_realignment.sh
```

#### 写入脚本

```bash
#!/usr/bin/env bash
set -euo pipefail

log() {
  echo -e "\033[1;32m[$(date +'%F %T')] $*\033[0m"
}

REF="modified_genome.fa"
READ1="zzq_1.trim.fq.gz"
READ2="zzq_2.trim.fq.gz"
OUT_SAM="realigned_new.sam"
OUT_BAM="realigned_new.sorted.bam"
THREADS=12

log "① 构建 BWA 索引..."
bwa index "$REF"

log "② 构建 faidx 索引..."
samtools faidx "$REF"

log "③ 比对生成 SAM 文件..."
bwa mem -t "$THREADS" "$REF" "$READ1" "$READ2" > "$OUT_SAM"

log "④ SAM 转 BAM + 排序..."
samtools view -@ "$THREADS" -Sb "$OUT_SAM" | samtools sort -@ "$THREADS" -o "$OUT_BAM" -

log "⑤ 构建 BAM 索引..."
samtools index "$OUT_BAM"

log "✅ 全部完成！输出文件: $OUT_BAM"
```

#### 保存并退出编辑器

在 `nano` 中：

- 按下 `Ctrl + O`（写出/保存）
- 然后按回车（确认保存）
- 然后按 `Ctrl + X`（退出）
