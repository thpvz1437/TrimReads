# TrimReads 用户手册

## 目录

1. [简介](#1-简介)  
2. [安装指南](#2-安装指南)  
3. [快速开始](#3-快速开始)  
4. [详细使用说明](#4-详细使用说明)  
5. [输出说明](#5-输出说明)  
6. [使用示例](#6-使用示例)  
7. [高级功能](#7-高级功能)  
8. [常见问题](#8-常见问题)  
9. [技术支持](#9-技术支持)  

---

## 1. 简介

TrimReads 是一个用于高通量测序数据质量修剪的 Python 工具包。它提供两种修剪算法：
- **逐个碱基修剪**：从 reads 两端移除低质量碱基
- **滑动窗口修剪**：基于窗口平均质量进行修剪

主要特点：
- 支持标准 FASTQ 格式（包括 gzip 压缩）
- 详细的统计输出
- 灵活的修剪参数配置
- 高效处理大型文件

## 2. 安装指南

### 系统要求
- Python 3.8+
- Linux/macOS/Windows (建议使用 Linux 环境)

### 安装方法

#### 通过 PyPI 安装
```bash
pip install TrimReads
```

#### 从源代码安装
```bash
git clone https://github.com/yourusername/TrimReads.git
cd TrimReads
pip install .
```

#### 验证安装
```bash
trimreads --version
fastq_utils --help
```

## 3. 快速开始

### 基本修剪
```bash
# 逐个碱基修剪 (Q≥25)
trimreads -i input.fastq -o trimmed.fastq --base_threshold 25

# 窗口修剪 (10bp窗口，平均Q≥20)
trimreads -i input.fastq -o trimmed.fastq --window_size 10 --window_threshold 20

# 组合修剪
trimreads -i input.fastq -o trimmed.fastq \
    --base_threshold 25 \
    --window_size 10 \
    --window_threshold 20 \
    --min_length 50
```

### FASTQ 工具
```bash
# 验证 FASTQ 文件
fastq_utils validate input.fastq

# 获取文件统计
fastq_utils stats input.fastq

# 过滤低质量 reads
fastq_utils filter input.fastq output.fastq --min_quality 25 --min_length 50
```

## 4. 详细使用说明

### trimreads 参数

| 参数                 | 缩写 | 默认值 | 描述                           |
| -------------------- | ---- | ------ | ------------------------------ |
| `--input`            | `-i` | 无     | 输入 FASTQ 文件路径 (支持 .gz) |
| `--output`           | `-o` | 无     | 输出 FASTQ 文件路径 (支持 .gz) |
| `--base_threshold`   | 无   | 无     | 碱基质量阈值 (Q值)             |
| `--window_size`      | 无   | 无     | 滑动窗口大小 (碱基数)          |
| `--window_threshold` | 无   | 无     | 窗口平均质量阈值 (Q值)         |
| `--min_length`       | 无   | 30     | 修剪后最小保留长度             |
| `--help`             | `-h` | 无     | 显示帮助信息                   |

### fastq_utils 子命令

#### validate
验证 FASTQ 文件格式是否正确
```bash
fastq_utils validate input.fastq
```

#### stats
获取 FASTQ 文件统计信息
```bash
fastq_utils stats input.fastq
```

#### filter
基于质量和长度过滤 reads
```bash
fastq_utils filter input.fastq output.fastq \
    --min_quality 20 \
    --min_length 50 \
    --max_length 200
```

| 参数            | 描述                |
| --------------- | ------------------- |
| `--min_quality` | 最小平均质量 (Q值)  |
| `--min_length`  | 最小序列长度        |
| `--max_length`  | 最大序列长度 (可选) |

## 5. 输出说明

### trimreads 输出

1. **控制台输出**：
   ```
   Starting quality trimming...
   
   === Trimming Summary ===
   Total reads processed: 1000000
   Reads passing filters: 956782 (95.68%)
   Reads discarded: 43218 (4.32%)
   Base-by-base threshold: Q25
   Window-based trimming: 10bp window, Q25 threshold
   Minimum length kept: 50bp
   Output saved to: trimmed.fastq
   ```

2. **输出文件**：
   - 修剪后的 FASTQ 文件
   - 保留原始格式和注释

### fastq_utils 输出

#### stats 输出
```
FASTQ File Statistics for input.fastq:
Total records: 1,000,000
Total bases: 150,000,000
Average sequence length: 150.0 bp
Min sequence length: 150 bp
Max sequence length: 150 bp
Average quality score: 28.45
Bases with Q<20: 12.34%
```

#### filter 输出
```
Filtering Statistics:
Total records: 1,000,000
Passed records: 850,000 (85.00%)
Low quality: 100,000
Too short: 50,000
Too long: 0
```

## 6. 使用示例

### 示例 1：基础修剪
```bash
trimreads -i raw_data.fastq.gz -o trimmed_data.fastq.gz \
    --base_threshold 20 \
    --min_length 40
```

### 示例 2：高级修剪
```bash
trimreads -i large_dataset.fastq -o processed.fastq \
    --base_threshold 25 \
    --window_size 15 \
    --window_threshold 20 \
    --min_length 60
```

### 示例 3：FASTQ 分析管道
```bash
# 步骤 1: 验证数据
fastq_utils validate raw_data.fastq

# 步骤 2: 获取统计信息
fastq_utils stats raw_data.fastq > raw_stats.txt

# 步骤 3: 质量修剪
trimreads -i raw_data.fastq -o trimmed.fastq \
    --base_threshold 25 \
    --window_size 10 \
    --window_threshold 20 \
    --min_length 50

# 步骤 4: 过滤短序列
fastq_utils filter trimmed.fastq final.fastq --min_length 70

# 步骤 5: 最终统计
fastq_utils stats final.fastq > final_stats.txt
```

## 7. 高级功能

### Python API 使用
```python
from trimreads import process_fastq

# 处理 FASTQ 文件
stats = process_fastq(
    "input.fastq",
    "output.fastq",
    base_threshold=25,
    window_size=10,
    window_threshold=20,
    min_length=50
)

print(f"保留的 reads: {stats['passed_reads']}")
```

### 自定义质量转换
```python
from trimreads import base_trim

# 自定义修剪函数
def custom_trim(sequence, quality):
    # 预处理逻辑
    processed_seq, processed_qual = base_trim(sequence, quality, threshold=25)
    # 后处理逻辑
    return processed_seq, processed_qual
```

### 流式处理大型文件
```python
from trimreads.fastq_parser import stream_fastq_records

def process_record(record):
    # 自定义处理逻辑
    return record

# 流式处理
for record in stream_fastq_records(
    "large_input.fastq.gz", 
    "processed_output.fastq.gz",
    process_func=process_record
):
    # 处理每个记录
    pass
```

## 8. 常见问题

### Q1：我应该使用哪种修剪方法？
- **逐个碱基修剪**：适用于两端有明显质量下降的数据
- **窗口修剪**：适用于有局部质量问题的数据
- **组合方法**：提供最全面的质量提升（推荐）

### Q2：如何选择阈值？
- Q20：可接受质量（错误率 1%）
- Q25：良好质量（错误率 0.3%）
- Q30：高质量（错误率 0.1%）

### Q3：处理大型文件时内存不足？
- 使用 gzip 压缩输入输出
- 确保有足够的磁盘空间
- 考虑分批处理

### Q4：如何验证结果？
- 使用 `fastq_utils stats` 比较修剪前后统计
- 使用 FastQC 等工具可视化质量

### Q5：支持双端测序数据吗？
当前版本需要分别处理 R1 和 R2 文件，未来版本将支持同步处理。

## 9. 技术支持

### 文档资源
- [在线文档](https://trimreads.readthedocs.io)
- [GitHub 仓库](https://github.com/yourusername/TrimReads)
- [示例数据集](https://github.com/yourusername/TrimReads/tree/main/data)

---

**版本**: 1.0.0  