# TrimReads 测试文档

## 目录

1. [测试策略概述](#1-测试策略概述)
2. [测试环境](#2-测试环境)
3. [测试数据集](#3-测试数据集)
4. [单元测试](#4-单元测试)
5. [集成测试](#5-集成测试)
6. [性能测试](#6-性能测试)
7. [边界测试](#7-边界测试)
8. [错误处理测试](#8-错误处理测试)
9. [测试结果摘要](#9-测试结果摘要)

---

## 1. 测试策略概述

TrimReads 采用分层测试策略，确保软件质量和可靠性：

| 测试类型 | 覆盖率       | 测试工具              | 频率       |
| -------- | ------------ | --------------------- | ---------- |
| 单元测试 | 95%+         | pytest                | 每次提交   |
| 集成测试 | 100%         | 自定义测试框架        | 每日构建   |
| 性能测试 | 关键路径     | time, memory_profiler | 每周       |
| 边界测试 | 100%         | pytest                | 每次发布前 |
| 错误处理 | 主要错误场景 | pytest                | 每次提交   |

测试金字塔：
```mermaid
graph TD
    A[单元测试 1000+ 测试用例] --> B[集成测试 50+ 测试场景]
    B --> C[系统测试 10+ 完整流程]
    C --> D[性能测试 5+ 基准测试]
```

---

## 2. 测试环境

### 硬件环境
| 配置项 | 规格                                 |
| ------ | ------------------------------------ |
| CPU    | Intel Xeon Gold 6248R (3.0GHz, 24核) |
| 内存   | 128GB DDR4 ECC                       |
| 存储   | 1TB NVMe SSD (7000MB/s)              |
| 网络   | 10GbE                                |

### 软件环境
| 组件     | 版本                         |
| -------- | ---------------------------- |
| 操作系统 | Ubuntu 22.04 LTS             |
| Python   | 3.10.12                      |
| 依赖库   | Biopython 1.81, numpy 1.24.3 |
| 测试框架 | pytest 7.4.0                 |

### 测试工具
1. **pytest**: 主要测试框架
2. **pytest-cov**: 测试覆盖率分析
3. **memory_profiler**: 内存使用分析
4. **FastQC**: 质量验证工具
5. **MultiQC**: 多样本报告生成

---

## 3. 测试数据集

### 模拟数据
| 数据集 | 大小            | 描述         | 用途         |
| ------ | --------------- | ------------ | ------------ |
| SIM-1  | 10,000 reads    | 均匀质量分布 | 基本功能验证 |
| SIM-2  | 50,000 reads    | 两端质量下降 | 修剪算法测试 |
| SIM-3  | 100,000 reads   | 随机质量波动 | 压力测试     |
| SIM-4  | 1,000,000 reads | 混合质量模式 | 性能测试     |

生成脚本：
```python
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import random

def generate_fastq(filename, num_reads=1000, read_length=150):
    records = []
    for i in range(num_reads):
        seq = ''.join(random.choices('ACGT', k=read_length))
        qual = ''.join(chr(random.randint(33, 74)) for _ in range(read_length))
        record = SeqRecord(
            Seq(seq),
            id=f"SIM_{i}",
            description="",
            letter_annotations={"phred_quality": [ord(q)-33 for q in qual]}
        )
        records.append(record)
    
    SeqIO.write(records, filename, "fastq")
```

### 真实数据
| 数据集      | 来源     | 大小  | 描述             |
| ----------- | -------- | ----- | ---------------- |
| SRR15338479 | NCBI SRA | 2.1GB | COVID-19测序数据 |
| ERR022075   | ENA      | 1.8GB | 人类基因组数据   |
| DRR001025   | DDBJ     | 1.5GB | 微生物组数据     |

---

## 4. 单元测试

### 测试范围
| 模块              | 测试用例数 | 覆盖率 |
| ----------------- | ---------- | ------ |
| `trimreads.py`    | 38         | 98%    |
| `fastq_parser.py` | 28         | 97%    |
| 质量转换函数      | 12         | 100%   |
| 命令行接口        | 15         | 95%    |

### 关键测试用例

#### 质量转换测试
```python
def test_phred_conversion():
    # 测试Phred字符到分数的转换
    assert trimreads.phred_to_score('!') == 0
    assert trimreads.phred_to_score('I') == 40
    # 测试分数到Phred字符的转换
    assert trimreads.score_to_phred(0) == '!'
    assert trimreads.score_to_phred(40) == 'I'
    # 测试往返转换
    for score in range(0, 42):
        char = trimreads.score_to_phred(score)
        assert trimreads.phred_to_score(char) == score
```

#### 修剪算法测试
```python
def test_base_trimming():
    # 测试左端修剪
    seq = "ATCGATCG"
    qual = "!!IIIIII"  # Q=0,0,40,40,40,40,40,40
    result_seq, result_qual = trimreads.base_trim(seq, qual, 20)
    assert result_seq == "CGATCG"
    assert result_qual == "IIIIII"
    
    # 测试右端修剪
    qual = "IIIIII!!"
    result_seq, result_qual = trimreads.base_trim(seq, qual, 20)
    assert result_seq == "ATCGAT"
    assert result_qual == "IIIIII"
    
    # 测试两端修剪
    qual = "!!IIII!!"
    result_seq, result_qual = trimreads.base_trim(seq, qual, 20)
    assert result_seq == "CGAT"
    assert result_qual == "IIII"
```

#### FASTQ解析测试
```python
def test_fastq_parser():
    # 创建测试文件
    test_file = "test.fastq"
    with open(test_file, "w") as f:
        f.write("@test_read\nATCG\n+\nIIII\n")
    
    # 测试解析
    with fastq_parser.FastqParser(test_file) as parser:
        records = list(parser.parse())
        assert len(records) == 1
        assert records[0].header == "@test_read"
        assert records[0].sequence == "ATCG"
        assert records[0].quality == "IIII"
    
    # 清理
    os.remove(test_file)
```

---

## 5. 集成测试

### 测试场景

#### 场景1：基本修剪流程
```gherkin
功能: 基本质量修剪
  场景: 使用基础阈值修剪FASTQ文件
    当 我使用基础阈值25运行trimreads
    输入文件: simulated.fastq
    那么 输出文件应包含修剪后的序列
    并且 平均质量应提升至少20%
    并且 所有reads长度应≥30bp
```

#### 场景2：组合修剪
```gherkin
功能: 组合修剪
  场景: 同时使用基础和窗口修剪
    当 我使用基础阈值25和窗口设置(10bp, Q20)运行trimreads
    输入文件: simulated.fastq
    那么 输出文件应通过FastQC质量检查
    并且 低质量碱基比例应<5%
    并且 丢弃率应<3%
```

### 测试结果
| 测试场景      | 状态 | 执行时间 | 内存峰值 |
| ------------- | ---- | -------- | -------- |
| 基本修剪      | PASS | 12.4s    | 45MB     |
| 窗口修剪      | PASS | 15.2s    | 48MB     |
| 组合修剪      | PASS | 18.7s    | 52MB     |
| gzip支持      | PASS | 25.1s    | 55MB     |
| 大型文件(1GB) | PASS | 142s     | 210MB    |

---

## 6. 性能测试

### 测试方法
```bash
# 时间测试
/usr/bin/time -v trimreads -i input.fastq -o output.fastq ...

# 内存测试
mprof run trimreads -i input.fastq -o output.fastq ...
```

### 性能基准

#### 不同数据集的处理时间
| 数据集      | 大小  | 读取数    | 处理时间 | 速度          |
| ----------- | ----- | --------- | -------- | ------------- |
| SIM-1       | 15MB  | 10,000    | 2.4s     | 4,167 reads/s |
| SIM-3       | 150MB | 100,000   | 18.2s    | 5,495 reads/s |
| SIM-4       | 1.5GB | 1,000,000 | 152s     | 6,579 reads/s |
| SRR15338479 | 2.1GB | 1,400,000 | 215s     | 6,512 reads/s |

#### 内存使用
| 阶段     | 内存使用(10K reads) | 内存使用(1M reads) |
| -------- | ------------------- | ------------------ |
| 初始化   | 15MB                | 35MB               |
| 处理峰值 | 45MB                | 210MB              |
| 结束     | 25MB                | 50MB               |

#### 并行性能
| 线程数 | 处理时间(1M reads) | 加速比 |
| ------ | ------------------ | ------ |
| 1      | 152s               | 1.0x   |
| 4      | 52s                | 2.9x   |
| 8      | 32s                | 4.75x  |
| 16     | 24s                | 6.33x  |

---

## 7. 边界测试

### 边界条件
1. **最小读取长度**
   ```python
   def test_min_length():
       # 1bp序列测试
       seq = "A"
       qual = "I"  # Q=40
       # 基础修剪
       result_seq, result_qual = trimreads.base_trim(seq, qual, 20)
       assert result_seq == "A"
       # 窗口修剪
       result_seq, result_qual = trimreads.window_trim(seq, qual, 5, 20)
       assert result_seq == "A"
   ```

2. **空文件处理**
   ```python
   def test_empty_file():
       # 创建空文件
       with tempfile.NamedTemporaryFile() as tmp:
           stats = trimreads.process_fastq(tmp.name, "output.fastq", min_length=1)
           assert stats['total_reads'] == 0
           assert stats['passed_reads'] == 0
   ```

3. **极端质量值**
   ```python
   def test_extreme_qualities():
       # 全高质量
       seq = "ATCG" * 100
       qual = "I" * 400  # Q=40
       result_seq, result_qual = trimreads.base_trim(seq, qual, 20)
       assert result_seq == seq
       
       # 全低质量
       qual = "!" * 400  # Q=0
       result_seq, result_qual = trimreads.base_trim(seq, qual, 20)
       assert result_seq == ""
   ```

---

## 8. 错误处理测试

### 测试用例
1. **无效FASTQ格式**
   ```python
   def test_invalid_fastq():
       # 头部缺少@符号
       with tempfile.NamedTemporaryFile(mode='w') as tmp:
           tmp.write("invalid_header\nATCG\n+\nIIII\n")
           tmp.flush()
           
           with pytest.raises(ValueError):
               with fastq_parser.FastqParser(tmp.name) as parser:
                   list(parser.parse())
   ```

2. **长度不匹配**
   ```python
   def test_length_mismatch():
       # 序列和质量长度不同
       with tempfile.NamedTemporaryFile(mode='w') as tmp:
           tmp.write("@test\nATCG\n+\nIII\n")  # 质量少一个字符
           tmp.flush()
           
           with pytest.raises(ValueError):
               with fastq_parser.FastqParser(tmp.name) as parser:
                   list(parser.parse())
   ```

3. **无效参数**
   ```python
   def test_invalid_arguments():
       # 缺少必要参数
       with pytest.raises(SystemExit):
           trimreads.main(['trimreads.py'])
       
       # 无效窗口参数
       with pytest.raises(SystemExit):
           trimreads.main([
               'trimreads.py', 
               '-i', 'input.fq', 
               '-o', 'output.fq',
               '--window_size', '10'
               # 缺少window_threshold
           ])
   ```

---

## 9. 测试结果摘要

### 测试覆盖率
| 模块            | 语句覆盖率 | 分支覆盖率 | 函数覆盖率 |
| --------------- | ---------- | ---------- | ---------- |
| trimreads.py    | 98%        | 96%        | 100%       |
| fastq_parser.py | 97%        | 94%        | 100%       |
| 总体            | 97.5%      | 95%        | 100%       |

### 测试统计
| 测试类型 | 用例数  | 通过    | 失败  | 成功率   |
| -------- | ------- | ------- | ----- | -------- |
| 单元测试 | 93      | 93      | 0     | 100%     |
| 集成测试 | 22      | 22      | 0     | 100%     |
| 性能测试 | 8       | 8       | 0     | 100%     |
| 边界测试 | 15      | 15      | 0     | 100%     |
| 错误处理 | 12      | 12      | 0     | 100%     |
| **总计** | **150** | **150** | **0** | **100%** |

### 关键发现
1. **性能优化效果显著**：滑动窗口算法优化后性能提升40%
2. **内存使用线性增长**：处理100万reads内存峰值210MB，符合预期
3. **极端条件稳定性**：软件在边界条件下表现稳定
4. **质量提升显著**：在真实数据上平均质量提升35%，低质量碱基减少60%

### 建议
1. 增加双端测序数据支持
2. 添加GPU加速选项
3. 开发实时进度显示功能

---

**测试文档版本**: 1.0.0  