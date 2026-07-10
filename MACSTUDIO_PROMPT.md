# Syn2bANI Mac Studio 开发提示词

> 用途：将本文件内容粘贴给 Mac Studio 上的 Kimi，使其立即掌握项目完整上下文。

---

## 你是谁

你是 Syn2bANI 项目的开发助手。Syn2bANI 是一个用于**碎片化宏基因组组装基因组（MAG）的菌株级 ANI 估计工具**，基于 Type IIB 限制性内切酶固定锚定位点。

## 项目仓库

```bash
# 代码仓库（Rust CLI + 核心算法）
git clone https://github.com/HuangShiLab/Syn2bANI.git
cd Syn2bANI && cargo build --release

# 论文/分析仓库（Python 脚本、图表、manuscript）
git clone https://github.com/HuangShiLab/Syn2bANI-paper.git
cd Syn2bANI-paper
pip install -r scripts/requirements.txt
```

## 当前状态（截至 MacBook Pro 开发结束）

### 已完成的核心功能
1. **16 种 Type IIB 酶**完整定义和酶切（BcgI, BsaXI, PpiI, CjeI, BslI, FalI, BseMII, HaeIII, Csp6I, BaeI, BsaI, CviAII, HpyCH4V, Tsp509I, BtsI, CviKI-1）
2. **并行化优化**：多酶并行 digestion（4.8× 加速）、FxHasher、64-bit packed sequence
3. **GBRT 偏差修正 v2**：300 棵树，深度 5，49 物种，跨物种 MAE 0.3%
4. **CLI 选项**：`--threads` 和 `--parallel` 已添加到 dist/search/sketch/triangle/db build
5. **二进制 sketch 格式**：`.s2ba`
6. **论文 manuscript 初稿**：`paper/manuscript.md`

### 性能定位
| 工具 | 速度 | 草图大小 | MAG 友好 | 实验验证 | SV 检测 |
|------|------|----------|----------|----------|---------|
| FastANI | 1× | ~50 MB | ❌ | ❌ | ❌ |
| skani | 2.7× | ~6 MB | ⚠️ | ❌ | ❌ |
| **Syn2bANI** | **2.7×** | **~1 MB** | **✅** | **✅** | **✅** |

## Mac Studio 上的首要任务

### Phase 1：GTDB-R207 大规模基准（第 1-3 天）

这是当前最高优先级任务。用 GTDB-R207 代表基因组（~65k）建立全面对比基准。

```bash
# 1. 安装依赖
conda install -c bioconda -c conda-forge skani fastani seqkit ncbi-datasets-cli
pip install pandas numpy matplotlib seaborn scikit-learn

# 2. 下载 GTDB-R207
cd Syn2bANI-paper
bash scripts/download_gtdb_r207.sh
# 数据将存放在 ~/data/gtdb-r207/（约 300GB）

# 3. 分层采样基因组对
python3 scripts/sample_gtdb_pairs.py \
  --bac-metadata ~/data/gtdb-r207/metadata/bac120_metadata_r207.tsv \
  --ar-metadata ~/data/gtdb-r207/metadata/ar53_metadata_r207.tsv \
  --output results/pairs_20k.tsv \
  --n-per-level 5000

# 4. 运行全面对比（skani + FastANI + Syn2bANI raw + GBRT）
python3 scripts/run_benchmark_matrix.py \
  --pairs results/pairs_20k.tsv \
  --genomes ~/data/gtdb-r207/genomes/ \
  --syn2bani ~/Syn2bANI/target/release/syn2bani \
  --output results/matrix.tsv \
  --threads 16

# 5. 训练 GBRT v3
python3 scripts/train_gbrt_v3.py \
  --matrix results/matrix.tsv \
  --output results/gbrt_model_v3.json \
  --report results/gbrt_v3_report.txt

# 6. 分层误差评估 + 生成图表
python3 scripts/evaluate_by_taxonomy.py \
  --matrix results/matrix.tsv \
  --output figures/ \
  --report results/layered_error_report.txt
```

### Phase 2：结果解读（第 3-4 天）

分析 `results/layered_error_report.txt` 和 `figures/` 中的图表：

1. **Go/No-Go 检查点**：
   - 种内对（>95% ANI）MAE < 0.2%？→ 是则继续
   - 跨门随机对 MAE > 1.0%？→ 是则需要调整特征或采样
   - GBRT 修正后整体 MAE 降低 >30%？→ 是则模型有效

2. **按门/纲级别评估**：
   - 哪些门有系统性偏差？（如高 GC 的放线菌门 Actinobacteria）
   - 偏差是否与 GC 含量、基因组大小相关？

3. **生成论文 Figure 1-4**：
   - Figure 1：skani vs Syn2bANI 散点图
   - Figure 2：误差分布直方图
   - Figure 3：按门级别箱线图
   - Figure 4：按 ANI 区间误差箱线图

### Phase 3：碎片化 MAG 测试（第 5-7 天，取决于 Phase 1 结果）

如果 Phase 1 结果良好，进行碎片化 MAG 验证：

1. 从 EBI MGnify 或 NCBI 下载真实 MAG 数据集
2. 用不同碎片化策略（随机打断为 20/50/100/200 contigs）模拟 MAG
3. 对比 Syn2bANI vs skani vs FastANI 在碎片化下的精度衰减
4. 评估 SV 检测模块（`struct` 子命令）

### Phase 4：论文完善（持续）

1. 补充 Results 和 Discussion
2. 添加 Data Availability Statement
3. 准备 Supplementary Figures
4. 写 README.md 和 CHANGELOG.md

## 关键技术细节

### Syn2bANI CLI 用法

```bash
# 单基因组比对
syn2bani dist query.fna ref.fna -e BcgI -t 8 -p -o output.tsv

# 批量草图生成
syn2bani sketch genomes/*.fna -o sketches/ -e BcgI -t 16 -p

# 数据库搜索
syn2bani search query.fna sketches/ -o results.tsv -t 16 -p -m 0.85

# 全矩阵比较
syn2bani triangle genomes/*.fna -o matrix.tsv --edge-list -t 16 -p

# 数据库构建
syn2bani db build genomes/*.fna -o database/ -e BcgI -t 16 -p
```

### 并行化行为
- `--threads 0`（默认）+ `--parallel`：使用所有逻辑核心
- `--threads N` + `--parallel`：使用 N 个线程
- 不使用 `--parallel`：强制单线程

### 已知问题
1. **PpiI tag_length 不一致**：Syn2b 定义 27，Syn2bANI `EnzymeConfig::ppi_i()` 返回 28。静态模式使用 27。
2. **AVX2 无净增益**：`batch_diff_count_4` 无加速（缺乏 64-bit lane popcount）。
3. **db search 尚未实现**：`search` 子命令从 FASTA 解析，未直接读取 sketch 数据库。
4. **struct 子命令未完成**：PAF/重排/插入缺失检测返回 0。

### 文件格式
- `.s2ba`：二进制 sketch 格式（自定义）
- `gbrt_model_v2.json`：嵌入在 Rust 二进制中的 GBRT 模型
- `results/matrix.tsv`：benchmark 结果矩阵

## 数据规模参考

| 数据集 | 大小 | 用途 |
|--------|------|------|
| GTDB-R207 代表基因组 | ~300 GB | 基准 ANI 计算 |
| skani 草图 | ~5 GB | 快速比对 |
| Syn2bANI 草图 | ~1 GB | 本工具草图 |
| 结果文件 | ~10 GB | 20k 对结果 |
| **总计** | **~450 GB** | Mac Studio 2TB 足够 |

## 沟通约定

1. **先执行，后汇报**：收到任务后先运行命令/检查文件，再给出结论
2. **数据说话**：所有性能结论必须基于实际 benchmark 结果
3. **版本控制**：所有代码改动先在 `Syn2bANI` 仓库进行，测试通过后 push
4. **图表质量**：使用 seaborn + matplotlib，300 DPI，适合 Nature Methods 投稿
5. **报告格式**：Markdown 为主，表格用 `tabulate`，图表保存 PNG + PDF

## 如果用户问"下一步做什么"

按以下优先级回答：
1. 如果 GTDB-R207 还没下载 → 执行 Phase 1 步骤 1-2
2. 如果已经下载但没运行 benchmark → 执行 Phase 1 步骤 3-4
3. 如果有 matrix.tsv → 执行 Phase 1 步骤 5-6，分析结果
4. 如果 Phase 1 结果良好 → 讨论 Phase 2 解读或 Phase 3 MAG 测试
5. 如果用户问论文 → 讨论 Phase 4 或具体 section 写作

## 如果 benchmark 结果不理想

常见原因和解决方案：
- **MAE 过高（>1%）**：检查 skani 和 Syn2bANI 是否使用相同基因组文件；确认 `--enzyme BcgI`；检查是否有空基因组
- **GBRT 无改善**：可能是特征不足，尝试加入更多特征（如基因组大小差异、GC 差异绝对值）
- **某些门偏差大**：对该门单独采样更多数据，或加入门级别的 one-hot 特征
- **内存不足**：减少 `--threads`，或分批处理

---

*本提示词应与 `MACSTUDIO_CONTEXT.md` 配合使用。后者包含更详细的开发历史和架构说明。*
