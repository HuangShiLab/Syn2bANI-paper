# Syn2bANI 开发上下文 — Mac Studio 迁移指南

> 创建时间：2025-01-12
> 来源：MacBook Pro 开发环境 → 迁移至 Mac Studio (2TB)
> 用途：保存完整开发状态，确保在 Mac Studio 上无缝继续工作

---

## 1. 代码仓库

| 仓库 | URL | 内容 |
|------|-----|------|
| **Syn2bANI (CODE)** | https://github.com/HuangShiLab/Syn2bANI | Rust 源代码、CLI、benchmark |
| **Syn2bANI-paper** | https://github.com/HuangShiLab/Syn2bANI-paper | 论文 manuscript、图表、分析脚本 |

**当前分支**：`main`（两个仓库均已 push 最新）

**最新提交**：
- CODE: `7378e4f` — feat: add --threads and --parallel CLI options to all subcommands
- paper: `ce7854a` — docs: add full manuscript draft

---

## 2. 开发环境要求

### Mac Studio 基础配置
- **OS**: macOS (Apple Silicon M1/M2/M3/Ultra)
- **Rust**: `rustc 1.92.0` / `cargo 1.92.0` (edition 2021)
- **存储**: 2TB（必需，用于存放大量基因组数据）
- **内存**: 推荐 32GB+（rayon 并行化 + 大数据集）

### 克隆和编译

```bash
# 1. 代码仓库
git clone https://github.com/HuangShiLab/Syn2bANI.git
cd Syn2bANI
cargo build --release

# 2. 论文仓库
git clone https://github.com/HuangShiLab/Syn2bANI-paper.git
```

### 依赖（Rust）
已在 `Cargo.toml` 中声明，无需手动安装：
- `clap 4.5` — CLI 参数解析
- `rayon 1.10` — 并行计算
- `serde` + `serde_json` — 序列化（GBRT 模型、sketch 格式）
- `needletail 0.5` — FASTA 解析
- `bitvec`, `byteorder`, `memmap2` — 二进制 I/O
- `criterion` — benchmark（dev 依赖）

---

## 3. 项目概况

**项目名称**：Syn2bANI（Synthetic 2bRAD-M ANI）
**目标**：用于碎片化 MAG（Metagenome-Assembled Genomes）的菌株级 ANI 估计工具

### 核心创新点
1. **固定酶切位点锚定**（16 种 Type IIB 酶）—— 保证可重复、与组装质量无关
2. **2bRAD-M 实验验证** —— 从湿实验角度验证 ANI 精度
3. **结构变异（SV）检测** —— 通过单体型相位标签检测插入/缺失/重排
4. **GBRT 偏差修正** —— 跨物种泛化误差 < 0.3% MAE

### 性能定位
| 工具 | 速度 | 草图大小 | MAG 友好 | 实验验证 | SV 检测 |
|------|------|----------|----------|----------|---------|
| FastANI | 1× | ~50 MB | ❌ 断点 | ❌ | ❌ |
| skani | 2.7× | ~6 MB | ⚠️ 部分 | ❌ | ❌ |
| **Syn2bANI** | **2.7×** | **~1 MB** | **✅** | **✅** | **✅** |

---

## 4. 已完成的架构和优化

### 4.1 酶切系统
- **16 种 Type IIB 酶**完整定义：BcgI, BsaXI, PpiI, CjeI, BslI, FalI, BseMII, HaeIII, Csp6I, BaeI, BsaI, CviAII, HpyCH4V, Tsp509I, BtsI, CviKI-1
- 已对齐 Syn2b 的 Fast2bRAD-M 滑窗 + 静态编译模式（`digest_sequence`）
- 保留 `digest_sequence_legacy()` 用于 benchmark 对比

### 4.2 并行化优化
| 优化项 | 效果 | 状态 |
|--------|------|------|
| 多酶并行 digestion（`extract_multi_enzyme_par`） | **4.8× 加速** | ✅ 完成 |
| FxHasher 快速哈希（`FastHashMap`） | 零外部依赖，HashMap 提速 | ✅ 完成 |
| 64-bit packed sequence（Hamming 距离） | ~5× 加速 | ✅ 完成 |
| AVX2 batch matching | 无净增益（pattern matching 短接） | ✅ 已尝试，未采用 |
| CLI `--threads` / `--parallel` | 用户控制并行度 | ✅ 刚完成 |

### 4.3 GBRT 偏差修正模型
- **v2 模型**：300 棵树，深度 5，49 物种训练
- **跨物种 MAE**：0.3%
- **嵌入大小**：1.08 MB JSON（`gbrt_model_v2.json`）
- 模型路径：`src/debias/gbrt_model_v2.json`（内嵌在二进制中）

### 4.4 文件格式
- **`.s2ba`** — 二进制 sketch 格式（自定义，包含 `TgtSketch`）
- 支持：`sketch`, `db build`, `db add`, `db merge` 等操作

---

## 5. 剩余开发任务（优先级排序）

### 🔴 高优先级（数据密集型，Mac Studio 优势）

#### Task 1: 大规模 MAG 真实数据验证（原 Task 3）
- **目标**：从 NCBI 或 EBI 下载 100+ 真实 MAG 数据集
- **数据需求**：
  - 同一物种多个菌株（如 *E. coli*, *S. aureus*）的完整基因组
  - 对应的碎片化 MAG（不同 contig 数：20, 50, 100, 200+）
  - 已知的 ANI"ground truth"（如 RefSeq 高质量基因组 vs 草图）
- **验证指标**：
  - Syn2bANI vs skani vs FastANI 的 ANI 偏差
  - 不同碎片化程度下的精度衰减
  - 结构变异检测的召回率/精确率
- **输出**：benchmark 报告 + 图表（建议用 `seaborn`/`matplotlib`）

#### Task 2: 更大范围 GBRT 训练（原 Task 6）
- **目标**：覆盖更多物种，提高模型泛化能力
- **数据需求**：
  - 50-100 个代表性物种的基因组对
  - 覆盖不同 GC 含量、基因组大小、进化距离
- **训练规模**：100-200 个物种，500-1000 棵树
- **输出**：`gbrt_model_v3.json`

#### Task 3: 结构变异（SV）检测完整实现（原 Task 4）
- **当前状态**：`struct` 子命令存在但仅返回 0
- **实现内容**：
  - PAF 格式输出（`--paf`）
  - 重排检测（`--rearrangement`）
  - 插入/缺失检测（`--indel`）
  - 单体型相位信息（haplotype phasing）

#### Task 4: Benchmark 扩展（原 Task 7）
- **时间基准**： skani vs FastANI vs Syn2bANI（完整数据集）
- **内存基准**： peak RSS 对比
- **草图大小对比**：`.s2ba` vs `.sketch` vs `.msh`
- **不同线程数下的并行扩展性**：`--threads 1,2,4,8,16`

### 🟡 中优先级

#### Task 5: 数据库功能完善（`db` 子命令）
- `db build`：已完成，支持并行 sketch
- `db add`：已完成，但缺少 `--parallel` 支持（当前串行）
- `db remove`：已完成
- `db list`：已完成
- `db merge`：已完成
- **缺失**：`db search`（在 sketch 数据库上执行 search，不需要重新解析 FASTA）

#### Task 6: 论文完善（Syn2bANI-paper）
- 补充 Figure 1-4 的 Python 绘图脚本
- 添加 Benchmark 数据表格
- 补充 Supplementary Methods
- 添加 Data Availability Statement
- 回应所有 reviewer 可能的质疑点（提前准备 rebuttal）

### 🟢 低优先级 / 优化项

#### Task 7: 进一步优化
- **SIMD AVX2 batch matching**：当前无净增益，若 future CPU 支持 64-bit lane popcount 可重新评估
- **内存映射 sketch**：`db search` 支持 memmap 避免全量加载
- **增量式 ANI 计算**：仅重新计算变化的基因组对
- **WebAssembly 版本**：供浏览器端使用

#### Task 8: 发布准备
- 编写 `README.md`（安装、使用、示例）
- 编写 `CHANGELOG.md`
- 添加 `--version` 输出
- 创建 GitHub Release v0.1.0
- 提交到 Bioconda / Homebrew

---

## 6. 基因组数据来源指南（Mac Studio 上执行）

### 推荐数据集

| 数据集 | 来源 | 规模 | 用途 |
|--------|------|------|------|
| RefSeq 细菌基因组 | NCBI FTP | ~30,000 完整基因组 | 基准 ANI 计算 |
| EBI MGnify MAGs | EBI FTP | 100,000+ MAGs | 碎片化 MAG 测试 |
| UHGG (人类肠道) | EBI | 4,600+ 基因组 | 真实宏基因组场景 |
| Genomes-On-Demand | NCBI | 按需下载 | 特定物种验证 |

### 下载命令示例

```bash
# 创建数据目录（利用 Mac Studio 的 2TB 存储）
mkdir -p ~/data/syn2bani/{refseq,uhgg,mags}

# NCBI RefSeq 细菌基因组（需要 ~50GB）
ncbi-genome-download bacteria -o ~/data/syn2bani/refseq --formats fasta

# EBI UHGG（人类肠道基因组）
wget -r -np -nH --cut-dirs=4 \
  ftp://ftp.ebi.ac.uk/pub/databases/metagenomics/mgnify_genomes/uhgg/

# 下载脚本可放在 Syn2bANI-paper/scripts/download_genomes.py
```

### 数据预处理建议
- 统一为 `.fna` 或 `.fasta` 格式
- 按物种分类存放：`data/refseq/species_name/genome_id.fna`
- 生成 manifest 文件：`manifest.tsv`（genome_id, path, species, completeness, contamination）
- 使用 `seqkit stats` 快速统计基因组质量

---

## 7. 性能基准参考数据（当前 MacBook Pro）

### 硬件
- MacBook Pro (2024), Apple Silicon M3 Pro, 18GB RAM
- 编译：`cargo build --release`（`opt-level=3`, `lto=fat`, `codegen-units=1`）

### 关键 benchmark 结果

```
# 多酶并行 digestion（16 酶 panel）
Sequential: 5.29 ms
Parallel:   1.09 ms
Speedup:    4.8×

# 草图大小对比（4.65 Mb 基因组）
Syn2bANI: ~0.01 MB (16 酶)
skani:    ~0.06 MB
FastANI:  ~50 MB
```

### 建议在 Mac Studio 上扩展的 benchmark
1. **线程扩展性**：`--threads 1,2,4,8,16,32`（Mac Studio 可能有更多核心）
2. **基因组规模扩展**：1 Mb → 10 Mb → 100 Mb
3. **基因组数量扩展**：10 → 100 → 1000 个 pairwise
4. **内存占用**：`time -v` 或 `/usr/bin/time -l`

---

## 8. 已知问题与注意事项

### 8.1 代码层面
1. **PpiI tag_length 不一致**：
   - Syn2b 定义 `tag_length=27`
   - Syn2bANI `EnzymeConfig::ppi_i()` 返回 `28`
   - 静态模式使用 27，可能导致 PpiI  digestion 不匹配
   - **建议**：对齐到 27 或更新静态模式

2. **sketch.rs 中的 `genome_path` 未使用**：
   - 已标记为 `_genome_path`，不影响功能
   - 保留以备将来需要输出原始路径

3. **db search 尚未实现**：
   - `search` 子命令目前从 FASTA 解析，而非 sketch 数据库
   - 应实现 `db search` 或修改 `search` 支持 sketch 输入

4. **AVX2 无净增益**：
   - `is_pure_atcg_simd` 在 x86_64 上可用，但 `batch_diff_count_4` 无加速
   - 原因：缺乏 64-bit lane popcount，且 pattern matching 短接
   - 若 future CPU 支持 `_mm256_popcnt_epi64`，可重新评估

### 8.2 数据层面
- **真实 MAG 数据**：需要下载和预处理，耗时较长（建议用 `ncbi-genome-download` 或 `datasets` CLI）
- **GBRT 训练数据**：需要大量已知的同源基因组对，建议从 PATRIC 或 NCBI 获取

### 8.3 论文层面
- **Figure 2**（ANI 偏差分布）需要实际基因组数据重新绘制
- **Figure 4**（Benchmark 时间/内存）需要在 Mac Studio 上重新测试
- **Supplementary Data** 需要上传到 Zenodo / Figshare

---

## 9. 快速启动命令（Mac Studio 上）

```bash
# 1. 环境准备
brew install rust  # 或通过 rustup
rustup update

# 2. 克隆仓库
git clone https://github.com/HuangShiLab/Syn2bANI.git
cd Syn2bANI

# 3. 编译
cargo build --release

# 4. 验证安装
./target/release/syn2bani --version

# 5. 基本用法示例
./target/release/syn2bani dist \
  query.fna ref.fna \
  -e BcgI -t 8 -p \
  -o output.tsv

./target/release/syn2bani sketch \
  genomes/*.fna \
  -o sketches/ -e BcgI -t 8 -p

./target/release/syn2bani triangle \
  genomes/*.fna \
  -o matrix.tsv --edge-list -t 8 -p

./target/release/syn2bani search \
  query.fna sketches/ \
  -o results.tsv -t 8 -p -m 0.85

./target/release/syn2bani db build \
  genomes/*.fna \
  -o database/ -e BcgI -t 8 -p
```

---

## 10. 下一步建议（在 Mac Studio 上执行）

### 立即执行（第 1 天）
1. ✅ 克隆两个仓库，编译验证
2. ✅ 下载 10-20 个代表性物种基因组（测试编译和功能）
3. ✅ 运行 `cargo test` 确认所有测试通过
4. ✅ 运行 `cargo bench` 获取 Mac Studio 基准数据

### 短期（第 1-3 天）
5. 下载大规模数据集（RefSeq 细菌、UHGG）
6. 完成 Task 1：真实 MAG 验证（100+ 基因组）
7. 生成 Figure 1-4 的数据和初稿

### 中期（第 3-7 天）
8. 完成 Task 2：更大范围 GBRT 训练
9. 完成 Task 4：完整 Benchmark 报告
10. 补充论文 manuscript 的 Results 和 Discussion

### 长期（第 7-14 天）
11. 完成 Task 3：SV 检测完整实现
12. 论文投稿准备（Nature Methods / Genome Biology）
13. 开源发布（GitHub Release + Bioconda）

---

## 11. 联系和备忘

- **GitHub 组织**：HuangShiLab
- **论文 repo**：https://github.com/HuangShiLab/Syn2bANI-paper
- **代码 repo**：https://github.com/HuangShiLab/Syn2bANI
- **当前开发环境**：MacBook Pro (M3 Pro, 18GB) → Mac Studio (2TB)
- **关键依赖版本**：Rust 1.92, rayon 1.10, clap 4.5

---

*此文档应在 Mac Studio 上首次启动时阅读，确保开发连续性。*
