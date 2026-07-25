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

**最新提交**（2026-07-25 更新）：
- CODE: `3048177` — feat: chain-restricted stratified MLE ANI (`syn2bani ani`)
- paper: 见本次提交

**当前推荐路径**：`syn2bani ani`（§4.6）。`dist` + GBRT 是旧路径，
已知缺陷记录在 `V8_MLE_VALIDATION.md` §5。下一步首要任务是 §5 Task 0。

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

### 4.5 v7 算法：TGT sparse chaining + chain 内 k-mer ANI + GBRT v7 校准模型
- **实现时间**：2026-07-25
- **核心改动**：
  - `TagSet` 携带原始 contig 序列 (`sequences: Vec<Vec<u8>>`)
  - `SyntenyBuilder` 改为按 `(q_contig_id, r_contig_id, orientation)` 分组，组内用 DP 稀疏链式算法（indel tolerance 5000 bp）构建共线性 block
  - `AniCalculator` 在每个 synteny block 的 query/reference 区间内用 canonical k-mer（默认 k=15）计算局部 ANI，按 anchor tag 数加权平均
  - `dist` / `struct` 默认输出改为 `chained_kmer_ani`
  - raw-features TSV 新增 `chained_kmer_ani` 列
  - 新增 GBRT v7 校准模型：features = `[raw_ani, mash_ani, chained_kmer_ani, shared_log, af_q, af_r]`，在 GTDB-R207 728 对有效 pair 上训练
- **对应脚本**：
  - 实现/验证：`benchmark_v7_mid_ani.py`、`slurm_v7_mid_ani.sh`
  - 特征提取/训练：`extract_syn2bani_features.py`、`train_gbrt_v7.py`、`submit_v7_feature_extraction_and_training.sh`
- **15 对 mid-ANI（85–95%）口腔/肠道验证结果 vs FastANI**：
  - skani：MAE 0.468%
  - Syn2bANI mash_ani：MAE 2.808%
  - Syn2bANI chained_kmer_ani（默认）：MAE 2.971%
  - Syn2bANI -0.028 经验校准：MAE 1.031%
  - **Syn2bANI GBRT v7：MAE 1.142%**
- **结论（2026-07-25 修订）**：这条路线已被 4.6 取代。`chained_kmer_ani` 的
  Pearson r = −0.107 不是"block 偏保守区"造成的（系统偏差会保留相关性），
  而是 block 区间本身算错了——详见 `V8_MLE_VALIDATION.md` §5。
  v2–v7 模型文件保留作为对照，不再继续扩大训练集。

### 4.6 v8 算法：chain-restricted stratified MLE（当前推荐路径）

- **实现时间**：2026-07-25
- **代码提交**：Syn2bANI `3048177`
- **详细文档**：代码仓库 `ALGORITHM_MLE.md`（推导 + CLI），
  本仓库 `V8_MLE_VALIDATION.md`（验证记录 + 论文影响）
- **核心思想**：`raw_ani` 和 `mash_ani` 不是两个待回归的特征，而是**同一个
  截断二项似然的两个矩**。拟合该似然即可，不需要学习组合，也不会在多酶精确
  匹配把 `raw_ani` 钉到 1 时退化。
- **新增模块**（全部是新增路径，`dist`/`sketch`/`db`/`search`/`struct` 未改动）：
  - `src/core/mle.rs` — 按酶分层的截断二项 MLE，外加 loss-only / histogram-only
    两个偏估计量和一个按存留率门控的一致性检查
  - `src/core/chain_ani.rs` — 鸽巢原理容忍种子、按
    `(q_contig, r_contig, orientation)` 分组的带 gap penalty + max_gap 的
    chaining DP、链内 seed-and-extend 局部填充
  - `tag_extractor.rs` — `revcomp_packed` / `canonical_packed` /
    `GenomeTag::canonical()`，修复倒位段和反向 draft contig 无法匹配的问题
  - `src/cli/ani.rs` — `syn2bani ani` 子命令
  - `prototype/` — Python ground-truth 工具链（生成物约 90 MB，不入库）
- **验证结果（精确构造真值，*E. coli* K-12）**：
  - 参考基因组固定为 ENA `U00096.3`（完整 MG1655，4,641,652 bp），
    用 `prototype/fetch_reference.sh` 下载，保证逐字节可复现
  - 12 个基因组，85–99.9% ANI：**MAE 0.053%**，零训练数据、零校准，全部 flag=ok
  - 6 个基因组，true ANI 固定 95%、accessory 0→50%：**MAE 0.114%**，
    估计值平坦且不随 accessory 单调漂移，AF 精确跟踪 `1−F`（误差 ≤ 0.004）；
    同批基因组上全基因组 containment 漂移 95.18 → 93.27
  - 37/37 lib 测试通过
- **链断裂判据（2026-07-25 修订）**：改为数"跳过多少个 query tag 位点"，
  阈值 `j* = ln(α)/ln(1−p)` 由拟合分歧度两趟确定。固定 bp 的 `max_gap` 在原理上
  不可行——详见 `V8_MLE_VALIDATION.md` §3.4。**真实 accessory 区段（前噬菌体、
  基因组岛）通常远小于 50 kb，旧默认值会把它们吞进链里系统性压低 ANI**，
  Task 0 里要专门核对 AF 列。
- **⚠️ 尚未验证**：indel、非 *E. coli* 物种、真实基因组对（未与 FastANI/skani
  直接比较）。~85% 以下一致性交叉检查自动关闭。小于 5 个 tag 位点的等长
  accessory 区段仍会被跨过。
- **用法**：
  ```bash
  syn2bani ani <queries> <reference> --verbose -p -o mle.tsv
  # 列：query reference ani af_query af_reference std_err
  #     + --verbose: ani_from_loss ani_from_hist n_anchors n_chains n_tags flag
  ```
  `--verbose` 的诊断列用于定位真实数据上的问题：`n_anchors`/`n_chains` 看
  chaining 是否成功，两个偏估计量 + `flag` 看模型是否拟合观测。

---

## 5. 剩余开发任务（优先级排序）

> 2026-07-25 修订：Task 0 新增并置顶；原 Task 2（扩大 GBRT 训练）已废弃。

### 🔴 高优先级（数据密集型，Mac Studio / HPC 优势）

#### Task 0: 在真实基因组上验证 v8 MLE 路径 ★部分完成 2026-07-25

- **已完成部分**：13 条肠杆菌科染色体对 *E. coli* K-12，三方对比 skani + FastANI。
  发现速率异质性是真实数据上的主要偏差来源，已用 gamma 混合模型修复
  （详见 `V8_MLE_VALIDATION.md` §3.6）。可报告的 8 对上 MAE 0.314（vs skani）
  / 0.256（vs FastANI）。用 `prototype/realgenome_bench.sh` 可一键复现。
- **仍待完成**：你那 15 对口腔/肠道基因组；draft assembly / MAG；高低 GC 类群；
  换用 ANIm 或 minimap2 作为真值而非 FastANI
- **目标**：确认 `ani` 在真实数据上是否复现仿真上的精度
- **步骤**：
  1. `git pull` 代码仓库，`cargo build --release`，`cargo test --release --lib`
     （应为 37/37）
  2. 先跑 `prototype/` 两组仿真复现 MAE 0.053% / 0.114%，确认环境一致
     （先 `bash fetch_reference.sh mg1655.fasta` 拿到固定的参考基因组）
  3. 在既有 15 对 mid-ANI（85–95%）口腔/肠道基因组上跑
     `syn2bani ani ... --verbose`，与 skani / GBRT v7 并列比较
  4. **换掉真值来源**：不要继续用 FastANI（< 92% 区间可靠性有争议，且拿它当
     真值又和 skani 比有方法学偏差）。改用 nucmer/MUMmer（ANIm）或
     minimap2 比对得到的真值
- **诊断优先于调参**：如果真实数据上明显变差，先看 `--verbose` 的
  `n_anchors` / `n_chains` / `flag` 三列定位是 chaining 失败还是模型不拟合，
  **不要加经验偏移**（见 `V8_MLE_VALIDATION.md` §1.2：常数偏移等价于硬编码
  一个固定的共享含量比例）
- **输出**：真实数据验证表，补进 `V8_MLE_VALIDATION.md` §3.3

#### Task 0b: 补齐仿真覆盖面

- **indel**：`simulate.py` 已支持 `indel_rate` 但当前跑的是 0，需要开启以压测
  gap 算术路径
- **多物种**：除 *E. coli* 外增加高 GC（如 *Streptomyces*）和低 GC（如
  *Staphylococcus*）各一个，检查 GC 偏好对 tag 密度和 MLE 的影响
- **碎片化**：把参考切成 20/50/100/200 contig 并随机翻转方向，验证
  `canonical_packed` 修复确实生效（旧代码在这个场景下会丢掉约一半共享 tag）

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

#### ~~Task 2: 更大范围 GBRT 训练~~（已废弃 2026-07-25）

被 4.6 的 MLE 路径取代。废弃理由不是"模型不够大"，而是这个建模框架本身有问题：

- GBRT v7 in-sample MAE 0.238% → out-of-sample 1.142%，**约 5 倍过拟合/分布外**
- 特征重要性最高的 `raw_ani`（0.3197）实测对真值斜率只有 **0.336**，被
  `near_match_tolerance = 2` 的截断钉死在 `≥ 1 − 2/32`
- 喂进去的 `containment` 特征是从 `raw_ani` 伪造的，实际是常数，零信息量
- 模型真正在学的是"用 `af_q` 反推 accessory 比例造成的偏差"，而这个偏差
  在链内计数下根本不会出现

扩大训练集只会让模型更精确地拟合一个错误参数化。v2–v7 模型文件保留在代码仓库
中作为论文对照。

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
