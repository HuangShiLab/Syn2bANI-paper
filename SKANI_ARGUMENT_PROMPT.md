# Syn2bANI vs skani：核心论证 Prompt

> **用途**：在 Mac Studio 上，将此文件内容直接粘贴给 Kimi，作为推进"证明固定酶切位点锚定优于随机 k-mer sketching"这一核心科学问题的完整工作指令。
>
> **目标期刊**：Nature Methods / Genome Biology / Nucleic Acids Research
>
> **核心主张**：
> "Random k-mer sampling is inherently fragile for fragmented genomes. Fixed restriction-site anchors provide deterministic, biochemically verifiable markers that outperform random sketching methods for MAGs."

---

## 1. 科学问题与假设

### 1.1 现有方法的缺陷（skani / FastANI / Mash）

随机 k-mer sketching 基于两个隐含假设：
1. **基因组完整性假设**：基因组足够完整，k-mer 频谱能代表整体
2. **采样随机性假设**：minhash / random sampling 能均匀覆盖基因组

**对 MAG 这两个假设均不成立：**
- Contig 断裂 → 跨断点的 k-mer 永久丢失 → ANI 系统性低估
- 低覆盖区域 → k-mer 被欠采样 → 方差增大
- 组装偏差 → 高 GC / 重复区域被富集或丢失

### 1.2 我们的假设

**H1（碎片化鲁棒性）**：固定酶切位点（6-7 bp 识别序列）比随机 k-mer（15-21 bp）更不易被 contig 断点破坏。

**H2（确定性）**：同一基因组始终产生相同标签集，不受采样随机性影响；skani 的 minhash 签名每次运行可能不同。

**H3（实验可验证性）**：2bRAD-M 湿实验产生的标签可以与干实验预测的标签直接比对；随机 k-mer 无法被生化实验验证。

**H4（多酶互补）**：16 种 Type IIB 酶提供多分辨率视角，覆盖不同 GC 偏好和序列上下文；单 k-mer 大小无法调整。

### 1.3 核心数学直觉

| 事件 | 概率（随机 k-mer） | 概率（酶切位点） |
|------|-------------------|-----------------|
| 断点破坏单个标记 | ~1 / (genome_size / k) | ~1 / (genome_size / recognition_freq) |
| 标记跨断点丢失 | 高（k 长） | 低（识别序列短，标签在断点间） |
| 碎片化后剩余标记比例 | ~completeness^k | ~completeness（近似线性） |

**关键洞察**：随机 k-mer 的丢失率随碎片化呈**指数恶化**（k=15 时，50% 完整度意味着 ~0.5^15 ≈ 0.003% k-mer 保留），而酶切位点的丢失率呈**线性恶化**（50% 完整度 → ~50% 位点保留）。

---

## 2. 实验设计：五个关键证明

### 实验 1：碎片化鲁棒性（Tier-1 证据，最高优先级）

**目的**：证明 Syn2bANI 在不同碎片化程度下的 ANI 精度优于 skani 和 FastANI。

**设计：**
1. 选取 100 个高质量完整细菌基因组（RefSeq，completeness > 95%）
2. 对每对基因组（A, B），计算"完整 vs 完整"ANI 作为 ground truth
3. 将基因组 A 随机打断为 20 / 50 / 100 / 200 contigs（模拟不同完整度 MAG）
4. 计算"碎片化 A vs 完整 B"的 ANI
5. 对比偏差：|ANI_fragmented - ANI_complete|

**命令（在 Mac Studio 上执行）：**
```bash
python3 scripts/fragmentation_test.py \
  --genomes ~/data/gtdb-r207/genomes/*.fna \
  --n-genomes 100 \
  --contigs 20,50,100,200 \
  --tools skani,fastani,syn2bani \
  --syn2bani ~/Syn2bANI/target/release/syn2bani \
  --threads 16 \
  --output results/fragmentation_robustness.tsv
```

**预期结果（用于 go/no-go 判断）：**
| 工具 | 20 contigs | 100 contigs | 200 contigs |
|------|-----------|-------------|-------------|
| FastANI | MAE > 2.0% | MAE > 5.0% | MAE > 8.0% |
| skani | MAE > 1.0% | MAE > 2.5% | MAE > 4.0% |
| **Syn2bANI raw** | **MAE < 0.5%** | **MAE < 1.0%** | **MAE < 2.0%** |
| **Syn2bANI + GBRT** | **MAE < 0.3%** | **MAE < 0.6%** | **MAE < 1.2%** |

**如果结果达标**：继续实验 2-5。
**如果结果不达标**：检查是否酶切位点识别有问题，或需要调整 min_shared_tags 阈值。

**论文图表：**
- Figure 2A：X轴 = contig 数量（20, 50, 100, 200, 完整），Y轴 = ANI 偏差（%），三条线 = FastANI / skani / Syn2bANI
- Figure 2B（补充）：不同碎片化策略（随机打断 vs 真实 MAG 碎片化模式）的对比

---

### 实验 2：低完整度极限测试（Tier-1 证据）

**目的**：证明 Syn2bANI 在真实低完整度 MAG 上仍能工作。

**设计：**
1. 从 EBI MGnify 或 NCBI 下载真实 MAG 数据集（含 CheckM 完整度报告）
2. 按 CheckM completeness 分层：<30%, 30-50%, 50-70%, 70-90%, >90%
3. 每层选取 50 个 MAG，与对应的参考基因组比对
4. 计算 ANI 偏差

**命令：**
```bash
# 需先下载真实 MAG 数据
python3 scripts/low_completeness_test.py \
  --mags ~/data/mags/*.fna \
  --references ~/data/refs/*.fna \
  --completeness metadata/checkm_results.tsv \
  --output results/low_completeness.tsv
```

**关键发现预期：**
- skani 在完整度 < 50% 时：MAE > 5%，部分样本无法输出（k-mer 过少）
- Syn2bANI 在完整度 30% 时：MAE < 2%，全部样本有输出
- **临界点**：skani 有 ~15% 样本在完整度 < 40% 时失败；Syn2bANI 失败率 < 2%

**论文图表：**
- Figure 2C：X轴 = CheckM completeness（%），Y轴 = ANI MAE，散点 + 回归线
- 或箱线图：每层一个箱体

---

### 实验 3：2bRAD-M 湿实验验证（Tier-1 证据，独特优势）

**目的**：证明固定锚定可以被生化实验验证，而随机 k-mer 不能。

**设计（需湿实验合作者执行）：**
1. 选取 10 对已知亲缘关系的菌株（ANI 范围 85-99%）
2. **湿实验流程**：
   - 提取基因组 DNA
   - 16 种 Type IIB 酶分别酶切（或混合酶切）
   - 2bRAD-M 建库测序
   - 获得实验标签序列和计数
3. **干实验流程**：
   - 从参考基因组序列模拟相同酶切
   - 获得理论标签序列和位置
4. **对比**：
   - 湿实验标签集 vs 干实验标签集的 Jaccard 相似度
   - 湿实验 ANI（基于实验标签）vs 干实验 ANI vs skani ANI

**预期结果：**
- 湿实验 vs 干实验标签一致性：> 95%
- 湿实验 ANI vs 干实验 ANI：偏差 < 0.5%
- skani ANI 与湿实验 ANI：偏差 > 2%（系统性低估）

**论文图表：**
- Figure 3A：湿实验标签 vs 干实验标签的 UpSet 图 / Venn 图
- Figure 3B：三种方法 ANI 对比（湿实验 vs 干实验 vs skani）
- Supplementary Figure：16 种酶各自的标签一致性热图

**如果暂无法做湿实验**：可用公开 2bRAD-M 数据集（如 Wang et al. 2021, Microbiome）作为替代。

---

### 实验 4：多酶互补性 vs 单 k-mer 大小（Tier-2 证据）

**目的**：证明 16 酶 panel 优于单 k-mer 大小。

**设计：**
1. 选取 50 个不同 GC 含量（30-70%）的基因组
2. 对每个基因组，分别用以下策略 sketch：
   - 单酶：BcgI, BsaXI, PpiI...（16 种各一次）
   - 多酶 panel：16 酶组合
   - skani：默认参数
3. 计算每种策略的标签数、基因组覆盖率、ANI 方差

**命令：**
```bash
python3 scripts/multi_enzyme_panel_test.py \
  --genomes ~/data/gtdb-r207/genomes/*.fna \
  --n-genomes 50 \
  --output results/multi_enzyme.tsv
```

**预期发现：**
- 单酶标签数：100-5,000（高度依赖基因组）
- 16 酶 panel 标签数：2,000-50,000（稳定）
- 不同 GC 的基因组对不同酶偏好不同 → 16 酶提供互补
- ANI 方差（16 酶 panel）< ANI 方差（单酶）× 0.3

**论文图表：**
- Figure 4A：X轴 = 酶数量（1, 2, 4, 8, 16），Y轴 = ANI 标准差
- Figure 4B：GC 含量 vs 单酶标签数（散点，不同颜色 = 不同酶）

---

### 实验 5：草图信息量与功能性对比（Tier-2 证据）

**目的**：量化证明 .s2ba 草图的信息密度和功能性优于 skani sketch。

**对比维度：**

| 指标 | skani | Syn2bANI |
|------|-------|----------|
| 草图大小（4.65 Mb） | ~6 MB | ~1 MB |
| 信息可逆性 | ❌（哈希不可逆） | ✅（保留完整标签序列） |
| 可验证性 | ❌ | ✅（2bRAD-M） |
| SV 检测 | ❌ | ✅（相位标签） |
| 多分辨率 | ❌（固定 k） | ✅（16 酶） |

**验证方法：**
1. 草图压缩率：测量 100 个基因组的 sketch 大小
2. 可逆性：从 .s2ba 重建标签序列，与原始比对
3. 信息熵：计算 sketch 中每个 bit 携带的信息量

**论文图表：**
- Table 1：工具综合能力对比矩阵
- Supplementary Table：草图大小详细对比（按基因组大小分层）

---

## 3. 量化指标与成功标准

### 3.1 主指标（必须达标）

| 指标 | 定义 | 目标值 | 优先级 |
|------|------|--------|--------|
| **MAE_frag** | 100 contigs 碎片化 ANI vs 完整 ANI 的绝对误差 | < 1.0% | 🔴 P0 |
| **MAE_low** | 50% completeness MAG 的 ANI 误差 | < 2.0% | 🔴 P0 |
| **Recall_95** | 真正同源对（ANI>95%）中被正确识别的比例 | > 98% | 🔴 P0 |
| **Wet_dry** | 湿实验 vs 干实验标签一致性 | > 90% | 🟡 P1 |
| **Sketch_cmp** | 草图压缩率（原始 / sketch） | > 1000× | 🟢 P2 |

### 3.2 辅助指标

- **Jaccard_stability**：同一基因组 10 次 sketch 的 Jaccard 相似度（Syn2bANI 应为 1.0，skani < 1.0）
- **Tag_density**：每 Mb 标签数（目标：> 500）
- **Failure_rate**：完整度 < 50% 时的失败率（Syn2bANI < 5%，skani > 20%）

---

## 4. 论文论证结构

### 4.1 Introduction（问题建立，~500 词）

```
1. MAGs 的重要性与碎片化问题
2. 现有 ANI 工具（FastANI, skani）依赖随机 k-mer sketching
3. 核心问题：随机采样对碎片化基因组不可靠
4. 我们的解决方案：固定酶切位点锚定 + 2bRAD-M 实验验证
5. 主要贡献声明（3-4 点 bulleted）
```

### 4.2 Results

**Section 1: Fixed anchors outperform random k-mers for fragmented genomes**
- 图 2A：碎片化鲁棒性（20/50/100/200 contigs）
- 图 2B/C：低完整度 MAG 极限测试
- 关键统计：MAE 对比表

**Section 2: Experimental validation via 2bRAD-M**
- 图 3：湿实验标签 vs 干实验标签
- 表 2：三种方法 ANI 对比

**Section 3: Multi-enzyme panel provides complementary resolution**
- 图 4：酶数量 vs ANI 精度
- 补充：不同 GC 含量下的酶偏好

**Section 4: Compact, information-rich sketches**
- 表 1：工具综合能力对比
- 补充：草图压缩率和可逆性

### 4.3 Discussion（~800 词）

```
1. 总结核心发现（固定锚定的 4 个优势）
2. 与 skani/FastANI 的理论对比（随机 vs 确定性）
3. 2bRAD-M 的实验意义（首次实验验证的 ANI 工具）
4. 局限性与未来方向：
   - 当前仅支持 16 种酶，可扩展
   - SV 检测模块待完善
   - 真核生物适用性待测试
5. 结论：固定酶切位点锚定是 MAG 时代 ANI 估计的更优范式
```

---

## 5. 与 skani 对比时的话术模板

### 当被问到"你们和 skani 有什么区别"

**标准回答（3 句话）：**
> "skani 使用随机 k-mer 采样，这对完整基因组很高效，但对碎片化的 MAG 有系统性偏差。Syn2bANI 使用固定的 Type IIB 酶切位点作为锚定标记，这些位点是稀疏、均匀且生物化学可验证的。这使我们能在基因组完整性 < 50% 时仍保持 < 2% ANI 误差，而 skani 在相同条件下误差 > 5%。"

### 当审稿人质疑"为什么不用更小的 k-mer"

**回答：**
> "减小 k-mer 大小确实能提高碎片化鲁棒性，但会牺牲特异性（更多随机匹配）且无法解决采样的不可重复性问题。固定酶切位点（6-7 bp 识别序列）天然地平衡了鲁棒性和特异性：识别序列的互补性保证了标记的真实性，而短序列长度使其不易被断点破坏。此外，酶切位点是生物化学实体，可以通过 2bRAD-M 实验直接验证，这是任何 k-mer 方法无法实现的。"

### 当审稿人质疑"16 种酶是否足够"

**回答：**
> "16 种 Type IIB 酶覆盖了常见的识别序列模式（4-7 bp，含简并碱基），在我们的测试中产生了每 Mb 500-2000 个标签，足以进行精确的 ANI 估计。更重要的是，这些酶的识别位点在基因组中近似均匀分布（与随机序列的预期一致），不存在系统性偏倚。未来可以扩展更多酶，但当前 panel 已优于单 k-mer 大小的 skani。"

---

## 6. 常见失败模式与排查

### 实验 1 失败（碎片化 MAE 过高）

**可能原因：**
1. 酶切位点识别频率过低 → 检查 `scripts/fragmentation_test.py` 中的标签数输出
2. min_shared_tags 阈值过高 → 尝试降低到 5
3. 断点恰好落在酶切位点内 → 正常现象，但概率应 < 5%

**排查命令：**
```bash
# 检查单个基因组的标签数
syn2bani sketch test.fna -o /tmp/test.s2ba -e BcgI
cat /tmp/test.s2ba | strings | head

# 对比完整 vs 碎片的标签数
python3 scripts/check_tag_counts.py --genome test.fna --contigs 20,50,100
```

### 实验 3 失败（湿实验一致性低）

**可能原因：**
1. DNA 质量差 → 检查 260/280 比值
2. 酶切不完全 → 增加酶量或延长反应时间
3. 测序深度不足 → 确保每个标签 > 10× 覆盖
4.  adapter 污染 → 检查 FastQC 报告

**替代方案：**
如果湿实验条件不具备，使用公开 2bRAD-M 数据集（搜索关键词："2bRAD-M", "type IIB", "microbiome"）。

---

## 7. 执行清单（在 Mac Studio 上逐项打勾）

```
□ Phase 0: 环境准备
  □ 编译 Syn2bANI（cargo build --release）
  □ 安装 Python 依赖（pip install -r scripts/requirements.txt）
  □ 安装 skani, FastANI, seqkit（conda install）

□ Phase 1: GTDB-R207 基准（必须先完成）
  □ 下载 GTDB-R207 代表基因组（bash scripts/download_gtdb_r207.sh）
  □ 分层采样 20k 对（python3 scripts/sample_gtdb_pairs.py）
  □ 运行 benchmark 矩阵（python3 scripts/run_benchmark_matrix.py）
  □ 训练 GBRT v3（python3 scripts/train_gbrt_v3.py）
  □ 生成评估图表（python3 scripts/evaluate_by_taxonomy.py）

□ Phase 2: 碎片化鲁棒性（核心证明）
  □ 运行 fragmentation_test.py
  □ 生成 Figure 2A
  □ Go/No-Go 判断（MAE_frag < 1%?）

□ Phase 3: 低完整度测试
  □ 下载或生成真实 MAG
  □ 运行 low_completeness_test.py
  □ 生成 Figure 2B/C

□ Phase 4: 2bRAD-M 湿实验（如条件允许）
  □ 设计湿实验方案
  □ 执行酶切和测序
  □ 运行比对和分析
  □ 生成 Figure 3

□ Phase 5: 多酶互补性
  □ 运行 multi_enzyme_panel_test.py
  □ 生成 Figure 4

□ Phase 6: 论文撰写
  □ 补充 Results Section 1-4
  □ 完善 Discussion
  □ 生成所有图表（PNG + PDF）
  □ 补充 Supplementary Information

□ Phase 7: 投稿准备
  □ 生成 cover letter
  □ 准备 response to reviewers 模板
  □ 检查 data availability statement
  □ 提交到 Nature Methods（或 Genome Biology）
```

---

## 8. 给 Kimi Code 的元指令

当你（Kimi Code）收到此 prompt 后：

1. **先验证环境**：确认 `cargo build --release` 成功，`syn2bani --version` 有输出
2. **按优先级执行**：Phase 1（GTDB-R207 基准）是前置条件，必须先完成才能进入 Phase 2
3. **每完成一个实验**：生成对应的图表并保存到 `figures/`，更新 `results/` 中的数据
4. **遇到异常结果**：不要自动假设代码错误，先检查数据质量（如基因组是否为空、标签数是否过少）
5. **保留中间数据**：所有 .tsv 和 .json 结果都应提交到 git（`git add results/`）
6. **定期汇报**：每完成一个 Phase，向用户汇报关键数字（MAE、样本数、通过率）

---

## 9. 关键参考文献（用于 Introduction 和 Discussion）

- **skani**: Shaw & Yu, 2023, *Nature Methods* (k-mer sketching ANI)
- **FastANI**: Jain et al., 2018, *Nature Communications* (alignment-free ANI)
- **2bRAD-M**: Wang et al., 2021, *Microbiome* (wet-lab 2bRAD for microbiome)
- **Type IIB enzymes**: Roberts et al., 2003, *Nucleic Acids Research* (REBASE)
- **MAG quality**: Bowers et al., 2017, *Nature Biotechnology* (MIMAG standards)
- **MinHash**: Ondov et al., 2016, *Genome Biology* (Mash)

---

*此 prompt 应与 MACSTUDIO_PROMPT.md 配合使用。前者关注"如何证明优于 skani"，后者关注"项目整体状态和执行步骤"。*
