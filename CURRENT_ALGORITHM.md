# Syn2bANI 当前核心算法公式梳理

> 基于代码逐行分析：`tag_extractor.rs` → `tag_matcher.rs` → `synteny_builder.rs` → `ani_calculator.rs`

---

## 一、算法流程图

```
Genome A (FASTA)          Genome B (FASTA)
    ↓                          ↓
Digestion                  Digestion
(酶切提取标签)              (酶切提取标签)
    ↓                          ↓
TagSet A: [t₁, t₂, ...]   TagSet B: [t₁', t₂', ...]
    │                          │
    └──────────┬───────────────┘
               ↓
        TagMatcher::match_tag_sets()
               ↓
    ┌──────────────────────┐
    │   1. 建索引：ref_tags  │
    │      按 packed_sequence│
    │      哈希索引          │
    │                      │
    │   2. 遍历 query_tags │
    │      找最佳匹配        │
    │      (Hamming distance)│
    └──────────────────────┘
               ↓
    MatchedPairs: [(qᵢ, rⱼ, hamming, local_ani)]
               ↓
    SyntenyBuilder::build_blocks()
               ↓
    SyntenyBlocks: [block₁, block₂, ...]
               ↓
    AniCalculator::calculate_ani()
               ↓
    AniResult { ani, raw_ani, weighted_ani, af_q, af_r, confidence }
```

---

## 二、各步骤的数学公式

### Step 1: Digestion（酶切提取）

```
对于每种酶 E（识别序列 pattern）:
    在基因组序列 S 中扫描 pattern
    对每个识别位点:
        提取两侧 flanking sequence（tag）
        长度 = enzyme.tag_length  (BcgI: 32bp, PpiI: 27bp, ...)
        方向 = '+' 或 '-'
        
    GenomeTag = {
        position: 位点在基因组中的绝对位置,
        sequence: [u8; 32],     // 实际碱基序列
        packed_sequence: u64,   // 2-bit 压缩 (A=00, T=01, C=10, G=11)
        seq_len: 实际长度,
        direction: '+'/'-',
        enzyme: "BcgI",
    }
```

**关键性质：**
- 标签数是稀疏的：~10,000-50,000 个 / 细菌基因组
- 位置是固定的：同一基因组始终产生相同标签集
- 序列是完整的：保留了实际 DNA 序列（不是哈希值）

---

### Step 2: Matching（标签匹配）

#### 2.1 索引构建

```
ref_index: HashMap<u64 → Vec<usize>>  // packed_sequence → 标签索引列表

for i, tag in reference.tags:
    ref_index[tag.packed_sequence].push(i)
```

#### 2.2 遍历匹配

```
for each q_tag in query.tags:
    
    // 第一层：精确匹配（packed_sequence 相同）
    if q_tag.packed_sequence in ref_index:
        candidates = ref_index[q_tag.packed_sequence]  // 候选列表
        
        // 在候选中找 Hamming 距离最小的
        best_dist = ∞
        best_idx = None
        for idx in candidates:
            if matched_ref_flags[idx] already used: skip
            dist = hamming_distance(q_tag, ref_tags[idx])
            if dist < best_dist:
                best_dist = dist
                best_idx = idx
        
        if best_idx exists:
            accept = (best_dist ≤ tolerance)  // tolerance = 2 (默认)
            if accept:
                create MatchedPair
                continue
    
    // 第二层：近匹配回退（若精确匹配失败）
    if allow_near_match:
        遍历所有未匹配的 ref_tags
        找全局 Hamming 距离最小的
        if best_dist ≤ tolerance:
            create MatchedPair
            continue
    
    // 未匹配
    unmatched_query.push(q_tag)
```

#### 2.3 Hamming distance 计算

```
hamming_distance(a: GenomeTag, b: GenomeTag):
    cmp_len = min(a.seq_len, b.seq_len)  // 比较长度取最小
    xor = a.packed_sequence XOR b.packed_sequence
    mask = (1 << (cmp_len * 2)) - 1      // 只比较有效位
    return popcount(xor & mask)          // 不同碱基数
```

**注意：** Hamming distance 只记录 **substitutions**（替换），不记录 indels（因为标签长度固定）。

---

### Step 3: Synteny Block 构建

```
输入：按 query 位置排序的 MatchedPairs
输出：SyntenyBlocks

初始化：
    current_start = 0
    orientation = pairs[0].query_tag.direction

for i = 1 to n-1:
    prev = pairs[i-1]
    curr = pairs[i]
    
    query_gap = curr.query.position - prev.query.position
    ref_gap = curr.ref.position - prev.ref.position
    gap_diff = |query_gap - ref_gap|
    
    // 检测断点条件：
    condition_1 = (curr.direction ≠ orientation)     // 方向翻转
    condition_2 = (query_gap > 10000)                // query 大间隙
    condition_3 = (ref_gap > 10000)                  // ref 大间隙
    condition_4 = (gap_diff > 5000)                  // 间隙差异大
    
    if any(condition_1, condition_2, condition_3, condition_4):
        // 结束当前 block
        block = create_block(pairs, current_start, i-1, orientation)
        blocks.push(block)
        
        // 开始新 block
        current_start = i
        orientation = curr.direction

// 最后一个 block
block = create_block(pairs, current_start, n-1, orientation)
blocks.push(block)
```

**Block 的属性：**
```
SyntenyBlock = {
    query_start: block 内最小 query 位置,
    query_end:   block 内最大 query 位置,
    ref_start:   block 内最小 ref 位置,
    ref_end:     block 内最大 ref 位置,
    matched_tags: block 内标签数,
    orientation:  '+' 或 '-',
    block_ani:    block 内所有 local_ani 的平均,
}
```

---

### Step 4: ANI 计算（核心公式）

#### 4.1 基础量定义

```
N_match = matched_pairs.len()          // 匹配上的标签对数
N_q_unmatch = unmatched_query.len()    // query 未匹配标签数
N_r_unmatch = unmatched_ref.len()      // ref 未匹配标签数
N_q_total = N_match + N_q_unmatch      // query 总标签数
N_r_total = N_match + N_r_unmatch      // ref 总标签数

AF_query    = N_match / N_q_total      // query 的 alignment fraction
AF_reference= N_match / N_r_total      // ref 的 alignment fraction

local_ani_i = 1.0 - (hamming_dist_i / tag_length_i)   // 第 i 对的局部 ANI
```

#### 4.2 Raw ANI（未修正）

```
RAW_ANI = mean(local_ani_i) = (Σ local_ani_i) / N_match

等价展开：
RAW_ANI = (1/N_match) × Σ_i [1 - (hamming_dist_i / tag_length_i)]
        = 1 - (1/N_match) × Σ_i (hamming_dist_i / tag_length_i)
```

**物理意义：** 匹配上的标签中，平均有多少比例的碱基是相同的。

#### 4.3 Weighted ANI（加权平均）

```
weights_i = weight_strategy(pairs_i, blocks)

WEIGHTED_ANI = Σ_i (local_ani_i × weights_i) / Σ_i weights_i
```

四种权重策略：

| 策略 | 公式 | 含义 |
|------|------|------|
| **Uniform** | `w_i = 1.0` | 所有对平等 |
| **Synteny** | `w_i = sqrt(block_length)` | 长 synteny block 内的标签权重更高 |
| **Position** | `w_i = 1 + sin(norm_position)` | 基因组中间位置权重略高 |
| **GapAdjusted** | `w_i = 1 - 0.5 × min(|gap_diff|, 10)/10` | 间隙差异大的标签权重降低 |

#### 4.4 Simple Debias（多项式修正）

```
ANI_simple = RAW_ANI + 0.02 × (1 - RAW_ANI) × (1 - min(AF_query, AF_reference))
```

**直觉：**
- 当 `AF` 低时（匹配标签少），ANI 可能被低估，需要向上修正
- 当 `RAW_ANI` 已经很高（接近 1）时，修正量小
- 系数 0.02 是经验值

#### 4.5 GBRT Debias（梯度提升树修正）

```
输入特征（7维）：
    x₁ = RAW_ANI                  // 原始 ANI
    x₂ = AF_query                 // query AF
    x₃ = AF_reference             // ref AF
    x₄ = shared_tags ≈ RAW_ANI × min(N_q_total, N_r_total)  // 共享标签数
    x₅ = containment = shared_tags / max(N_q_total, N_r_total) // 包含度
    x₆ = div_proxy = 1 - RAW_ANI  // 分歧代理
    x₇ = ref_gc = 0.5             // 参考 GC（默认，特征重要性≈0）

GBRT_ANI = GBRT_model.predict(x₁, x₂, x₃, x₄, x₅, x₆, x₇)
         = init_value + learning_rate × Σ_tree(leaf_value_traversed)
```

模型结构：
- 300 棵决策树，深度 5
- 每棵树由一系列 split nodes + leaf nodes 组成
- 预测时从根节点遍历到叶子，累加 leaf value
- 最终输出经过 clamp（限制在 [0, 1] 内）

#### 4.6 Confidence（置信度）

```
confidence = (1 - exp(-N_match / 100)) × sqrt(min(AF_query, AF_reference))
```

**直觉：**
- 匹配标签越多 → 置信度越高（指数饱和）
- AF 越高 → 置信度越高（平方根）

---

## 三、当前输出的所有字段

### TSV 输出（默认）

```
query_file    ref_file    ani      af_q    af_r    query_name    ref_name    shared_tags    sv_count
```

其中：
- `ani` = `AniResult.ani` = GBRT 修正后的最终 ANI（或 simple debias）
- `af_q` = `AniResult.af_query` = query 的 alignment fraction
- `af_r` = `AniResult.af_reference` = ref 的 alignment fraction
- `shared_tags` = `local_ani_profile.len()` = matched_pairs 数量

### 内部 AniResult 结构

```rust
AniResult {
    ani: f64,                    // 最终输出（debiased）
    raw_ani: f64,                // RAW_ANI（未修正）
    af_query: f64,               // N_match / N_q_total
    af_reference: f64,           // N_match / N_r_total
    weighted_ani: f64,           // 加权平均（当前默认 Uniform，所以 ≈ raw_ani）
    confidence: f64,             // (1 - e^(-N/100)) × sqrt(min(AF))
    local_ani_profile: Vec<f64>, // 每个 matched_pair 的 local_ani 列表
}
```

### JSON 输出（更详细）

```json
{
    "query": "genomeA",
    "reference": "genomeB",
    "ani": 0.9876,               // debiased ANI
    "af_query": 0.85,
    "af_reference": 0.82,
    "weighted_ani": 0.9876,
    "confidence": 0.91,
    "shared_tags": 15234,
    "synteny_blocks": [
        {
            "query_start": 1000,
            "query_end": 50000,
            "ref_start": 2000,
            "ref_end": 51000,
            "matched_tags": 340,
            "orientation": "+",
            "block_ani": 0.989
        }
    ],
    "structural_variations": []
}
```

---

## 四、关键问题的公式级回答

### Q1: Syn2bANI 的 ANI 与 FastANI 有什么本质不同？

| 维度 | FastANI | Syn2bANI |
|------|---------|----------|
| **采样** | 密集随机（~数百万 k-mer） | 稀疏固定（~10k-50k 标签） |
| **匹配单元** | k-mer（15-21 bp） | 酶切标签（27-32 bp） |
| **匹配方式** | MinHash 近似 | Hamming 距离精确 |
| **位置信息** | ❌ 无 | ✅ 有（固定锚点位置） |
| **ANI 公式** | alignment 区域内平均 identity | 匹配标签的平均 identity |

**FastANI ANI** ≈ (所有可比对区域的 identity 之和) / (可比对区域总数)

**Syn2bANI RAW_ANI** = (Σ 匹配标签的 identity) / (匹配标签数)

**差异来源：**
1. **采样密度**：FastANI 密集 → 包含更多低保守区域；Syn2bANI 稀疏 → 只采样酶切位点附近
2. **匹配严格度**：FastANI 允许 gaps/indels（alignment）；Syn2bANI 纯 Hamming（仅 substitutions）
3. **位置约束**：Syn2bANI 要求标签在对应位置附近匹配（窗口约束）；FastANI 无此约束

### Q2: 当前公式有什么问题？

**问题 1：ANI 是混合信号**
```
RAW_ANI = mean(local_ani_i)  // local_ani_i 来自所有匹配的标签
```
这些标签中：
- 有些在 synteny block 内（保守区域）
- 有些是孤立匹配（可能是假阳性、重复序列、HGT）
- 全部被平等对待

**问题 2：GBRT 修正目标不明确**
```
GBRT_ANI = model.predict(RAW_ANI, AF_q, AF_r, shared, containment, div_proxy, ref_gc)
```
训练时的 y 是什么？如果是 `FastANI - Syn2bANI_raw`，那就是在"追逐 FastANI"。但 FastANI 本身也只是密集采样的估计量，对碎片化基因组也有偏差。

**问题 3：没有输出结构信息**
SyntenyBlock 已经计算了，但没有作为独立指标输出。用户只能看到一个 ANI 数字，丢失了"基因组骨架是否保守"的信息。

### Q3: 如何分离两个信息维度？

**维度 A：序列一致性（ANI_seq）**
```
只使用 synteny block 内的 matched pairs
ANI_seq = mean(local_ani_i for i in conserved_blocks)
```
这更接近 FastANI 的概念：在"可信任的同源区域"内计算平均 identity。

**维度 B：结构分歧（SDI）**
```
SDI = 1 - (conserved_tags / total_matched_tags)
    = (isolated_or_rearranged_tags) / total_matched_tags
```
或更精细：
```
SDI = weighted_mean(1 - anchor_adjacency_i)
其中 anchor_adjacency_i = (邻居匹配数) / (邻居总数)
```

---

## 五、建议的公式重构

### 新的数据结构

```rust
pub struct AniResultV2 {
    // 第一部分：序列层面（对标 FastANI）
    pub ani_seq: f64,              // 仅 synteny block 内的平均 identity
    pub ani_seq_corrected: f64,    // GBRT 修正后的序列 ANI
    pub n_seq_pairs: usize,        // 参与计算的标签对数
    
    // 第二部分：结构层面（Syn2bANI 独有）
    pub sdi: f64,                  // Structural Divergence Index
    pub n_conserved_blocks: usize,
    pub n_breakpoints: usize,
    
    // 辅助信息
    pub af_query: f64,
    pub af_reference: f64,
    pub confidence: f64,
}
```

### 新的计算流程

```
Phase 1: 匹配（不变）
    → MatchedPairs + SyntenyBlocks

Phase 2: 分离两个维度
    conserved_pairs = filter(pairs, |p| p in synteny_block)
    isolated_pairs  = filter(pairs, |p| p not in synteny_block)

Phase 3: 计算 ANI_seq
    RAW_ANI_seq = mean(conserved_pairs.local_ani)
    
    // GBRT 只修正这个
    if use_gbrt:
        ANI_seq = gbrt_correct(RAW_ANI_seq, AF_q, AF_r, ...)
    else:
        ANI_seq = simple_debias(RAW_ANI_seq, AF_q, AF_r)

Phase 4: 计算 SDI
    SDI = 1 - (conserved_pairs.len() / total_matched_pairs.len())
    
    // 或用更精细的邻居得分
    for each pair in total_matched_pairs:
        anchor_adjacency = count_neighbor_matches(pair) / neighbor_window_size
    SDI = mean(1 - anchor_adjacency)

Phase 5: 输出
    TSV: query, ref, ani_seq, ani_seq_corrected, sdi, af_q, af_r, ...
```

### TSV 输出格式（建议）

```
query_file  ref_file  ani_seq  ani_seq_corr  sdi  af_q  af_r  n_conserved  n_break  shared_tags  confidence
```

### 何时使用哪个指标？

| 场景 | 关注指标 | 原因 |
|------|---------|------|
| 菌株分型（>99% ANI） | `ani_seq_corrected` | 对标 FastANI，精确 SNP 距离 |
| 物种边界判定（~95%） | `ani_seq_corrected` + `sdi` | 序列 + 结构综合判断 |
| MAG 质量评估 | `sdi` | 高 SDI 暗示组装错误或真实重排 |
| 进化关系推断 | `ani_seq` + `sdi` | 两个维度都重要 |
| 论文对标 FastANI | `ani_seq_corrected` | 可比性验证 |

---

*此文档应与 `ALGORITHM_REFACTOR.md` 配合使用，后者包含更详细的 Rust 代码改造方案。*
