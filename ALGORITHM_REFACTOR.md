# Syn2bANI 核心算法梳理与"三部分"输出架构

> **目标**：理清当前算法的计算流程，设计"ANI_seq + SDI + Composite"的三指标输出体系，明确 GBRT 的修正边界。

---

## 一、当前算法流程（逐步拆解）

### Step 1: Digestion（酶切提取标签）

```
Input:  genome sequence (Vec<u8>), enzyme config
Output: Vec<GenomeTag>

for each recognition site in enzyme.pattern:
    extract flanking sequence (tag)
    verify purity (only A/T/C/G)
    create GenomeTag {
        position: absolute_pos,
        sequence: [u8; 32],       // 实际序列
        packed_sequence: u64,      // 2-bit packed (SIMD加速)
        seq_len: actual_length,    // BcgI=32, PpiI=27, etc.
        direction: '+' or '-',
        enzyme: "BcgI",
    }
```

**当前局限**：
- 只提取了序列和位置，没有记录"上下文"
- 没有区分"核心区标签"和"边缘区标签"

---

### Step 2: Matching（标签匹配）

```
Input:  Vec<GenomeTag> (query), Vec<GenomeTag> (reference)
Output: Vec<MatchedPair>

for each q_tag in query_tags:
    // 在 reference 中寻找候选
    candidates = r_tags.filter(|r| {
        r.position in [q.position - window, q.position + window]
        && r.enzyme == q.enzyme
    })
    
    for each candidate in candidates:
        // 序列相似度
        similarity = hamming_similarity(q.packed, r.packed, q.seq_len)
        
        if similarity > threshold (e.g., 0.85):
            create MatchedPair {
                q_tag, r_tag,
                similarity,  // ← 目前唯一的量化信息
            }
```

**当前局限**：
- `MatchedPair` 只记录了 `similarity`（标量）
- 没有利用"邻居标签是否也匹配"这一关键结构信息
- 所有 matched pairs 被平等对待（uniform weighting）

---

### Step 3: ANI Calculation（当前实现）

```
Input:  Vec<MatchedPair>, AniConfig
Output: AniResult

// 当前：简单平均
let total_similarity = matched_pairs.iter().map(|p| p.similarity).sum();
let ani = total_similarity / matched_pairs.len() as f64;

// 可选：GBRT 修正（当前做法）
if config.use_gbrt_debias {
    let features = extract_features(matched_pairs, q_tag_set, r_tag_set);
    let bias = gbrt_model.predict(features);
    ani += bias;  // ← 问题：修正的是"总体ANI"
}
```

**当前问题**：
- `similarity` 是混合信号：既包含 SNP 信息，也可能包含 small indel
- GBRT 修正的是最终输出，但修正目标不明确（FastANI？MUMmer？）
- 没有输出"为什么 ANI 是这个值"的可解释性

---

## 二、核心问题：MatchedPair 到底携带了什么信息？

每个 `MatchedPair` 在概念上包含**两个独立维度**的信息：

### 维度 A：序列一致性（Sequence Identity）
- **来源**：`hamming_similarity(q.sequence, r.sequence)`
- **含义**：两个标签的 DNA 序列有多相似
- **对应生物学**：SNP 密度、small indel（< tag_length）
- **对标**：FastANI 的局部 alignment identity

### 维度 B：结构保守性（Synteny Conservation）
- **来源**：`q_tag` 的邻居标签是否也在 `r_tag` 附近有匹配
- **含义**：这个匹配是"孤立的"还是"在保守区块中"
- **对应生物学**：
  - 邻居都匹配 → 该区域 synteny 保守
  - 邻居不匹配 → 可能是重排断点、大 indel、或假阳性匹配
- **对标**：无（FastANI/skani 完全不输出此信息）

**关键洞察**：维度 A 和维度 B 是正交的：
- 高序列相似 + 高结构保守 → 真正的同源保守区
- 高序列相似 + 低结构保守 → 可能是重复序列、HGT、或假阳性
- 低序列相似 + 高结构保守 → 高 SNP 但基因组骨架保守
- 低序列相似 + 低结构保守 → 非同源区域

---

## 三、改进后的"三部分"输出架构

### 3.1 新的数据结构

```rust
// 改进的 MatchedPair：显式分离两个维度
pub struct MatchedPair {
    pub q_tag: GenomeTag,
    pub r_tag: GenomeTag,
    
    // 维度 A：序列层面
    pub sequence_identity: f64,      // 纯 Hamming similarity
    pub n_mismatches: u32,           // 明确的 mismatch 数
    
    // 维度 B：结构层面
    pub synteny_score: f64,          // [0, 1], 邻居匹配比例
    pub is_conserved_block: bool,    // synteny_score > threshold
    pub block_size: usize,           // 连续保守标签数
}

// 改进的 AniResult：三部分输出
pub struct AniResult {
    // 第一部分：ANI_seq（对标 FastANI）
    pub ani_seq: f64,                // 仅序列身份的平均
    pub ani_seq_corrected: f64,      // GBRT 修正后的序列 ANI
    pub n_seq_pairs: usize,          // 参与计算的 matched pairs
    
    // 第二部分：SDI（Syn2bANI 独有）
    pub sdi: f64,                    // Structural Divergence Index [0, 1]
    pub n_conserved_blocks: usize,   // 保守区块数
    pub n_breakpoints: usize,        // 结构断点数
    
    // 第三部分：Composite（综合指标）
    pub syn2bani_composite: f64,     // 加权组合
    pub weight_seq: f64,             // 序列权重
    pub weight_struct: f64,          // 结构权重
}
```

### 3.2 计算流程（改进后）

```rust
impl TagMatcher {
    pub fn match_with_structure(q: &TagSet, r: &TagSet, config: &MatchConfig) -> MatchResult {
        let mut pairs = Vec::new();
        
        // Phase 1: 基础匹配（同当前）
        for q_tag in &q.tags {
            for r_tag in find_candidates(q_tag, r) {
                let seq_id = sequence_identity(q_tag, r_tag);
                if seq_id > config.min_identity {
                    pairs.push((q_tag, r_tag, seq_id));
                }
            }
        }
        
        // Phase 2: 结构分析（新增）
        let mut enriched_pairs = Vec::new();
        for (q_tag, r_tag, seq_id) in &pairs {
            // 检查 q_tag 的邻居在 reference 中是否也有匹配
            let neighbors_q = get_neighbor_tags(q_tag, &q.tags, window=3);
            let neighbors_r = get_neighbor_tags(r_tag, &r.tags, window=3);
            
            let mut synteny_matches = 0;
            for nq in neighbors_q {
                if let Some(nr) = find_match_in_pairs(nq, &pairs) {
                    if is_position_consistent(nq, nr, q_tag, r_tag) {
                        synteny_matches += 1;
                    }
                }
            }
            
            let synteny_score = synteny_matches as f64 / neighbors_q.len().max(1) as f64;
            
            enriched_pairs.push(MatchedPair {
                q_tag: *q_tag,
                r_tag: *r_tag,
                sequence_identity: *seq_id,
                n_mismatches: count_mismatches(q_tag, r_tag),
                synteny_score,
                is_conserved_block: synteny_score > 0.6,
                block_size: estimate_block_size(q_tag, &pairs),
            });
        }
        
        MatchResult { pairs: enriched_pairs }
    }
}

impl AniCalculator {
    pub fn calculate_three_metrics(result: &MatchResult, config: &AniConfig) -> AniResult {
        let pairs = &result.matched_pairs;
        
        // ===== 第一部分：ANI_seq =====
        let ani_seq_raw = pairs.iter()
            .map(|p| p.sequence_identity)
            .sum::<f64>() / pairs.len() as f64;
        
        // GBRT 仅修正 ANI_seq
        let ani_seq_corrected = if config.use_gbrt_debias {
            let features = extract_gbrt_features(result);
            let bias = gbrt_model.predict(features);
            (ani_seq_raw + bias).clamp(0.0, 1.0)
        } else {
            ani_seq_raw
        };
        
        // ===== 第二部分：SDI =====
        let n_conserved = pairs.iter().filter(|p| p.is_conserved_block).count();
        let sdi = 1.0 - (n_conserved as f64 / pairs.len() as f64);
        
        let n_breakpoints = count_breakpoints(pairs);
        
        // ===== 第三部分：Composite =====
        // 加权组合：序列 + 结构
        // 权重可配置：高序列相似但高 SDI → 可能是假同源
        let w_seq = config.weight_sequence;      // e.g., 0.7
        let w_struct = config.weight_structure;   // e.g., 0.3
        
        // 结构惩罚：SDI 高 → composite ANI 降低
        let structural_penalty = sdi * w_struct;
        let syn2bani_composite = ani_seq_corrected * w_seq 
                               + ani_seq_corrected * (1.0 - sdi) * w_struct;
        // 简化为：composite = ani_seq * (1 - sdi * w_struct)
        
        AniResult {
            ani_seq: ani_seq_raw,
            ani_seq_corrected,
            n_seq_pairs: pairs.len(),
            sdi,
            n_conserved_blocks: n_conserved,
            n_breakpoints,
            syn2bani_composite,
            weight_seq: w_seq,
            weight_struct: w_struct,
        }
    }
}
```

---

## 四、GBRT 修正边界（关键决策）

### 4.1 GBRT 修正什么？

**只修正 ANI_seq 的技术性偏差**：

| 偏差来源 | 特征 | 是否修正 | 原因 |
|----------|------|---------|------|
| GC 含量影响酶切效率 | 高 GC 基因组标签数偏少 | ✅ | 技术性采样偏差 |
| 标签长度差异 | BcgI(32bp) vs PpiI(27bp) 的 identity 基准不同 | ✅ | 技术性测量偏差 |
| 稀疏采样方差 | 标签数 < 100 时估计不稳定 | ✅ | 统计技术性偏差 |
| 酶识别位点分布 | 某些基因组中特定酶位点过少 | ✅ | 可通过多酶 panel 补偿 |

### 4.2 GBRT 不修正什么？

**不应修正定义性差异**：

| 差异来源 | 是否修正 | 原因 |
----------|---------|------|
| 稀疏采样 vs 密集采样 | ❌ | 这是估计量的固有特性 |
| 固定锚点 vs 随机 k-mer | ❌ | 采样策略不同 |
| 结构信息（SDI） | ❌ | 这是 Syn2bANI 独有的价值 |
| 与 FastANI 的系统性差异 | ❌ | 两个估计量估计同一参数但方式不同 |

### 4.3 Ground Truth 选择（重新明确）

**GBRT 训练的 y 不再是 "FastANI - Syn2bANI"，而是 "技术性偏差"**

```python
# 修正后的 GBRT 目标
y = ani_mummer - ani_seq_raw  # MUMmer ANI 作为金标准
# 或
y = ani_mlst - ani_seq_raw    # MLST 距离转换的 ANI

# 如果 MUMmer/MLST 不可得，用 leave-one-species-out 验证：
# 同一物种内完整基因组对的 Syn2bANI 自身一致性
```

---

## 五、论文中的表述框架

### 5.1 三指标的定义

```markdown
**ANI_seq**: Sequence-level average nucleotide identity estimated from 
tag sequence matches. Directly comparable to FastANI and skani, 
with technical biases (GC-dependent cleavage, tag-length effects) 
corrected by GBRT.

**SDI (Structural Divergence Index)**: Fraction of matched tags located 
in synteny-conserved blocks. SDI = 0 indicates perfect collinearity; 
SDI > 0.2 suggests substantial rearrangement or large indels. 
This metric is unique to anchor-based methods and invisible to 
sequence-only approaches.

**Syn2bANI_composite**: A weighted integration of ANI_seq and (1-SDI), 
reflecting that high sequence similarity in rearranged regions may 
not indicate true overall genome relatedness. Default weights: 
0.7×ANI_seq + 0.3×(1-SDI).
```

### 5.2 与 skani/FastANI 的对比表（更新）

| 能力 | FastANI | skani | Syn2bANI |
|------|---------|-------|----------|
| ANI 估计 | ✅ | ✅ | ✅ (ANI_seq) |
| 碎片化鲁棒 | ⚠️ | ⚠️ | ✅ |
| 结构分歧检测 | ❌ | ❌ | ✅ (SDI) |
| 湿实验验证 | ❌ | ❌ | ✅ (2bRAD-M) |
| 可解释输出 | ❌ (单一数字) | ❌ (单一数字) | ✅ (三部分) |

### 5.3 典型使用场景

**场景 1：菌株分型（需要精确的 SNP 距离）**
```bash
# 用户关注：ANI_seq（修正后）
syn2bani dist strainA.fna strainB.fna --output-seq-only
# 输出：ANI_seq_corrected = 99.2%（对标 FastANI）
```

**场景 2：基因组比较（关注整体进化关系）**
```bash
# 用户关注：Composite（综合）
syn2bani dist magA.fna refB.fna
# 输出：
#   ANI_seq = 98.5%
#   SDI = 0.15（15% 的标签在重排区域）
#   Composite = 97.8%（比 ANI_seq 低 0.7%，因为结构分歧）
```

**场景 3：SV 检测（关注结构变异）**
```bash
syn2bani struct magA.fna refB.fna --paf --rearrangement
# 输出：SDI = 0.35，断点列表，PAF 文件
```

---

## 六、代码实现优先级

### Phase A（立即）：数据结构改造
1. 修改 `MatchedPair`，增加 `synteny_score` 和 `is_conserved_block`
2. 修改 `MatchResult`，增加结构分析阶段
3. 修改 `AniResult`，改为三部分输出
4. 更新 TSV 输出格式（增加列）

### Phase B（接下来）：GBRT 修正范围收紧
1. 修改 `AniCalculator`，GBRT 只修正 `ani_seq`
2. 重新训练 GBRT（目标改为 MUMmer/MLST，而非 FastANI）
3. 如果没有 MUMmer 数据，用同物种完整基因组自对比作为近似 ground truth

### Phase C（随后）：Composite 权重优化
1. 在 GTDB-R207 数据集上测试不同权重组合
2. 找到使 Composite 与系统发育树最一致的权重
3. 默认权重设为 `w_seq=0.7, w_struct=0.3`

---

## 七、关键决策点（需要用户确认）

1. **SDI 的计算方式**：
   - 选项 A：邻居匹配比例（当前建议）
   - 选项 B：基于最大保守区块的比率
   - 选项 C：与参考基因组的共线性得分（需全局比对）

2. **Composite 的加权策略**：
   - 选项 A：线性加权 `w1*ANI_seq + w2*(1-SDI)`
   - 选项 B：乘法惩罚 `ANI_seq * (1 - SDI * penalty)`
   - 选项 C：不输出 Composite，只输出 ANI_seq + SDI，让用户自行判断

3. **GBRT 的 ground truth**：
   - 选项 A：MUMmer/NUCmer（最精确但数据量受限）
   - 选项 B：同物种完整基因组对的 Syn2bANI 自一致性（数据充足）
   - 选项 C：FastANI（数据充足，但承认 FastANI 也有偏差）

我的建议：
- SDI：选项 A（邻居匹配比例，计算简单，直觉清晰）
- Composite：选项 C（先不输出，让用户/审稿人判断是否需要）
- GBRT：选项 B（用同物种完整基因组作为 self-consistency ground truth）

你觉得这个方向对吗？需要我先把 Phase A 的数据结构改造写出来吗？
