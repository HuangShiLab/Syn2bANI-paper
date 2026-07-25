# v8：chain-restricted stratified MLE — 算法与验证记录

> 日期：2026-07-25
> 代码：Syn2bANI `3048177`，子命令 `syn2bani ani`
> 设计文档：代码仓库根目录 `ALGORITHM_MLE.md`（含完整推导）

---

## 1. 为什么放弃 GBRT 校准路线

不是调参问题，是建模问题。三条独立原因：

### 1.1 `raw_ani` 被 mismatch 预算钉死

`local_ani = 1 − hamming / tag_len` 在 `near_match_tolerance = 2`、32 bp tag 下
只能取 `{1.0, 0.969, 0.938}`，所以 `mean(local_ani) ≥ 1 − 2/32 = 0.9375`，
与两个基因组多分歧无关。

实测它对已知真值的**斜率只有 0.336**（85–99.9% ANI 区间）：信号被压缩约 3 倍，
反解时噪声被放大 3 倍。而 GBRT v7 给这个特征的重要性最高（0.3197）。

### 1.2 全基因组 containment 把分歧和共享含量混在一起

```
C_genome ≈ (1 − F) · ANI^k        F = accessory / 非同源比例
```

所以 `C^(1/k)` 返回 `ANI · (1−F)^(1/k)`，偏差大小取决于 `F`，而 `F` 每一对都不同。

**关键推论**：在 `1 + ln(C)/k` 上加一个常数偏移，代数上等价于给所有基因组对
硬编码同一个 `F`。它在共享含量相近的验证集上"有效"，在共享含量变化时必然失效。
这正是 `−0.028` 经验校准在做的事情。

### 1.3 多酶混用单一 `k` 有系统偏差

酶 panel 的 `tag_length` 从 25（BslFI）到 33（CspCI）。98% ANI 时各酶精确匹配
存留率跨 0.51–0.60。混池后用 `k = 32` 反解，98% 处偏高约 0.18 个百分点、
99.5% 处约 0.04 个百分点——绝对值不大，但在菌株分辨率上是决定性的，而且
**换酶组合时会漂移**。

### 1.4 真正的关系

`raw_ani` 和 `mash_ani` 不是两个独立特征，而是**同一个似然的两个矩**：
mismatch 直方图和丢失率由同一个每碱基分歧度生成。拟合这个似然就不需要学习组合，
也不会退化——直方图饱和时，权重自动转移到丢失率上，由 Fisher information
决定，而不是由训练数据决定。

---

## 2. 估计量

Type IIB tag 内部含识别位点，所以 query tag 存在的前提是该位点在 query 中完好。
链内一个 query tag，令 `k` = tag 长度、`s` = 识别位点长度、`b = k − s` = 可变体、
`tol` = mismatch 预算：

```
被找到、body 有 m 个 mismatch    P_m(a) = C(b, m) · (1−a)^m · a^(k−m)    m ≤ tol
未找到                            P_miss(a) = 1 − Σ_{m≤tol} P_m(a)
```

mismatch 只能在 body 上被观测到：位点上的突变直接把 tag 从那个基因组里删除。
所以位点贡献 `a^s`（已折进 `a^(k−m)`）但永远不会表现为观测到的 mismatch。
这也给出一个硬上限——**任何 mismatch 预算都无法保留超过 `a^s` 的 tag**，
95% ANI 时约 0.72，90% 时约 0.50。

对数似然按酶分层求和（每个酶用自己的 `k`、`b`），对单个标量 `a` 最大化。
一维、光滑，golden-section 收敛。

**计数必须来自链内**：这从构造上把 accessory 排除在分母之外，于是 ANI 只度量分歧、
AF 单独度量共享含量。这是对 §1.2 的结构性修复。

**内建质控**：`ani_from_loss`（仅丢失率）和 `ani_from_hist`（仅直方图，按
`P(found)` 重归一化）是同一个量的两个独立估计。差异超过约 5 个标准误说明模型假设
被违反（重复序列、污染、链拉错），标记 `INCONSISTENT`。该检查在期望存留率
< 0.20 时自动关闭：低于该值时存活的直方图只是被截断的尾部，不再是独立估计量。

---

## 3. 验证：精确已知 ground truth

真值由构造给出：对 *E. coli* K-12 施加计数好的替换突变，
所以 true ANI = `1 − n_subs / length`。工具在代码仓库 `prototype/`，
每套仿真约 90 MB（生成物不入库），**不需要 GTDB**。

```bash
cd Syn2bANI/prototype
python3 simulate.py /path/to/e_coli_k12.fasta sim
python3 simulate_accessory.py /path/to/e_coli_k12.fasta simacc 0.95
../target/release/syn2bani ani sim/q_*.fasta sim/ref.fasta --verbose -p
../target/release/syn2bani ani simacc/acc*.fasta simacc/ref.fasta --verbose -p
```

### 3.1 ANI 梯度（12 个基因组，85 → 99.9%，含 400 kb 倒位，无 accessory）

| true ANI | 估计 | 误差 | AF |
|---|---|---|---|
| 85.000 | 85.143 | +0.143 | 0.959 |
| 88.000 | 88.062 | +0.062 | 0.987 |
| 90.000 | 90.027 | +0.027 | 0.995 |
| 92.000 | 92.170 | +0.170 | 0.996 |
| 94.000 | 94.109 | +0.109 | 0.997 |
| 95.000 | 94.986 | −0.014 | 0.998 |
| 96.000 | 96.006 | +0.006 | 0.997 |
| 97.000 | 97.074 | +0.074 | 0.998 |
| 98.000 | 97.982 | −0.018 | 0.998 |
| 99.000 | 98.935 | −0.065 | 0.998 |
| 99.500 | 99.476 | −0.024 | 0.998 |
| 99.900 | 99.885 | −0.015 | 0.998 |

**MAE 0.061%**，零训练数据、零校准模型、零经验偏移。12 对基因组 2.6 秒。

### 3.2 accessory 混淆（6 个基因组，true ANI 固定 95.000，accessory 0 → 50%）

| accessory | 估计 | 误差 | AF | chains |
|---|---|---|---|---|
| 0%  | 95.060 | +0.060 | 1.000 | 1 |
| 10% | 95.077 | +0.077 | 0.898 | 6 |
| 20% | 95.002 | +0.002 | 0.797 | 6 |
| 30% | 95.014 | +0.014 | 0.698 | 6 |
| 40% | 95.058 | +0.058 | 0.597 | 6 |
| 50% | 95.147 | +0.147 | 0.497 | 6 |

**MAE 0.060%。** 估计值保持平坦，AF 精确跟踪 `1 − F`（误差 ≤ 0.003）。
两条信号解耦。作为对照，同一批基因组上全基因组 containment 从 95.18 漂移到
93.27——**这个漂移正是 GBRT 一直在通过 `af_q` 隐式反推的东西**。

### 3.3 与既有数字的关系

⚠️ **不是同一把尺子，不能直接比较。** 下表仅供定位：

| 方法 | MAE | 真值来源 | 数据 |
|---|---|---|---|
| skani | 0.468% | FastANI | 15 对真实口腔/肠道基因组，85–95% |
| Syn2bANI mash_ani | 2.808% | FastANI | 同上 |
| Syn2bANI chained_kmer_ani (v7) | 2.971% (r = −0.107) | FastANI | 同上 |
| Syn2bANI −0.028 经验校准 | 1.031% | FastANI | 同上 |
| Syn2bANI GBRT v7 | 1.142% | FastANI | 同上 |
| **Syn2bANI MLE (v8)** | **0.061%** | **构造真值** | 12 个仿真基因组，85–99.9% |

两点必须在论文里说清楚：

1. v8 的数字来自**仿真**，v7 的数字来自**真实基因组 + FastANI 当真值**。
   在真实数据上重跑 v8 是下一步的首要任务。
2. FastANI 在 < 92% ANI 区间本身可靠性有争议（skani 论文即以此为主要批评点），
   拿它当真值又和 skani 比较，方法学上有偏。真实数据验证应改用
   nucmer/MUMmer（ANIm）或 minimap2 比对得到的真值。

另外，v7 的 `r = −0.107` 不应解释为"block 只覆盖保守区所以系统性偏高"：
**系统偏差会保留相关性**，`r ≈ 0` 说明区间本身算错了。见 §5。

---

## 4. 尚未验证（论文不能声称的部分）

- **没测 indel**。梯度实验只有替换突变加一个 400 kb 倒位；`simulate.py` 支持
  `indel_rate` 但本次跑的是 0，所以 gap 算术路径未被压测。
- **accessory 是 shuffle 造的**：保留碱基组成、彻底破坏同源性，但真实 accessory
  基因 GC 不同，且可能带部分旁系同源。
- **只用了一个物种**（*E. coli* K-12）。GC 含量和重复结构会变。
- **没有真实基因组对**，没有和 FastANI / skani 直接比较。
- **~85% ANI 以下**期望存留率跌破 0.20，一致性交叉检查自动关闭；联合拟合仍然工作，
  但少了一层质控。

---

## 5. v7 代码中发现的缺陷（记录用；`dist` 路径仍然存在这些问题）

审阅 HEAD `5148d88` 时发现。`ani` 子命令是新增路径，不含这些问题。

| 位置 | 缺陷 |
|---|---|
| `ani_calculator.rs::gbrt_debias_ani` | 用 `raw_ani` 伪造 `shared` 和 `containment`（`shared = raw_ani × total.min()`）。因为 `raw_ani ≈ 0.95` 恒定，`containment` 实际是常数，零信息量。真实 shared 数就在旁边未被使用。 |
| `tag_matcher.rs` near-match fallback | 扫描**全部**参考 tag 找全局最小 Hamming，无任何位置约束——O(n²)，而且会配对位置完全无关的 tag。 |
| `synteny_builder.rs::sparse_chain` | 重建出 path 后只保留首尾，再把整段连续切片 `indices[start..=end]` 交给 `create_block`，把 DP 刚拒绝的非共线 anchor 全部重新放行。剩余 anchor 还被逐个发射为单点 block。 |
| `synteny_builder.rs::sparse_chain` | 只有 `gap_diff ≤ INDEL_TOL`，**没有 max_gap**，一条链可横跨整个基因组，于是"chain 内"实际等于全基因组。与上一条叠加，是 `chained_kmer_ani` 得到 r = −0.107 的最可能原因。 |
| `tag_extractor.rs` | tag 从未 canonical 化，倒位段和反向 draft contig（约占任意 draft assembly 的一半）的 tag 完全无法匹配。已由 `canonical_packed` 修复。 |
| v7 `TagSet.sequences` | 为算 `chained_kmer_ani` 而携带原始 contig 序列，使方法**必须有完整装配**才能运行，丧失了在真实 2bRAD 测序数据上的适用性。`ani` 路径只用 tag，保住了这个性质。 |

---

## 6. 论文写作影响

- **核心卖点应改为"闭式统计模型 + 内建质控"**，而不是"更好的校准"。
  卖点是：不需要训练数据、不需要参考数据库、换酶 panel 不用重训、
  每个估计自带标准误和一致性标记。这比"MAE 比 skani 低"更难被质疑。
- **识别位点作为不可消除的截断机制**（`a^s` 上限）是一个新的、可发表的定量结论：
  它解析地给出了 2bRAD 类方法在低 ANI 端的天花板，与 mismatch 预算无关。
- **AF 与 ANI 解耦**可以画成图 3.2 那张表的形式（真值恒定、accessory 变化），
  这是对 containment 类方法最直观的对照实验。
- **GBRT 章节需要重写或降级**为"我们最初尝试的路线及其为何不可行"。
  v2–v7 的模型文件保留在仓库中作为对照。
