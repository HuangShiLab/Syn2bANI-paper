# ISCB-Asia 2026 — revised 20-minute talk script
**Title:** One coordinate system, three kingdoms: type IIB restriction landmarks from community composition to genome architecture
**Deck:** `ISCB_Asia_2026_Syn2b_TGT_Landmarks.pptx`
**Narrative:** Type IIB landmarks are identity tags and genomic coordinates. One coordinate system moves from community composition to genome architecture.
**Time allocation:** Introduction and Part I ≈ 8 minutes; Part II ≈ 10–11 minutes; closing ≈ 1 minute.

---

## Slide 1 — Title
**Time:** 00:00–00:30

Good morning. Today I will talk about one representation rather than one isolated tool. Type IIB restriction landmarks let us compare bacteria, archaea and fungi in a shared coordinate system. The same data model starts with community composition and ends with genome architecture.

---

## Slide 2 — Opening narrative
**Time:** 00:30–01:30

The premise is simple. Type IIB restriction enzymes cut fixed, reproducible positions across bacterial, archaeal and fungal genomes. The resulting landmarks are homologous. Each one is both an identity tag and a genomic coordinate.

Because identity and position are retained together, we do not need one method for composition, another for strain resolution, and another for structural comparison. We can move from “what is in the community” to “how genomes are arranged” in one coordinate system.

---

## Slide 3 — Part I divider
**Time:** 01:30–01:45

The first part treats landmarks as identity tags: who is present, which strain is present, and what functions may follow.

---

## Slide 4 — Why a shared coordinate system is needed
**Time:** 01:45–02:45

Current representations diverge because each discards something important.

16S and ITS are different molecules, so they cannot align across bacteria and fungi. Shotgun analysis is powerful, but reference databases are often kingdom-partitioned. K-mer sketches are fast, but they discard genomic coordinates.

If a landmark retains both sequence identity and position, one search can be reused across these tasks rather than switching methods by kingdom or question.

---

## Slide 5 — Fast2bRAD-M concept
**Time:** 02:45–03:30

Fast2bRAD-M uses type IIB cut sites to define that shared search space.

Fixed cut sites produce reproducible tags. Because the same tag logic applies across kingdoms, classification becomes a search in one landmark space. The coordinate is retained even though this first layer focuses on identity.

---

## Slide 6 — Fast2bRAD-M benchmark
**Time:** 03:30–05:00

The expanded database contains 394,921 bacterial genomes, 7,777 archaeal genomes and 11,184 fungal genomes.

Relative to 2bRAD-M, Fast2bM reduces runtime by 81 percent and memory by 64 percent while retaining 99.2 percent cross-kingdom accuracy.

This is the identity-tag benchmark. It shows that the shared landmark scan remains accurate while becoming substantially more practical across kingdoms.

---

## Slide 7 — One phylogeny across kingdoms
**Time:** 05:00–05:50

Because the landmarks are homologous, they support one phylogeny across 19,709 taxa spanning bacteria, archaea and fungi.

That matters because beta-diversity can now be evolutionary rather than only compositional. The comparison does not require separate marker systems for 16S, ITS and prokaryotic genomes.

---

## Slide 8 — Strain2bScan key result
**Time:** 05:50–06:50

The strain layer is Strain2bScan.

It converts type IIB landmarks into strain-cluster-resolved profiles for native and conventional data. For novel strains, it recovers the correct closest relative in 83 percent of cases.

This is the step where species-level profiles often stop. Strain resolution lets functional interpretation follow the strain rather than only the species label.

---

## Slide 9 — Strain2bScan evidence
**Time:** 06:50–07:50

The supporting evidence has two parts.

Reduced-representation markers preserve strain relationships better than 16S. And the profiler remains computationally lightweight. Together, these results show that the strain layer is biologically informative and scalable.

---

## Slide 10 — Meta2bAnalyst
**Time:** 07:50–08:20

Meta2bAnalyst closes the identity-tag loop.

It connects landmark profiles to downstream statistics, ecological visualization and AI-agent-assisted exploration. I mention it briefly here; it provides the bridge from profiling results to interpretable ecological analysis.

---

## Slide 11 — Part II divider
**Time:** 08:20–08:40

The second part is the core of the talk.

We now stop treating the landmark only as an identity tag and start using its position.

---

## Slide 12 — From tag identity to tag order
**Time:** 08:40–09:40

Here is the hypothesis.

Represent each genome as a set of fixed type IIB anchors. The order of those anchors is a compact macrogenomic signature.

If two genomes share landmark order, they share architecture. If landmark order breaks, those breaks point to rearrangement. This is the property that ANI alone cannot provide.

---

## Slide 13 — TGT structure
**Time:** 09:40–10:20

The TGT table is the structural bridge.

Each record contains the canonical tag and its coordinate. The tag defines the landmark identity. The coordinate preserves order and enables breakpoint detection.

That is why this is not just a classification hash. It is a genomic anchor.

---

## Slide 14 — Syn2b / Syn2bANI method
**Time:** 10:20–11:30

Syn2b and Syn2bANI compare landmark order directly.

Shared anchors are chained collinearly. Tags inside chains estimate ANI. Anchor order and missing anchors report synteny, inversions, translocations and indels.

There is no alignment-first step. The coordinate-carrying anchors provide both sequence divergence and architectural information.

---

## Slide 15 — Structural benchmark
**Time:** 11:30–13:30

This is the key structural benchmark.

For E. coli K-12 versus W3110, Syn2b localizes the approximately 780 kilobase rrn-mediated inversion within 5.5 kilobases. I want to emphasize that this limit comes from the repeats at the inversion junctions; it is not method error.

For indels, Syn2b calls 184 of 201 events correctly, with a median size ratio of exactly 1.000. The landmark-based structural comparison takes about 50 milliseconds, compared with 8.2 to 9.6 seconds for dnadiff.

---

## Slide 16 — Contrast with ANI-only sketches
**Time:** 13:30–14:20

This leads to the central contrast.

K-mer sketches can estimate ANI quickly, but they discard coordinates. A faster sketch still cannot localize an inversion.

Architecture may also be nearly invisible at the ANI level. Syn2bANI fills that missing layer by retaining landmark order.

---

## Slide 17 — ANI accuracy at GTDB-R207 scale
**Time:** 14:20–15:30

The coordinate layer is added without losing ANI performance.

On the GTDB-R207 benchmark, Syn2bANI was evaluated on 43,334 genome-disjoint pairs, with calibrated calls on 39,903. It reached 0.619 ANI points MAE and showed a genome-level bootstrap improvement over skani.

---

## Slide 18 — Catalog-scale search and SV follow-up
**Time:** 15:30–16:20

At catalog scale, skani remains the fastest ANI-only search. Syn2bANI is slower but returns coordinate-aware architecture in the same framework.

For top-hit structural follow-up, the landmark route is substantially lighter than full alignment.

---

## Slide 19 — B. longum abfA
**Time:** 16:20–17:00

The first biological example is Bifidobacterium longum.

Although pairwise ANI is high, local landmark coverage resolves abfA-cluster state. Species-level identity therefore does not guarantee conservation of a functionally important locus.

---

## Slide 20 — H. pylori cagPAI
**Time:** 17:00–18:00

The second example is the Helicobacter pylori cagPAI cohort.

The coordinate layer reveals island states and lineage-associated architecture, while disease associations remain lineage-confounded. These are exactly the distinctions that ANI alone cannot provide.

---

## Slide 21 — Census scale
**Time:** 18:00–18:50

The coordinate layer also extends to catalogs.

Across GTDB, breakpoint prevalence declines in the very-high-ANI range. This is a descriptive census result. It does not infer migration direction or mutation rate.

---

## Slide 22 — How to use the two layers
**Time:** 18:50–19:30

The practical message is to match the tool to the question.

Use Fast2bM or Meta2bSeek for cross-kingdom community profiling. Use Strain2bScan for strain-cluster and functional questions. Use Syn2bANI when both ANI and genome architecture are needed. For MAGs, mid-ANI pairs and flagged comparisons, use raw gated estimates. Full alignment remains appropriate for base-level confirmation.

---

## Slide 23 — Conclusion
**Time:** 19:30–20:00

To conclude, return to the opening idea.

Type IIB landmarks are identity tags and genomic coordinates. That single property lets one assay move from who is present to how genomes are arranged.

Thank you.

---

## Timing checkpoint

- Slide 1: 0:00
- Slide 2: 0:30
- Slide 3: 1:30
- Slide 4: 1:45
- Slide 5: 2:45
- Slide 6: 3:30
- Slide 7: 5:00
- Slide 8: 5:50
- Slide 9: 6:50
- Slide 10: 7:50
- Slide 11: 8:20
- Slide 12: 8:40
- Slide 13: 9:40
- Slide 14: 10:20
- Slide 15: 11:30
- Slide 16: 13:30
- Slide 17: 14:20
- Slide 18: 15:30
- Slide 19: 16:20
- Slide 20: 17:30
- Slide 21: 18:20
- Slide 22: 19:10
- Finish: 20:00

## If running short

- Compress slide 4 to two sentences.
- Present slide 5 as a transition and move directly to slide 6.
- Limit slide 17 to the final MAE comparison.
- Present slide 18 only as the practical operating-point slide.
