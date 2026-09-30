# OAPGC within-species census plan

Generated: `2026-09-30 09:28:35`

## Selection

- annotation rows: **149,921** HQMAGs
- non-redundant strain representatives: **99,215**
- multi-genome SGB clusters: **1,488**
- genomes in structural census: **98,229**
- unique within-species pairs: **30,560,247**

## Cost

- batch size: **250**
- structural tasks: **2,120**
- Syn2b estimate: **405.6** core-h
- digestion estimate: **0.001** core-h
- at 16 local cores: **25.4** node-h
- gzipped structural output estimate: **0.92 GB**

This estimate assumes OAPGC MAG lengths and assembly fragmentation are broadly
similar to HROM. Run the first 20 largest clusters as a smoke test before the
full run; revise constants if measured core-h/pair differs.

## Largest clusters

| SGB | genomes |
|---|---:|
| SGB0467 | 2146 |
| SGB0588 | 2125 |
| SGB1734 | 1661 |
| SGB0325 | 1511 |
| SGB0150 | 1486 |
| SGB0623 | 1480 |
| SGB0583 | 1446 |
| SGB1611 | 1427 |
| SGB0976 | 1409 |
| SGB0692 | 1392 |
| SGB0582 | 1387 |
| SGB1086 | 1368 |
| SGB0558 | 1352 |
| SGB0604 | 1259 |
| SGB1731 | 1196 |
| SGB0640 | 995 |
| SGB0883 | 995 |
| SGB0805 | 908 |
| SGB0009 | 899 |
| SGB0328 | 862 |
