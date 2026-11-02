# linkerdb 数据获取层设计（v0 草稿）

> 范围：本文件只覆盖**数据获取层**——从 PDB 下载结构、解析结构域边界、提取 domain linker、
> 构建非冗余数据集并输出 parquet。下游的特征计算 / 建模不在本会话范围内。

## 1. 数据来源

| 用途 | 来源 | 具体文件 / 接口 |
|---|---|---|
| 结构坐标 | RCSB PDB（mmCIF 优先） | `https://files.rcsb.org/download/{ID}.cif.gz`；镜像回退 PDBe / PDBj；最后回退 PDB 格式 `{ID}.pdb.gz` |
| 结构域划分（主） | CATH v4.4.0 | `cath-domain-boundaries.txt`（CDF 2.0，作者编号 + 插入码）、`cath-domain-boundaries-seqreschopping.txt`（SEQRES 编号）、`cath-domain-list.txt`（CATH 分类码） |
| 结构域划分（回退 1） | SCOPe 2.08 | `dir.cla.scope.2.08-stable.txt` |
| 结构域划分（回退 2） | Pfam（经 PDBe/SIFTS） | `https://www.ebi.ac.uk/pdbe/api/mappings/pfam/{pdb}` |
| 分辨率 / 实验方法 | wwPDB derived data | `resolu.idx`、`pdb_entry_type.txt` |
| 链 → 实体映射 | RCSB GraphQL | `polymer_entities { auth_asym_ids, uniprot_ids }` |
| 去冗余 | RCSB 预计算序列聚类 | `clusters-by-entity-40.txt`（MMseqs2，40% 序列一致性，实体级） |

所有下载都经过同一个下载器：本地缓存（`data/raw/...`）、`.part` 临时文件 + HTTP Range 断点续传、
gzip 完整性校验、原子重命名、失败重试（指数退避）和镜像回退。

## 2. 编号体系

内部所有计算统一用 **SEQRES 索引**（mmCIF 的 `label_seq_id`，1 起始，包含未建模残基）。

* `pdbx_poly_seq_scheme` 给出每个 SEQRES 位置的 `pdb_strand_id`（作者链）、`pdb_seq_num`、
  `pdb_ins_code`（作者编号，**未建模残基也有**）以及是否被建模（`auth_mon_id != '?'`）。
* CATH CDF 文件的作者编号范围 `(链, 起始号, 插入码) – (链, 终止号, 插入码)` 通过
  `(pdb_strand_id, pdb_seq_num, pdb_ins_code) → label_seq_id` 映射到 SEQRES 索引。

## 3. 结构域边界的定义

* 一个 CATH domain 由 1 个或多个连续片段（segment）组成；domain ID 形如 `1cnsA01`，`00` 表示整链单结构域。
* 只保留**同一条链上 ≥ 2 个 domain** 的链（多结构域链）。
* CDF 中的 fragment（`Fxx`，CATH 未分配给任何 domain 的小片段）不属于任何 domain。

## 4. Linker 判定规则

1. 把一条链上所有 domain 的所有 segment 按起点排序。
2. 对每一对**相邻 segment**：若它们属于**不同 domain**，就构成一个 *junction*；
   linker 为两者之间的残基 `end(s_i)+1 … start(s_{i+1})-1`（可能为 0 个残基）。
3. 相邻 segment 属于同一个 domain（domain 内部的未分配缺口）——不是 linker。
4. 链两端未分配的残基——不是 linker（是末端尾巴）。
5. 每个 linker 记录 N 端 domain、C 端 domain、两侧各 5 个 flanking 残基。

## 5. 特殊情况

| 情况 | 处理 |
|---|---|
| 缺失残基（未建模） | 序列取自 SEQRES，记录 `n_missing` 与 `frac_observed` |
| 链断裂 | 相邻已建模残基 C(i)–N(i+1) > 2.0 Å（缺骨架原子时用 CA–CA > 4.2 Å）记为断裂，记录 `n_breaks` |
| 多模型 NMR | 只用第 1 个模型；记录 `n_models` |
| 非标准氨基酸 | 用 gemmi 的残基表映射到母体氨基酸（MSE→M 等），记录 `n_nonstandard` |
| altloc | 取第一个构象 |

## 6. 数据集构建

* **小样本**：从 CATH 多结构域链中随机抽 50 个 PDB entry（固定随机种子），跑通全流程并人工检查。
* **完整集**：
  1. 取所有 PDB 实体，按 RCSB 40% 聚类每簇取分辨率最好的一个代表（类似 PISCES）；
  2. 保留分辨率 ≤ 2.5 Å、且在 CATH 中为多结构域的代表链；
  3. 下载、解析、提取 linker，输出 `linkers.parquet`。

## 7. 边界歧义

* CATH 边界本身有 ±几个残基的不确定性：输出两侧 flanking 序列，供下游自行放宽。
* 非连续 domain（如 `1cnsA01 = 1-87,146-243`）：按 segment 排序后自然产生两个 junction（A→B、B→A）。
* 作者编号在当前 mmCIF 中找不到时，该链丢弃。

## 8. 输出

`linkers.parquet`，每行一个 linker：`pdb_id, chain, domain_n, domain_c, start, end, length, sequence,
n_missing, n_breaks, n_nonstandard, resolution, method, ...`
