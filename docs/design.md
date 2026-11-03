# linkerdb 设计

> 范围：从 PDB 下载结构、解析结构域边界、提取 domain linker，并构建一个面向**融合蛋白 / 生物药设计**的
> linker 数据库（第 13 节）。下游的建模不在范围内。
>
> 版本记录：v0 初稿（git 历史）；v1 自审修订（第 10 节）；v2 按 50 结构小样本人工检查改规则（第 11 节），
> 并跑通 ≤ 2.5 Å、40% 去冗余的完整提取（第 11.5 节）；**v3 把产物从"单一非冗余集"改为"两层数据库"（第 13 节）**,
> 目标用户是做融合蛋白 / 双抗 / Fc 融合的工程人员。第 12 节是已知局限（已随 v3 更新）。
>
> §1–§7 是提取规则，§8–§11 是第一阶段的非冗余数据集构建和验证，§12 是已知局限，§13 是 v3 的数据库层。

---

## 1. 数据来源

| 用途 | 来源 | 具体文件 / 接口 |
|---|---|---|
| 结构坐标 | RCSB PDB（**mmCIF 优先**） | `https://files.rcsb.org/download/{ID}.cif.gz`；依次回退 wwPDB / PDBe / PDBj 镜像；全部失败才回退 PDB 格式 `{ID}.pdb.gz` |
| 结构域划分（主） | CATH v4.4.0 | `cath-domain-boundaries.txt`（CDF 2.0，**作者编号 + 插入码**，含 fragment）、`cath-domain-boundaries-seqreschopping.txt`（**SEQRES 编号**）、`cath-domain-list.txt`（CATH 分类码、domain 长度） |
| 结构域划分（回退 1） | SCOPe 2.08 | `dir.cla.scope.2.08-stable.txt` |
| 结构域划分（回退 2） | Pfam（经 PDBe / SIFTS） | `https://www.ebi.ac.uk/pdbe/api/mappings/pfam/{pdb}` |
| 分辨率 / 实验方法 | wwPDB derived data | `resolu.idx`、`pdb_entry_type.txt` |
| 链 → 实体、UniProt | RCSB GraphQL | `polymer_entities { auth_asym_ids, uniprot_ids }`（批量查询 + 本地缓存） |
| 去冗余 | RCSB 预计算序列聚类 | `clusters-by-entity-40.txt`（MMseqs2，40% 序列一致性，实体级） |
| 天然性判定 / UniProt 映射 | SIFTS（主）+ mmCIF 自身（回退） | `pdb_chain_uniprot.tsv.gz`（每周更新）；`_struct_ref`、`_struct_ref_seq`、`_struct_ref_seq_dif` |

**下载器**（所有网络数据共用）：

* 本地缓存：`data/raw/<类别>/...`，PDB 结构按 `mmcif/<中间两位>/<id>.cif.gz` 分目录。已存在且校验通过的文件直接复用。
* 断点续传，分两层：
  * **字节级**：写入 `<file>.part`，中断后重新请求时带 `Range: bytes=<已下载大小>-`。服务器返回 206 时续写，
    返回 200（不支持 Range）时从头写。完成后先校验，再原子 `rename` 成正式文件。
  * **任务级**：每个 entry 的处理结果写成单独的 checkpoint（`data/work/<run>/entries/<id>.json`）。
    重跑时跳过已完成的 entry，只处理剩下的和之前失败的。
* 校验：gzip 能完整解压、mmCIF 以 `data_` 开头。校验失败的文件删除后重下。
* 重试：网络错误或 5xx 用指数退避（1, 2, 4, 8 s），404 不重试，直接换镜像。

---

## 2. 编号体系（核心约定）

内部**所有**残基位置统一使用 **SEQRES 索引**（mmCIF `label_seq_id`，1 起始，包含未建模残基）。

* `_pdbx_poly_seq_scheme` 给出每个 SEQRES 位置的作者链 `pdb_strand_id`、作者编号 `pdb_seq_num`
  和插入码 `pdb_ins_code`（**未建模残基也有作者编号**），以及该残基是否被建模（`auth_mon_id != '?'`）。
* CATH CDF 的作者编号范围通过 `(pdb_strand_id, pdb_seq_num, pdb_ins_code) → label_seq_id` 映射到 SEQRES 索引；
  seqreschopping 文件本身就是 SEQRES 索引。两者**互相校验**，规则见 §7。
* 微异质性（同一个 `seq_id` 有多行，`hetero = y`）：取第一行，并记录 `microheterogeneity` 标志。
* 回退到 PDB 格式时没有 `pdbx_poly_seq_scheme`，改用 gemmi 把已建模残基比对到 SEQRES，并标记 `format = pdb`。

---

## 3. 结构域边界的定义

* 一个 CATH domain 由 1 个或多个连续片段（segment）组成。domain ID 形如 `1cnsA01`，`00` 表示整条链就是一个 domain。
* 只保留**同一条链上 ≥ 2 个 domain** 的链，即多结构域链。
* CDF 中的 fragment（`Fxx`）是 CATH 未分配给任何 domain 的小片段，不属于任何 domain。
  它们夹在两个 domain 之间时，自然成为 gap linker 的一部分；落在 domain 内部或链末端时不算 linker。
* 每个 domain 带上 CATH 分类码（C.A.T.H），来自 `cath-domain-list.txt`，用于后续去冗余和统计。

---

## 4. Linker 判定规则

### 4.1 Junction（结构域交界）

1. 把一条链上所有 domain 的所有 segment 按 SEQRES 起点排序。
2. 两个**相邻 segment 属于不同 domain** 时，构成一个 junction：N 侧 domain 的 segment 末端 `e`，C 侧 domain 的 segment 起点 `s`。
3. 相邻 segment 属于**同一个** domain（domain 内部缺口）不构成 junction。
4. 链两端未分配的残基是末端尾巴，不构成 junction。
5. Junction 分类，按两个 domain 的**序列跨度**（第一个 segment 起点到最后一个 segment 终点）判断：
   * `sequential`：两个跨度不重叠（domain 内部有缺口不影响）；
   * `insertion`：一个 domain 的跨度完全包含另一个（如 `1cnsA01 = 1-87,146-243` 中间插着 `1cnsA02 = 88-145`），
     会产生 A→B 和 B→A 两个 junction；
   * `interleaved`：跨度有重叠但不是包含关系（多个 domain 互相穿插）。
6. **无二级结构的小片段不参与 junction 判定**：一个 segment 内如果一个规则二级结构残基都没有（纯 coil），
   就视为未分配残基（同 CATH fragment），然后再找 junction。原因见 §11（例：CATH 把激酶活化环中的 8 个残基
   划给了 N 叶，会在 C 叶内部产生两个互相重叠的伪 linker）。这类链打 `minor_segment_dropped` 标志。

### 4.2 两种 linker 定义（同时输出）

**数据事实**（CATH v4.4.0，全部 358,437 个 junction）：**86.0% 的 junction 两侧 domain 直接相接（gap = 0）**，
gap 为 1–10 的占 9.4%，≥ 30 的占 1.3%。CATH 会把绝大多数残基都分配给某个 domain，
所以按字面意思取"两个 domain 之间的残基"，大多数 linker 都是空的。为此同时输出两种定义：

1. **gap linker（严格定义）**：`e+1 … s-1`，即 CATH 未分配给任何 domain 的残基，长度 ≥ 0。
2. **结构 linker（主定义，`ss_*` 列）**：从 junction 出发，向 N 侧和 C 侧各自延伸，经过的残基
   都是**非规则二级结构**，遇到第一个规则二级结构残基就停。
   直观理解：结构 linker 是 N 侧 domain 最后一个二级结构元件和 C 侧 domain 第一个二级结构元件之间的那段。
   * 二级结构由内置的 DSSP 实现（Kabsch–Sander 氢键能量，含 β-bulge 合并）根据第 1 个模型的骨架坐标计算。
     在测试结构上与 mkdssp 4.2.2 的三态结果完全一致（§11.2）。
   * **规则二级结构 = DSSP 的 H（α）、I（π）、E（β 链）**。G（3₁₀）、B（孤立桥）、T、S 和 coil 都算非规则：
     linker 中常见单圈 3₁₀ 和孤立桥，以它们为界会把 linker 切碎。未建模残基没有坐标，DSSP 记为 `?`，按非规则处理。
   * 延伸只能在两侧 segment 内部进行，不跨进第三个 domain。每侧最多延伸 `MAX_EXT = 15` 个残基
     （由小样本标定：不设上限时单侧延伸的 P95 为 15，见 §11），到达上限打 `capped_n` / `capped_c` 标志。
   * gap 内部的残基无论二级结构如何都算 linker。
   * 结构 linker 长度为 0 时，用 `junction_ss` 说明原因：
     `helix_spanning`：junction 两侧残基在同一段连续螺旋里（"螺旋型 linker"，如天冬氨酸转氨酶 C 端长螺旋）；
     `strand_spanning`：在同一条 β 链里（如 α-裂解蛋白酶两个 β 桶之间共享的链）；
     `sse_abutting`：一侧螺旋结束，另一侧链立即开始。
     这几类没有不依赖人为阈值的 linker 边界，因此不强行划定。长度 > 0 时 `junction_ss = 'coil'`。

### 4.3 每个 linker 记录的内容

位置（SEQRES 和作者编号）、长度、规范一字母序列、原始三字母残基、两侧各 10 个 flanking 残基及其二级结构、
DSSP 二级结构字符串、缺失残基数、链断裂、非标准残基、UniProt 区间、两侧 domain 的 CATH 分类码等。
完整 schema 见 §9。

---

## 5. 特殊情况

| 情况 | 处理 |
|---|---|
| 缺失残基（SEQRES 中有、坐标中没有） | 序列始终取自 SEQRES，记录 `n_missing` / `frac_observed`。未建模残基的 DSSP 状态记为 `?`（按非规则处理）。linker 全部未建模时打 `fully_disordered` 标志，保留（柔性 linker 正是研究对象），由下游决定是否使用 |
| 链断裂 | 相邻 SEQRES 位置都已建模，但 `C(i)–N(i+1) > 2.0 Å`（缺骨架原子时改用 `CA–CA > 4.2 Å`）时记为断裂，统计 linker 区间及两端锚定残基内的 `n_breaks`，打 `chain_break` 标志。未建模残基造成的缺口由 `n_missing` 描述。DSSP 在断裂处（C–N > 2.5 Å）不放置酰胺 H，也不跨断点判定 turn/bridge |
| 多模型 NMR | 坐标、建模状态、二级结构都只用**第 1 个模型**；记录 `n_models`。全量数据集按分辨率筛选时 NMR 自然被排除；小样本里保留 NMR，用来验证流程 |
| 非标准氨基酸 | 一字母映射分层：① wwPDB 的 `entity_poly.pdbx_seq_one_letter_code_can`（逐位置、长度与 SEQRES 一致时才用）；② gemmi 残基表的母体氨基酸；③ `_pdbx_struct_mod_residue.parent_comp_id`；④ 仍无法确定时记为 `X`。同时保留原始三字母码，并记录 `n_nonstandard` 及其列表（例如 MSE 是硒代标记，SEP 是磷酸化，下游可以区分） |
| altloc | gemmi 只保留第一个构象（`remove_alternative_conformations`） |
| 非 L 型多肽链 | `entity_poly.type != 'polypeptide(L)'` 的链直接排除 |
| CATH 链号 `0` | 早期无链号条目在重整后的 mmCIF 中通常是 `A`。只有当该 entry 恰好只有一条蛋白链时才接受这种映射，并打标志 |
| 负编号 / 插入码 | CDF 解析时完整保留，映射时按 `(编号, 插入码)` 精确匹配 |

---

## 6. 天然性过滤（"天然多结构域蛋白"）

PDB 里有大量**工程化融合构建体**，比如 GPCR–T4 溶菌酶融合体。CATH 会把它们正常拆分成多结构域：
`2rh1A`（β2AR–T4L，2.4 Å）被拆成 `2rh1A01 = 36-237,398-477`（GPCR）和 `2rh1A02 = 238-397`（T4L）。
这类 linker 是人工设计的，必须排除。判定基于残基到 UniProt 的映射：优先用 **SIFTS**（`pdb_chain_uniprot.tsv`，PDBe 每周更新，SEQRES 区间 ↔ UniProt 区间），
SIFTS 没有收录的链回退到 mmCIF 的 `_struct_ref_seq`。之所以 SIFTS 优先，是因为有些天然蛋白在 mmCIF 里只引用了
`PDB` 自身（小样本中的 `7kw9`），SIFTS 却给出了 UniProt 映射。

* **链级（不通过则整条链排除）**
  * 该链映射到的 UniProt accession **恰好 1 个**。
    0 个记为 `no_uniprot`（从头设计蛋白、合成构建体，以及没有 UniProt 条目的抗体/TCR）；
    ≥ 2 个记为 `chimera`（融合体或嵌合体）。
* **linker 级（打标志，最终集排除）**
  * linker（含两侧各 1 个锚定残基）中每个残基都能映射到 UniProt，
    并且 UniProt 编号连续（SEQRES 位置差 = UniProt 位置差），否则记为 `unp_discontinuous` 或 `unp_unmapped`。
    这样可以抓到"把天然区段删掉、换成 GS linker"这类构建体。
  * linker 区间内的 `_struct_ref_seq_dif` 如果是 linker、insertion、expression tag、cloning artifact、deletion 之类，
    记为 `engineered_linker`；点突变（engineered mutation、conflict、variant）只记录数量，不排除。

---

## 7. 边界歧义的处理

1. **CATH 边界本身的不确定性（±几个残基）**：同时输出 gap linker 和结构 linker，并提供两侧各 10 个 flanking
   残基及其二级结构。下游可以按需重新划定边界，不需要重新下载数据。
2. **零长度 gap**：用结构 linker 解决；仍为 0 时说明 junction 落在连续二级结构中，用 `junction_ss` 标明。
3. **不连续 domain / 插入型 domain**：`junction_type = insertion`（或 `interleaved`），并记录 junction 序号和 domain 对。
   同一对 domain 可能有两个 junction，下游统计时要注意。没有二级结构的小片段先剔除（§4.1 第 6 条）。
4. **编号不一致**：CDF（作者编号）映射到 SEQRES 后与 seqreschopping 逐个比对 segment 端点：
   * 全部一致 → `numbering_check = agree`；
   * 不一致时，利用"CATH 只按坐标切分，所以正确映射的每个 segment 端点都必须是已建模残基"这一事实裁决
     （小样本中一致的 57 条链 312/312 个端点都已建模）：只有一方的端点全部已建模 → 采用该方
     （`resolved_cdf` / `resolved_seqres`，打 `numbering_resolved` 标志）；两方都满足或都不满足 → `numbering_disagree`，整链排除；
   * 只有其中一个可用 → `cdf_only` / `seqres_only`；
   * 对 CATH，最终采用的映射如果有端点落在未建模残基上，整链按 `mapping_failed` 排除。
5. **domain 长度校验**：映射后每个 domain 内已建模残基数与 CATH 列表（含 unclassified 列表）中的 domain 长度比较
   （允许 ±2）。不一致记为 `domain_length_mismatch`，仅作提示（CATH 的长度统计口径与"已建模"不完全一致）。
6. **几何合理性校验**：链上每跨一个残基最多 ~3.8 Å。对 junction 两侧残基 `e`、`s`，以及 linker 加锚定残基范围内
   **每一对相邻的已建模残基** `i < j`（中间可以隔着未建模残基），要求 `d(CA_i, CA_j) ≤ 3.8 Å × (j − i) + 0.5 Å`。
   违反说明模型与 SEQRES 不自洽（编号错误，或者像切割后的 serpin 那样环已被切断），该 linker 打 `geometry_violation`，最终集排除。
   0.5 Å 容差来自小样本：真实结构中相邻残基的 CA–CA 在 3.75–3.95 Å 之间（§11.3）。扩展到"每对相邻已建模残基"是完整运行后的修订（§11.5 第 2 条）。
   C–N 轻微拉长、但 CA–CA 仍在范围内的情况只记 `chain_break`。
7. **超长 gap**（> 30 个残基）：可能是 CATH 未分类的结构域或大段无序区，打 `long_gap` 标志，最终集排除。
8. **同一 entry 中同一实体的多条链**（A、B、C…）：每个实体只取一条链（§8 的代表选择会自动保证），避免重复。

---

## 8. 数据集构建

> **v3 note**:§8–§11 记录的是第一阶段的"非冗余集"流程(代表链选择 + linker 级去冗余),作为 provenance 保留。
> v3 起产物改为数据库(§13),**linker 级去冗余的代码(`dedup_linkers` 及其 k-mer 预筛、序列比对)已删除**,
> `results/full/` 不再保留——非冗余改由实例上的 `cluster40` 标记表示。`results/sample50/` 仍用于验证提取规则。

### 8.1 小样本（50 个结构）

从 CATH 多结构域链对应的 PDB entry 中，用固定随机种子抽 **40 个衍射结构 + 10 个 NMR 结构**（共 50 个），
不加分辨率和冗余限制。跑通全流程后打印统计摘要，再人工逐条检查若干 linker。
检查结论和据此修改的规则记录在 §11。

### 8.2 完整数据集

1. **候选链**：CATH 多结构域链 ∩ `pdb_entry_type` 为 `prot`/`prot-nuc` ∩ 分辨率 ≤ 2.5 Å
   （X-ray、EM、中子等**任何报告了分辨率的方法**；NMR 没有分辨率，自然排除）。
2. 用 RCSB GraphQL 把链映射到实体，再用 `clusters-by-entity-40.txt` 映射到 40% 序列聚类。不在任何聚类中的实体单独成簇。
3. **先筛选、后聚类**：只在候选链之间选代表。每簇排序规则：
   GraphQL 中恰好 1 个 UniProt 的优先 → 分辨率高的优先 → domain 数多的优先 → PDB ID / 链号字典序（保证确定性）。
4. 处理排名第一的候选。如果它在链级 QC 中失败（下载失败、非天然、编号不一致、映射失败），
   换同簇下一个候选，每簇最多尝试 3 个。
5. **linker 级去冗余**（弥补实体级聚类的覆盖度漏洞，见 §10 问题 3）。linker 按（分辨率、PDB ID、链、junction 序号）排序，
   依次与已保留的 linker 比较；**同一条链上的不同 junction 永远不互相比较**（串联重复蛋白如锌指的多个 linker 都是独立样本）：
   * 同一个 UniProt accession，且 UniProt 区间（含锚定残基）重叠 → 冗余（`same_uniprot_region`）；
   * 两侧 CATH 超家族对相同（未分类记为 `?`），**整条链**序列做两端空位免罚的全局比对（BLOSUM62，开/延 −10/−0.5），
     一致性 ≥ 40%（分母取较短的链，**不设覆盖度要求**），并且比对能把一个 junction 映射到另一个 junction（±3 残基）
     → 冗余（`homologous_chain`）。
     最初用 linker ± 20 残基的短窗口比较，在冒烟测试中出现了误判（见 §11.4），所以改为整链比对。
     整链比对之前先做一个廉价的预筛：两条链共享的 3-mer 中，落在同一组 21 条相邻比对对角线内的最大数目
     （`diagonal_hits`）不少于 5 才做比对。阈值的标定见 §11.5。
6. 输出：
   * `linkers.parquet`：最终非冗余、通过 QC 的 linker；
   * `linkers_all.parquet`：所有处理过的链上提取出的全部 linker，含 QC 标志和排除原因；
   * `chains.parquet`：每条处理过的链的状态、失败原因、映射校验结果。

---

## 9. 输出 schema（`linkers*.parquet`）

| 列 | 含义 |
|---|---|
| `linker_id` | `{pdb}_{chain}_{junction序号}` |
| `pdb_id`, `auth_asym_id`, `label_asym_id`, `entity_id` | 结构标识 |
| `method`, `resolution`, `n_models` | 实验信息 |
| `domain_source` | 例如 `CATH v4.4.0` |
| `n_domains_chain`, `junction_index`, `junction_type` | junction 信息 |
| `domain_n`, `domain_c`, `cath_n`, `cath_c` | 两侧 domain ID 和 CATH 分类码 |
| `gap_start`, `gap_end`, `gap_length`, `gap_sequence` | gap linker（SEQRES 索引，闭区间；长度为 0 时 start = end + 1） |
| `ss_start`, `ss_end`, `ss_length`, `ss_sequence`, `ss_dssp` | 结构 linker |
| `ss_start_auth`, `ss_end_auth` | 作者编号（含插入码） |
| `ext_n`, `ext_c`, `junction_ss` | 向两侧的延伸长度；零长度 junction 的类别 |
| `flank_n_seq`, `flank_c_seq`, `flank_n_dssp`, `flank_c_dssp` | 结构 linker 两侧各 10 个残基及其 DSSP |
| `window_seq` | linker 加两侧各 20 个残基的序列（上下文） |
| `n_missing`, `frac_observed`, `n_breaks` | 缺失和断裂 |
| `n_nonstandard`, `nonstandard`, `residues_3` | 非标准残基（列表）和 linker 的原始三字母码 |
| `ca_junction`, `ca_span`, `end_to_end` | CA 距离：junction 两侧残基之间 / linker 两端锚定残基之间 / linker 首尾残基之间（Å） |
| `uniprot`, `unp_start`, `unp_end`, `unp_lo`, `unp_hi`, `n_mutations` | UniProt 映射（`unp_lo/hi` 含锚定残基） |
| `numbering_check`, `qc_flags`, `qc_pass` | QC；`qc_flags` 是列表，排除型标志见下 |
| `description`, `taxid`, `cluster`, `dedup_dropped` | 蛋白名称、物种、40% 聚类 ID、去冗余时被谁替代（仅 `linkers_all`） |

**排除型标志**（出现任一个则 `qc_pass = False`）：`long_gap`、`unp_unmapped`、`unp_discontinuous`、`engineered_linker`、
`geometry_violation`、`overlapping_segments`、`no_uniprot`、`chimera`、`not_protein`。
**提示型标志**（不影响 `qc_pass`）：`capped_n/c`、`missing_residues`、`fully_disordered`、`chain_break`、
`nonstandard_residue`、`domain_length_mismatch`、`cath_unclassified`、`minor_segment_dropped`、`numbering_resolved`、`chain0_remap`、`segment_consumed`。

---

## 10. 自审：v0 草稿的问题和修正

逐条检查 v0，下面每一条都已在 v1 正文中修正。

| # | v0 的问题 | 为什么会出错（证据） | v1 的修正 |
|---|---|---|---|
| 1 | linker 只按"两个 domain 之间的残基"定义 | 用 CATH v4.4.0 实际统计：**86.0%** 的 junction gap = 0。v0 得到的数据集绝大多数 linker 是空串，基本没用 | 增加基于 DSSP 的**结构 linker**作为主定义，严格 gap 定义同时保留（§4.2） |
| 2 | 去冗余顺序是"先在全 PDB 上按 40% 聚类选代表（类似 PISCES），再保留多结构域链" | 簇代表通常是分辨率最高的成员，而高分辨率结构经常是**单个结构域的截短构建体**。这样整簇的多结构域成员都会被丢掉，系统性地损失全长多结构域结构 | **先筛选（多结构域、分辨率），后在候选内部选代表**，代表失败时同簇顺延（§8.2 第 3–4 步） |
| 3 | 只做实体级 40% 聚类 | RCSB 实体聚类带覆盖度要求。同一蛋白的"2 结构域构建体"和"3 结构域构建体"覆盖度不够，会落在不同簇，同一个 linker 会在最终集中出现多次 | 增加 **linker 级去冗余**：按 UniProt 区间重叠，以及同超家族对的窗口序列一致性（§8.2 第 5 步） |
| 4 | 没有处理"天然"的要求 | 实例：`2rh1A`（β2AR–T4L 融合，2.4 Å）在 CATH 中是标准的 2-domain 链，它的两个 junction 都是人工的。v0 的过滤条件会把它收进来 | 增加基于 `_struct_ref` / `_struct_ref_seq` / `_struct_ref_seq_dif` 的天然性过滤（§6） |
| 5 | 作者编号找不到就丢弃，其余情况默认映射正确 | 更危险的是**静默的错误映射**：CATH 拆分时依据的旧版文件后来被重新编号；链号 `0`；插入码。找不到只是最容易发现的那种错误 | 用 CDF 和 seqreschopping 双源互校，加上 domain 长度校验和 CA 距离几何校验（§7 第 4–6 条） |
| 6 | 同一 entry 同一实体的多条链没有规定取哪条 | 同源寡聚体（A/B/C/D 链）会产生完全相同的 linker | 每个实体只取一条链；去冗余步骤保证每簇一条（§7 第 8 条） |
| 7 | 非标准残基只用 gemmi 表 | gemmi 表覆盖不全；而且 MSE（实验标记）和 SEP（真实修饰）语义不同，不应混为一谈 | 分层映射，并保留原始三字母码（§5） |
| 8 | 超长 gap 当普通 linker | gap ≥ 30 的 junction 占 1.3%，其中不少是 CATH 未分类的结构域或大段无序区，不是 linker | `long_gap` 标志，最终集排除（§7 第 7 条） |
| 9 | "分辨率 ≤ 2.5 Å"没说是否包含 2.5、是否只限 X-ray | 有歧义，会影响结果是否可复现 | 明确为 **≤ 2.5 Å（含）**，接受任何报告了分辨率的方法（§8.2）。可通过 `--max-resolution` 修改 |

---

## 11. 实际采用的方案和运行记录

### 11.1 实际使用的数据与版本

| 项 | 实际使用 |
|---|---|
| 结构域划分 | **CATH v4.4.0**（latest-release，文件日期 2024-12-16）。CATH 可正常下载和解析，所以**没有触发回退**；SCOPe 2.08 和 Pfam（PDBe API）回退已实现并有测试，但不参与本次输出。没有开启按链回退（`ChainedProvider`），避免把不同来源的 domain 定义混在一个数据集里 |
| CATH 文件 | CDF（作者编号）+ seqreschopping（SEQRES 编号）互校；`cath-domain-list.txt` 提供 C.A.T.H 和长度；`cath-unclassified-list.txt` 提供已切分、未分类 domain 的长度（这些 domain 的 `cath_n/cath_c` 为空，并打 `cath_unclassified` 标志） |
| 结构 | RCSB mmCIF（`.cif.gz`），2026-11-02 下载；没有条目需要回退到 PDB 格式 |
| UniProt 映射 | SIFTS `pdb_chain_uniprot.tsv.gz`（2026-10-25 版），SIFTS 未收录的链回退到 mmCIF `_struct_ref_seq` |
| 分辨率 / 方法 | wwPDB `resolu.idx`（2026-10-23）、`pdb_entry_type.txt` |
| 去冗余 | RCSB `clusters-by-entity-40.txt`（2026-11-02 下载）+ §8.2 第 5 步的 linker 级去冗余 |
| 二级结构 | 内置 numpy DSSP（`linkerdb.dssp`）。没有选 pydssp（它依赖 torch）。mkdssp 不随 Python 包分发，只用于验证 |

### 11.2 DSSP 实现的验证

在 9 个结构（13PK、1CNS、1CRN、1CUK、1TIM、2PTC、2RH1、3PGK、4HHB，共 3,766 个残基）上与 mkdssp 4.2.2 比较：

* 第一版没有实现 β-bulge 合并，三态一致率 99.5%，差异全部是 E 被判成 coil（例如 2PTC 有 10 个）。
* 加入 bulge 合并后，一开始反而**降到 97.8%**，E 过度分配。原因是 DSSP 用无符号整数计算两段 ladder 之间的间隔，
  负间隔会溢出成很大的数、从而不合并；而我的实现用了有符号整数。改成"负间隔不合并"之后，**三态一致率 100%**
  （五态 H/G/I/E/B 一致率 98.3%，差异只在 G/I/B，不影响"规则二级结构"判定）。
* 测试 `test_dssp_matches_mkdssp_reference` 用存储的 mkdssp 结果持续回归（1CNS、2DCH、5NJE）。

### 11.3 小样本（50 个 entry）的检查与规则修订

样本：固定种子抽取 40 个衍射结构 + 10 个 NMR 结构，对应 93 条 CATH 多结构域链（同实体只保留一条后为 57 条）。
第一次运行得到 99 个 junction，86 个通过 QC。逐项检查后发现 6 个问题，修改规则并重跑：

| # | 现象（第一次运行） | 检查结论 | 修改 |
|---|---|---|---|
| 1 | 7 个 `geometry_violation` | 全部是相邻残基（s − e = 1），CA–CA 为 3.85–3.89 Å，链上也没有断裂。是正常的几何波动，阈值 3.85 Å 没有留容差 | 改为 `3.8 Å × (s − e) + 0.5 Å`。重跑后 0 个 |
| 2 | `6r61` 报 `numbering_disagree` | 该条目 2025-03 做过大版本修订（晚于 CATH v4.4.0），SEQRES 前面加了 62 个未建模残基（作者编号从 −62 开始）。CATH 的 SEQRES 切分因此整体偏移，而作者编号映射仍然正确：前者的端点落在未建模残基上，后者的端点全部已建模 | 新增端点裁决规则（§7 第 4 条）；验证了一致链上 312/312 个端点都已建模。`6r61` 变为 `resolved_cdf` |
| 3 | `7kw9`（tRNA 2'-磷酸转移酶，天然蛋白）被判 `no_uniprot` | mmCIF 的 `_struct_ref` 只引用了 `PDB` 自身，但 SIFTS 把它映射到 UniProt A0A7U3ZK32 | UniProt 映射改为 SIFTS 优先（§6）。其余 4 条 `no_uniprot` 是 TCR（6v0y）和 Fab（8smt），SIFTS 中也没有，确认属于"无法验证为天然"，保持排除 |
| 4 | 15 个 linker 的 `cath_n/cath_c` 为空，domain 长度校验也无法进行 | 这些 domain 已切分但尚未分类，只出现在 `cath-unclassified-list.txt` 中。该文件是定宽列，计数大时相邻列会粘连（如 `7ppeA01 … 1110732`） | 读入 unclassified 列表（按列位置从右取长度），打 `cath_unclassified` 标志 |
| 5 | 27% 的非空结构 linker 达到延伸上限（当时 `MAX_EXT = 10`） | 把上限取消后统计单侧延伸：中位数 2，P90 = 12，**P95 = 15**，最大 31。超过 15 的是长 coil 区，例如 PKA 的 C 端尾巴（2f7z） | `MAX_EXT` 改为 15；达到上限的比例降到 11% |
| 6 | `5t6i`（CDPK1 激酶）报出两个互相重叠的 linker，都打了 `segment_consumed` | CATH 把活化环里 8 个残基（214–221，纯 coil，含 DFG 后段）划给了 N 叶，于是在 C 叶内部产生两个伪 junction。PKA 的 C 端尾巴（2f7z，31 个残基纯 coil）也是同一类问题 | 新增"无二级结构的小片段不参与 junction 判定"（§4.1 第 6 条）；junction 类型改为按跨度判断。修改后 `5t6i_A_1` 正好是激酶 hinge `EVYTGGE`，与已知的 hinge 基序一致 |

**人工逐条检查的 linker**（`linkerdb show <id>` 输出序列 / CATH / DSSP / 建模状态四条轨道）：

| linker | 蛋白 | 结论 |
|---|---|---|
| `1x9n_A_4` | DNA 连接酶 I | 螺旋末端和 β 链之间的 `GVG`，合理 |
| `7eb9_A_1` | 线性四聚泛素 | gap 6 = `RLRGG`+`M`，即 Ub C 端尾巴接下一个 Ub 的 M1，是天然的线性泛素连接；作者编号从 72 跳到 101 也正确处理 |
| `5e71_A_2` | tRNA 甲基转移酶 | gap 12 + 两侧 coil，共 23 个残基，gap 中含一段 5 残基螺旋（按规则计入），合理 |
| `5t6i_A_1` | CDPK1 | `EVYTGGE`，激酶 hinge，正确（修订 6 之后） |
| `1yzb_A_2`、`2eza_A_1`、`2eza_A_2` | ataxin-3、酶 I（NMR） | 两侧二级结构元件之间的 coil，合理；NMR 只用第 1 个模型 |
| `1tar_A_2`、`5t6i_A_2` | 天冬氨酸转氨酶、CDPK1 | `helix_spanning`：junction 在 30 残基连续螺旋 / CAD 螺旋中，标注正确 |
| `4pro_A_1` | α-裂解蛋白酶 | `strand_spanning`：CATH 在两个 β 桶共享的一条链中间切开，标注正确 |
| `2dch_X_1`、`5nje_A_2` | 归巢内切酶、PLK1 polo-box | linker 内含 3 / 8 个未建模残基，序列取自 SEQRES，DSSP 为 `?`，`missing_residues` 标志正确 |
| `2rh1_A_*` | β2AR–T4L 融合体 | `excluded_chimera`，所有 junction 都带 `unp_discontinuous`，不会进入最终集 |
| `2f7z_E_3` | PKA | 修订 6 之后不再报出（C 端尾巴是 coil 片段，不是 domain 间 linker） |

**修订后的小样本统计**（`results/sample50/`）：57 条链中 53 条 `ok`、4 条 `excluded_no_uniprot`（TCR/Fab）。
97 个 junction 中 92 个通过 QC。gap 长度为 0 的占 84.8%；结构 linker 长度中位数 7（q25 4、q75 14、q95 23、最大 32），
长度为 0 的占 13.0%（其中 3 个 `helix_spanning`、6 个 `strand_spanning`、3 个 `sse_abutting`）。
linker 组成相对整链：Pro ×1.87、Gly ×1.34、Ser ×1.11 富集；Trp ×0.39、Ile ×0.70、Met ×0.80、Val ×0.81、Leu ×0.82 贫化。
这与文献中 linker 富含 Pro/Gly/Ser、贫疏水残基的结论一致，可以作为整体的合理性检查。

### 11.4 冒烟测试（完整流程，前 60 个候选 entry）中发现的去冗余问题

| 现象 | 结论 | 修改 |
|---|---|---|
| Zif268（1a1h）第 2 个锌指 linker 被判为与同链第 1 个"同源"而删除；色氨酸合酶（1a50）两个相邻 linker 的窗口互相重叠，也被删除 | 同一条链上的不同 junction 是不同的样本，不应互相去冗余 | 同链不比较 |
| 凝血酶重链与胰蛋白酶因为 linker ±20 残基窗口局部一致性 ≥ 40% 被合并 | 短窗口的局部一致性波动大，不能代表"链之间 40% 一致性" | 改为整链比对（不设覆盖度要求），并要求两个 junction 在比对中互相对应（§8.2 第 5 步）。凝血酶/胰蛋白酶整链一致性 41%（按较短链计），仍判为冗余，这符合 40% 的口径 |

### 11.5 完整数据集（≤ 2.5 Å，40% 去冗余）

**流程与规模**（`results/full/`，`run_info.json` 记录全部参数）：

| 步骤 | 数量 |
|---|---|
| 候选链（CATH 多结构域链 ∩ 蛋白 entry ∩ 分辨率 ≤ 2.5 Å） | 102,283 条链，52,514 个 entry |
| 40% 聚类（RCSB 实体级） | 8,125 簇 |
| 下载并处理的 entry | 7,980（第 1 轮）+ 45（第 2 轮）+ 20（第 3 轮，同簇顺延）；下载 736 s，全部 mmCIF，没有触发 PDB 格式回退 |
| 链级结果 | `ok` 7,985；`excluded_chimera` 115；`excluded_no_uniprot` 77；`mapping_failed` 13；没有出现 `numbering_disagree` |
| 编号互校 | `agree` 7,942；`resolved_cdf` 43（均为 CATH 快照之后 SEQRES 变化的条目） |
| 代表链上的 junction | 15,363 个，其中 15,121 个通过 QC。排除原因：`long_gap` 161、`unp_unmapped` 70、`engineered_linker` 36、`overlapping_segments` 8、`geometry_violation` 6 |
| linker 级去冗余 | 删除 1,199 个：`homologous_chain` 1,150、`same_uniprot_region` 49 |
| **最终 `linkers.parquet`** | **13,922 个 linker**，来自 7,207 条链 / 7,094 个 entry / 7,173 个 UniProt / 1,826 个物种；X-ray 13,821、EM 101；分辨率中位数 1.9 Å |

**最终数据集的分布**

* CATH gap 为 0 的占 89.3%；gap > 0 的 1,493 个 linker 平均 gap 8.0。
* 结构 linker：非空的 12,315 个，长度均值 8.8、中位数 7（q25 4、q75 12、q90 18、最大 46）；
  长度为 0 的 1,607 个，其中 `helix_spanning` 740、`strand_spanning` 739、`sse_abutting` 128。
  `strand_spanning` 和螺旋型一样多，看起来反常，检查后确认主要来自 3.40.190.10（II 型周质结合蛋白，739 个中 115 个），
  这类蛋白的两个叶由跨越两侧的 β 链相连，属于已知的铰链结构。
* junction 类型：`sequential` 6,194、`insertion` 5,838、`interleaved` 1,890。
* 4.8% 的非空 linker 含未建模残基；3.0% 含非标准残基（MSE 336、MLY 14、LLP 4、CSO 2、CME 2，主要是晶体学修饰）。
* 组成（相对代表链全序列）：Pro ×2.05、Gly ×1.21、Asn ×1.20、Asp ×1.20、Ser ×1.15、Thr ×1.08 富集；
  Ile ×0.70、Val ×0.72、Leu ×0.79、Met ×0.80、Ala ×0.83 贫化。与小样本和文献中 linker 的组成特征一致。
* 两侧都已分类的 linker 有 10,512 个，覆盖 3,843 种（N 端超家族，C 端超家族）组合；另有 3,410 个至少一侧是 CATH 未分类 domain。

**对输出的检查**

1. `scripts/check_dataset.py`：14 项不变量全部通过，包括：最终集没有排除型标志；分辨率 ≤ 2.5 Å；没有 NMR；
   结构 linker 在 CATH gap 之外不含 H/E/I；linker 两端 CA 距离在物理上可达；每个 40% 簇只有一条代表链；
   不同链之间没有重叠的 UniProt linker 区间；被判为嵌合体的链没有进入最终集。
2. **第一次完整运行发现的问题**：上述检查中"两端 CA 距离物理可达"有 1 例失败：`3nda_A_3`，8 个残基的 linker
   首尾相距 59.5 Å。该结构是**切割后的 serpin**：反应中心环被切开，N 端一半插入 β 折叠 A，
   SEQRES 中相邻的 341 和 345 号残基实际相距 65 Å。原来的几何校验只检查 junction 两侧残基，而它们恰好未建模，所以漏掉了。
   修改：几何校验扩展到 linker 及锚定残基内**每一对相邻的已建模残基**（§7 第 6 条）。重跑后又排除了另外 2 个同类 linker，
   最终集从 13,925 变为 13,922。
3. **DSSP 的扩展验证**：随机抽取 50 条代表链与 mkdssp 4.2.2 比较（1vq8 是核糖体条目，mkdssp 本身处理失败），
   49 条链、17,932 个残基的三态一致率为 **99.87%**，最差的链为 98.3%。
   `dssp_vs_annotation` 中有 26 条链低于 0.8，最低 0.39（4fa8）。检查后确认这些都是富含 β 的结构，其 mmCIF 中**没有 sheet 注释**，
   而我们的 DSSP 与 mkdssp 一致（4fa8 为 98.9%、4ru4 为 100%）。所以低值反映的是沉积注释不完整，不是 DSSP 的错误。
4. **domain 长度提示**：1,057 条代表链（13.2%）有 `domain_length_mismatch`。抽查 200 条：已建模残基数比 CATH 的长度多出中位数 3 个
   （P95 为 8），相对差异超过 10% 的 domain 只占 0.2%。没有看到系统性的映射错误，更像是统计口径不同（或之后的重新精修补建了残基），
   因此这个标志保持为提示型。
5. **去冗余的审计**：被 `homologous_chain` 删除的 1,150 个 linker，整链一致性中位数为 0.44（范围 0.40–1.00），
   两条链的长度比中位数为 0.95，只有 9.4% 低于 0.8。也就是说，第 2 步主要补的不是覆盖度漏洞，而是 RCSB 贪心聚类的非传递性：
   两条链各自属于不同的簇代表，但彼此之间一致性 ≥ 40%。最终集因此满足更严格的"两两之间 < 40%（在同一 junction 上）"。
6. **人工抽查**：除 §11.3 的样本外，又随机检查了 6 个最终 linker，全部合理。其中两个与已知基序吻合：
   原钙黏蛋白 EC4→EC5 的 `DENDNAP` 正是钙黏蛋白重复之间的钙结合基序 DXNDNXP；
   GST 的 `EGTPWLPPPG` 是富含 Pro 的结构域间 linker。另外 MscS（EM，2.5 Å）TM3→β 结构域为 `LFRPFKVGD`；
   II 型 SBP 的 `strand_spanning` 铰链标注也正确。
7. **去冗余预筛的标定**：整链比对约 3.5 ms/对，未经筛选的上限是 265 万对（主要来自 2,185 个两侧都未分类的 linker 组成的大组）。
   在 25,000 对同超家族链上比较：普通的 3-mer 重叠率区分度很差（< 40% 的链对中位数为 0.09，≥ 40% 的最低为 0.099）；
   改用"同一组 21 条相邻对角线内的共享 3-mer 数"，在 351 对 ≥ 40% 的链中最小值为 5。取阈值 5，不漏掉任何一对。
   完整运行中做了 377,767 次比对，跳过 350,175 对，去冗余耗时约 30 分钟。比对结果写入 `data/work/full/alignments.jsonl`，
   中断后可以续算。
8. **任务级续跑**：完整运行期间容器被重启过两次。重跑时 8,190 个链级 checkpoint 全部直接读回（"0 to process"），
   7,980 个结构从缓存读取只用了 4.5 s。

## 12. 已知局限

**覆盖面**

1. **CATH 的覆盖滞后**：CATH v4.4.0 是 2024-12 的快照。分辨率 ≤ 2.5 Å 的 165,297 个蛋白 entry 中，只有 78.5% 至少有一条链进入 CATH；
   其余主要是较新的条目，它们不在数据集里。CATH 的链名只有 5 个字符，**多字符链号**（只有 mmCIF 的大型条目，如核糖体、病毒衣壳）也无法表示，被整体排除。
2. **结构域定义只来自 CATH**。CATH 与 SCOPe 对同一条链的切分可能不同（例如 1cns 在 CATH 中是 2 个 domain，在 SCOPe 中是 1 个）。
   数据集反映的是 CATH 的 domain 观念；SCOPe / Pfam 回退只在 CATH 不可用时启用，本次没有使用。
3. **"天然"的判定依赖 UniProt 映射**：没有 UniProt 条目的天然序列会被排除。最主要的是抗体和 TCR（体细胞重排序列），
   其次是部分宏基因组蛋白。它们的 linker（如 V–C 之间的 elbow）不在本数据集中。
4. **只收晶体 / 电镜结构**（按分辨率筛选）。柔性 linker 在晶体中常常无序（未建模），或者因为太柔而不利于结晶，
   因此数据集**偏向有序、较短、较刚性的 linker**。NMR 结构只在小样本中用于验证流程。

**linker 定义**

5. **严格 gap 定义大多为空**（CATH 把 85% 左右的 junction 切成紧邻），主定义因此依赖 DSSP 和 `MAX_EXT`。
   `MAX_EXT = 15` 是用 50 个 entry 的小样本标定的，超长 coil 区会被截断（带 `capped_*` 标志）；linker 长度的统计分布受这个参数影响。
6. **螺旋型 linker 长度记为 0**（`helix_spanning`）。junction 落在连续螺旋中时，没有客观的边界可以划，所以没有划。
   需要研究螺旋型 linker 的用户应当用 `flank_*_dssp` 和 `junction_ss` 自行处理。长度统计中螺旋型 linker 因此偏短。
7. **DSSP 只在单条链上计算**：链间 β-sheet（domain 交换、二聚体界面）和配体、金属的影响都没有考虑，链界面附近的残基可能被误判为 coil。
   NMR 只用第 1 个模型，没有利用构象系综。
8. **CATH 边界本身有 ±几个残基的不确定性**，未分类（unclassified）domain 的切分没有经过人工审核。
   "无二级结构的 segment 不参与判定"这条规则能去掉一类伪 linker，但不能纠正 CATH 本身的错误切分。
9. **只取第一个 altloc 和微异质性的第一个变体**；修饰残基映射到母体氨基酸（MSE→M、SEP→S）。原始三字母码保留在 `residues_3` 中，
   但没有区分实验引入的修饰（MSE）和天然修饰（磷酸化等）。

**质量与去冗余**

10. **没有使用残基级的结构质量指标**（B 因子、RSCC、Ramachandran 离群等）。linker 处的局部电子密度通常比 domain 核心差，
    已建模的 linker 也可能不可靠。下游如果需要，应结合 wwPDB validation report。
11. **天然性检查有盲区**：同一 UniProt 的人工串联（例如同一蛋白两个 domain 之间换成同源的天然片段）、环化排列体，
    以及 SIFTS 自身映射错误都无法识别。linker 中的点突变只计数（`n_mutations`），不排除。
12. **去冗余口径**：先用 RCSB 实体级 40% 聚类（MMseqs2，带覆盖度要求）每簇取一条链，再做 linker 级去冗余。
    第 2 步用的一致性定义（整链比对、按较短链计、不设覆盖度）与 MMseqs2 不同，实际上比 RCSB 聚类更严格，
    多删除了 7.6% 的 linker（§11.5 第 5 条）；如需严格复现 RCSB 口径，可只用 `linkers_all.parquet` 中代表链的 `qc_pass` linker。
    同源比较只在两侧 CATH 超家族对相同的 linker 之间做，不同超家族标注（例如一侧已分类、另一侧未分类）的同源链不会被比较；
    k-mer 预筛按 25,000 对标定，不能保证零漏检；代表链按分辨率优先选，不一定是最完整的构建体。
13. **CATH 快照与当前 PDB 文件之间的版本差异**：§7 第 4 条的端点裁决只能处理"整体偏移"一类的变化；
    两种编号都自洽但都与当前模型不符的情况（例如重新精修后改了建模范围）无法发现，只能依赖 domain 长度提示和几何校验。

---

## 13. 数据库层（v3，面向融合蛋白 / 生物药设计）

### 13.1 为什么从"非冗余集"改成"数据库"

第一阶段（§8–§11）产出的是**非冗余集**：每个 40% 聚类只留一条代表链。这适合做统计或当训练集,
但对设计不够用——判断一个 linker 是**刚性还是柔性**,必须看它在**多个结构**里的几何变化,而非冗余集恰好把这些副本删掉了。
所以 v3 改成两层:

* **实例层 `linkers_instances.parquet`**:每条链上的每个结构域交界(§1–§7 的规则不变),外加每个实例的
  结构域间几何(`dom_dist` = 两侧结构域 CA 质心距离,`elbow_angle` = 以 linker 处为顶点的夹角)。
* **典型层 `linkers_canonical.parquet`**:把实例按 UniProt 区间归并,每个"典型 linker"一行,汇总刚性、几何、可开发性。
  每行有 `linker_class`:`coil`(有独立 linker 区段)或 `helix_spanning` / `strand_spanning` / `sse_abutting`
  (零长度——两个结构域由一段连续二级结构直接相连,没有独立 linker)。**零长度交界不丢弃**:连续螺旋型连接正是
  刚性 linker 设计里要的(EAAAK 的天然对应),这类条目 `sequence` 存的是跨交界的 flanking SSE 窗口(`coil_sequence` 为空)。
* **非冗余不再靠删除**,而是作为实例上的 `cluster40` 标记保留。要非冗余子集就按标记筛,数据本身不丢。

处理范围也从"每簇一条代表"扩大到**所有多结构域链**(每个 entry 的每个实体取一条,避免把晶体学副本当成独立证据)。

### 13.2 归并:实例 → 典型 linker

按 UniProt accession 分组,组内把 UniProt 区间**重叠**的实例并成一个连通分量(union-find),每个分量是一个典型 linker,
ID 形如 `P12345:185-198`(区间取成员中位数)。不同结构域交界落在不同 UniProt 区段,不会重叠,自然分开。
没有 UniProt 映射的实例无法跨结构归并,单独成为 `grouped = False` 的条目(n_structures = 1)。

### 13.3 刚性判定

关键发现(在 120 个冗余簇的子样本上验证时暴露):**不能用 end-to-end 距离判断刚性**。同一个交界在不同结构里,
CATH/DSSP 画的 linker 边界会差几个残基,导致 end-to-end 在 7 Å 和 25 Å 之间跳(例:AIF `O95831` 的一个交界,
`ss_length` 在 3 和 9 之间变),这是**边界噪声**,不是结构域运动。

改用两个与 linker 边界无关、只由两侧完整结构域决定的量:**结构域质心距离 `dom_dist` 的跨结构变异系数(CV)**
反映"两个结构域是否保持固定距离";**肘角 `elbow_angle` 的跨结构稳健幅度(p90−p10)**反映"结构域是否相对转动"。
后者是必要的,因为质心距离看不出"距离不变、纯粹绕铰链转动"的运动(典型如 MBP 的捕蝇草式开合)——这一点由金标准集暴露(第 14 节)。
肘角的标准差会被占多数的那一种构象稀释,所以用幅度而不是标准差。

* `frac_helix ≥ 0.5` → `rigid_helical`(连续螺旋型 linker);
* 否则,`n_structures ≥ 3` 且有 ≥ 3 个测量:
  `dom_dist_cv ≤ 0.05` 且 `elbow_angle_range ≤ 12°` → `rigid`;
  `dom_dist_cv ≥ 0.15` 或 `elbow_angle_range ≥ 20°` → `flexible`;
  其间 → `intermediate`;置信度 `high`;
* 否则(单 / 少结构):若多数结构 linker 无序 → `flexible_disordered`;再否则按共识二级结构猜(`rigid_helical` /
  `strand_spanning` → 偏刚,`coil` → 偏柔),置信度 `low`。

阈值是启发式的,用第 14 节的种子金标准集标定(金标准集一致率 56% → 70%)。两条局限仍在:
(1)**采样**——若某蛋白在 PDB 里只有一种构象(如我们收到的柠檬酸合酶 13 个结构都是同一态),就测不出它的柔性,
由 `n_structures` / `rigidity_confidence` 提示;(2)肘角幅度能部分补救转动,但真正干净的办法是做结构叠合(roadmap 第 4 条)。

### 13.4 可开发性标记(`developability.py`)

对 linker 序列扫描已知的化学 / 生物液相不稳定位点,作为**提示**而非预测(是否真成问题取决于溶剂暴露和表达体系,
但 linker 恰好是最暴露、最可能反应的地方):N-糖基化序列子 `N-X-[S/T]`、脱酰胺 `N-G`(及较弱的 `N-[S/T/H/A]`)、
天冬氨酸异构化 `D-G`、氧化敏感的 Met/Trp、半胱氨酸(二硫键状态未知)、furin 样切位点 `R-X-[R/K]-R`、赖氨酸数。
每个典型 linker 带这些计数、`liability_tags` 和一个只用于排序的 `liability_score`。

### 13.5 检索(`linkerdb query` / `retrieve.py`)

按设计者的真实问法检索典型 linker:目标结构域间距离(`--distance`,匹配 `dom_dist_median`)、刚性(`--rigidity`)、
是否人源(`--human`)、长度范围、最少结构数(`--min-structures`,用来只取高置信条目)、排除糖基化 / 半胱氨酸、
可开发性上限。输出一个排好序的候选短名单,附刚性、`dom_dist_cv`、证据结构数和液相标记。

**边界**:这个库能缩小候选、给出带实验构象的起点,**不能预测**某个 linker 在你的具体融合蛋白里好不好用——那仍要靠
结构建模加实验。库让第一轮候选更聪明,不替代验证。

### 13.6 输出文件

| 文件 | 内容 |
|---|---|
| `linkers_instances.parquet` | 每条链每个交界一行(§9 的列 + `dom_dist`、`elbow_angle`、`cluster40`、`cath_chain`) |
| `linkers_canonical.parquet` | 每个典型 linker 一行:`canonical_id`、`uniprot`、`unp_lo/hi`、`grouped`、`linker_class`、`n_structures`、`n_entries`、`rep_*`、`best_resolution`、`sequence`、`coil_sequence`、`length_*`、`rigidity`、`rigidity_confidence`、`dom_dist_median`、`dom_dist_cv`、`elbow_angle_median`、`elbow_angle_std`、`elbow_angle_range`、`end_to_end_*`、`superpose_c_rmsd`、`superpose_hinge`、`superpose_n_used`、`superpose_rigidity`、`rigidity_final`、`rigidity_source`、`frac_disordered`、`frac_helix`、`junction_ss/type`、`cath_n/c`、`taxid`、`is_human`、可开发性计数 + `liability_tags` + `liability_score`、`n_clusters40` |
| `chains.parquet` | 每条处理过的链的状态(同 §9) |
| `candidates.parquet` | 候选链及其 40% 簇 |
| `run_info.json` | 版本、规则、计数 |

### 13.7 没做 / 下一步(roadmap)

1. **人工金标准集**:已建**种子集**(~12 条,第 14 节),并据此加了肘角幅度项、把刚性一致率从 56% 提到 70%。
   下一步扩到 100–200 条**专家 / 文献核定**的条目(含核定的残基边界),用来重新标定延伸上限和刚性阈值,并决定螺旋型 linker 的边界。
2. **和人工 linker 对照**:在相同距离需求下,天然 linker 与 (G4S)n、EAAAK 的长度 / 刚性 / 液相标记对比,作为库的价值证明。
3. **扩覆盖面**:加 ECOD 作第二个结构域来源(补 CATH 未覆盖的 22% 和新结构),抗体 / TCR 的 elbow 单列一类,
   AlphaFold 结构作单独一层(带 pLDDT)。
4. **更强的刚性度量**:做结构叠合,捕捉"质心距离不变的纯转动",区分真·刚性和采样不足。
5. **交付**:版本化 parquet + DuckDB/SQLite + 文档;定期增量收录新 PDB 条目;发布前核对各来源许可(PDB CC0,UniProt/CATH 需署名)。

---

## 14. 金标准集(种子)

### 14.1 这是什么、不是什么

`golden/golden_set.csv` 是一个**种子**金标准集(~12 个蛋白),用来检验——并最终标定——数据库的刚性判定。

**诚实声明**:它来自公认的结构生物学常识,**不是带文献引用的人工核定语料,也没有逐残基核定的边界**。
每个条目的 UniProt 都已核对确实在数据库里,但 `expected_rigidity`(rigid / flexible)是定性判断。
它足以暴露系统性错误(已经暴露了"质心距离看不出铰链转动"),但**不足以严格标定阈值**——要用它调参数之前,
必须由专家扩充并核定。CSV 里的 `note` 是行为说明,不是引用。

### 14.2 打分方法(`golden.py` / `linkerdb score-golden`)

对每个金标准蛋白,取它在数据库里结构数最多的那个典型 linker,把我们的刚性标签粗化到 rigid / flexible 轴
(`rigid_helical`、`strand_spanning` → rigid;`flexible_disordered` → flexible;`intermediate` 不计分),
再和期望比较。输出逐条对错 + 总一致率。

### 14.3 结果与由此做的修改

种子集上:**决出胜负的 10 条里 7 条一致(70%)**。

* **对**:钙调蛋白、MBP、磷酸甘油酸激酶(flexible);拓扑异构酶 V、LGP2(flexible);钙粘蛋白、GST(rigid)。
* **仍错(3 条)**:柠檬酸合酶、丙酮酸激酶、肌钙蛋白 C——都被判 rigid,但其实有铰链运动。检查发现这不是度量的错,
  而是**采样**:我们收到的这些蛋白的结构几乎全是同一种构象(柠檬酸合酶 13 个结构肘角幅度仅 1°;肌钙蛋白 C 只有 3 个结构)。
  没有第二种构象,任何度量都测不出柔性。由 `rigidity_confidence` 提示,等 PDB 里有更多构象时自然改善。
* **中间**:GalNAc 转移酶 2、蛋白二硫键异构酶判为 `intermediate`(不计分),与"柔性"方向一致但不确定。

**关键修改**:最初只用 `dom_dist_cv`,MBP(捕蝇草铰链)和 PGK 被判成 rigid——因为它们的结构域是**绕铰链转动**,
质心距离几乎不变(MBP `dom_dist_cv` 0.05),但肘角幅度有 20–22°。于是加了肘角幅度项(§13.3),一致率 56% → 70%,
全局只有约 4% 的 rigid(8,026 → 7,714)被重新分类,没有过拟合式的大摆动。

### 14.4 局限

* 种子集偏"柔性"(大多数多结构域蛋白都有些结构域间活动;刚性的天然 linker 更少)。
* 阈值是在这 12 个点上定的,属于**临时**值。
* `intermediate` 不计分,所以 70% 只是决出胜负那部分的一致率。
* 要真正用金标准集调参(尤其是延伸上限 `max_ext` 和边界),需要带核定残基边界的更大集合——这是 roadmap 第 1 条。

---

## 15. 结构叠合刚性度量(roadmap 第 4 条)

§13.3 的 `dom_dist` + 肘角幅度是**廉价代理**。真正干净、从根上处理转动的度量是**结构叠合**
(`superpose.py`):把一个蛋白的多个结构按 **N 侧结构域**叠合,再看 **C 侧结构域**移动了多少。

### 15.1 方法

对同一个 linker 的两个结构:按 **UniProt 位置**对应残基(不同 PDB 条目的 SEQRES 编号不同)→ Kabsch 叠合 N 侧结构域的
Cα → 把该变换作用到 C 侧结构域 → 量 C 侧 Cα 的 RMSD(`c_rmsd`,同时包含平移和转动)→ 另算两个结构域坐标系之间的相对转角
(`hinge`)。一组的"活动度"取其余结构对最佳分辨率参考结构的最大值。分类:`c_rmsd ≤ 1.5 Å` 且 `hinge ≤ 10°` → rigid;
`c_rmsd ≥ 3 Å` 或 `hinge ≥ 20°` → flexible。每侧至少 5 个共享 Cα 才算,否则 `undetermined`。

### 15.2 金标准集上的验证(`scripts/superpose_golden.py`)

| 蛋白 | 期望 | 代理 | 叠合 | c_rmsd(Å) | hinge(°) |
|---|---|---|---|---|---|
| Calmodulin | flexible | flexible | flexible | 46.7 | 178 |
| Maltose-binding protein | flexible | flexible | flexible | 18.5 | 47 |
| Phosphoglycerate kinase | flexible | flexible | flexible | 13.7 | 35 |
| GalNAc-transferase 2 | flexible | **intermediate** | **flexible** | 20.7 | 70 |
| Protein disulfide-isomerase | flexible | **intermediate** | **flexible** | 13.5 | 35 |
| Topoisomerase V | flexible | flexible | flexible | 3.9 | 17 |
| Glutathione S-transferase P | rigid | rigid | rigid | 0.6 | 2 |
| Citrate synthase | flexible | rigid | rigid | 1.0 | 3 |
| Pyruvate kinase | flexible | rigid | rigid | 1.2 | 3 |
| Troponin C | flexible | rigid_helical | rigid | 0.4 | 0 |
| Cadherin-3 | rigid | rigid | flexible | 8.8 | 21 |
| LGP2 | flexible | flexible | undetermined | — | — |

要点:

* 叠合**原生**判对了 MBP / PGK / 钙调蛋白(靠 hinge 转角,不需要 §13.3 的肘角幅度补丁——它本来就是那个补丁的"正版")。
* 它把代理判不定的两个 `intermediate`(GalNAc-T2、PDI)**解析成了正确的 flexible**,比代理更强。
* 仍判错的 3 个(柠檬酸合酶、丙酮酸激酶、肌钙蛋白 C)叠合后 `c_rmsd` 都 ≈ 1 Å、`hinge` ≈ 0——
  **确认是采样问题**:我们收到的这些蛋白的结构全是同一构象,任何度量都测不出柔性。
* Cadherin-3 叠合判成 flexible,其实不算错:它的结构里含**无钙(apo)**形式,那时 linker 确实是柔的——
  是我金标准里"rigid"这个粗标签把 Ca 结合态当成了全部。说明金标准标签本身也需要分条件核定。
* LGP2 `undetermined`:跨结构共享的核心残基不足 5 个,诚实地不给结论。

### 15.3 已并入数据库

叠合已对全库跑完(`scripts/superpose_db.py`,约 38 分钟,每条链只读一次):对 `n_structures ≥ 3` 的 7,914 条典型 linker
算了 `max_c_rmsd` / `max_hinge`,其中 **7,911 条**算出结果(3 条共享核心残基不足 5 个,`undetermined`)。
`linkers_canonical.parquet` 新增列:`superpose_c_rmsd`、`superpose_hinge`、`superpose_n_used`、`superpose_rigidity`,
以及 **`rigidity_final`**(能算叠合的用叠合,否则回退代理)和 `rigidity_source`(`superposition` / `proxy`)。
`linkerdb query`、金标准打分、`check_database.py` 都改用 `rigidity_final`。

结果:7,911 条用叠合判定、19,355 条回退代理;**3,523 条(约 13%)的最终标签与代理不同**。
`rigidity_final` 分布:flexible 18,027、rigid 5,565、rigid_helical 1,184、intermediate 1,714、flexible_disordered 776
(rigid 从代理的 7,714 降到 5,565——很多"刚性"其实是转动铰链或单构象采样,被叠合纠正)。
金标准集改用 `rigidity_final` 后 **12/12 全部决出胜负、8 条正确**(代理是 10 决出、7 对、2 条不定);
GalNAc-T2、PDI 由 intermediate 正确解析为 flexible,cadherin-3 翻成 flexible(§15.2 的 apo 柔性,属金标准标签过粗)。

`rigidity_confidence`(代理用)对叠合判定不再适用;叠合判定本身即高置信(`n_structures ≥ 3` + 跨构象实测),
由 `rigidity_source = superposition` 标识。单 / 少结构的 linker 仍只有代理判定,置信度低。

---

## 16. 回溯验证:对照 linker-swap 实验(`swaps.py` / `linkerdb score-swaps`)

### 16.1 和金标准集的区别

§14 的金标准集问的是:**我们给天然多结构域蛋白打的刚性标签,对不对得上它公认的构象行为**。
本节问的是另一个问题,也是做融合蛋白的人真正在意的:**库赖以给建议的两个机制模型,对不对得上
"同一融合蛋白换了 linker、测了活性 / 表达 / 寡聚态"的真实结果**。

库只有两类可被 linker-swap 实验证伪的声明:

* **跨距 / 几何**——linker 能把两个结构域撑多远。螺旋沿轴每残基伸展 `HELIX_RISE`(≈1.5 Å);
  柔性 coil 没有固定长度,完全伸展的**上界**是每残基 `EXTENDED_COIL`(≈3.5 Å)。两个常数都在 `artificial.py`。
  由此库可以预测,比如,要跨过某个间距所需的最短柔性 linker 长度。
* **刚性**——linker 是把两个结构域固定在某个相对几何(螺旋 / EAAAK 型),还是让它们采样一大片构象(coil / GS 型)。

很多 linker-swap 的成败取决于库完全不涉及的因素——表达量、蛋白酶敏感性、被融合伙伴本身的生物学。
这类条目记为 `no_geometry_prediction`,**不计入符合率**,和金标准集对 `intermediate` 不给分是同一个道理。
假装库能裁决它们是不诚实的;把库"沉默"的地方标出来,本身就是验证的一部分。

### 16.2 非自证的设计

`golden/linker_swap_set.csv` **只存实验结果**(被测体系、换了哪些 linker、观测量、方向性结论)。
库对每一行的预测由 `swaps.py` 用**库自己的常数和分类约定**现算,再和实验比。
CSV 里绝不出现库的答案,所以这不是循环论证——改动 `HELIX_RISE` / `EXTENDED_COIL` 会直接改变符合率,
这正是我们要的耦合。

每行一个 `expected_kind`,`swaps.py` 据此分支打分:

| `expected_kind` | 库的预测(由模型算出) | 判对条件 |
|---|---|---|
| `min_len_for_span` | 跨过 `span_target_A` 所需的最短 coil 残基数 = ⌈target / 3.5⌉ | 与实验阈值之差 ≤ `tol` |
| `span_insufficient` | `la_len` 的 coil 上界(3.5×len)是否 < `span_target_A` | 预测"够不着"且实验确实寡聚 |
| `span_sufficient` | 同上,是否 ≥ target | 预测"够得着"且实验确实单体 |
| `length_distance_coupling` | 螺旋→`controls`(刚性杆,长度定距离),coil→`weak`(柔性,长度不定距离) | 与实验观察到的耦合方向一致 |
| `more_rigid` | 按类型刚性秩(helix/strand > coil)判哪一侧更刚 | 与实验指认的更刚一侧一致;同秩=平局不计分 |
| `no_geometry_prediction` | 无(超出几何范围) | 不计分(coverage) |

### 16.3 诚实声明(与 §14.1 一致)

随附语料是一个**种子集**。`reference` 列现在是**真实引用**,已对照摘要 / 开放获取全文核实,核实状态记在
`verified` 列(`yes` / `partial` / `no`);但测得的距离和阈值长度仍按**教科书量级**录入(如 ~35 Å 的 VH–VL 间距、
~12 残基的 scFv 单体阈值),**没有**从每篇论文的图里逐点数字化,也没有逐残基核定的边界。它足以显示模型在
充分表征的案例上指没指对方向、以及打分流程是否成立;它**不是**图级精度的基准集。要下定量结论之前,
需扩充更多构建体并把测量值数字化。

种子集上 `score-swaps` 得到 **9 条中 8 条有核实引用;6 条决出、6 条全部一致;3 条超出几何范围(不计分)**。
这个 6/6 是**刻意挑了模型本就该搞定的、充分表征的案例**的结果,它说明的是"模型方向对 + 打分流程可用",
**不是**准确率的证据——就像金标准集的"决出 8/12"一样,样本太小,别过度解读。

* **跨距轴(最硬的检验,给出可证伪的数字)**
  * scFv 单体阈值(Holliger 1993 PNAS 90:6444;Kortt 1997 Protein Eng 10:423,**已核实**):
    VH C 端↔VL N 端间距 ≈ 32–35 Å(scFv 4H11 结构里直接量到 32 Å)。Holliger 原文就同测了 5 和 15 残基 linker
    (15 残基→单体、5 残基→diabody),Kortt 补上 10 残基仍为二聚体、15 残基为单体——所以单体阈值落在 10–15 之间,
    常用口径 ~12。库的 coil 模型给下界 ⌈35/3.5⌉ = 10 残基,与阈值差 2 以内,**且偏小方向正确**:
    3.5 Å/res 是完全伸展上限,真实 coil 没那么伸展,所以实际需要的残基数比下界略多(这正解释了 10 残基为何仍二聚)。
  * (G4S)₁(5 残基,上界 17.5 Å < 35 Å)→ 预测够不着 → 实验为 diabody;
    (G4S)₃(15 残基,上界 52.5 Å ≥ 35 Å)→ 预测够得着 → 实验为单体 scFv(均 Holliger 1993,**已核实**)。
* **刚性轴**
  * EBFP–EGFP FRET 传感器(Arai 2001 Protein Eng 14:529,**已核实**;构建体 A(EAAAK)ₙA,n=2–5):
    FRET 随螺旋加长**单调下降**(即间距增大)、CD 螺旋含量随长度上升;(GGGGS)ₙ 则对长度不敏感(柔性塌缩)。
    库把螺旋判为 `controls`、coil 判为 `weak`,方向一致。**注意**:逐残基 ~1.5 Å 是通用的 α-螺旋轴向上升
    (5.4 Å/圈 ÷ 3.6 残基),属教科书常数、与本文无关,所以数值本身不算独立验证;被验证的是"长度→距离"的方向性耦合。
  * EAAAK 比 GS 更刚:检验库的前提(helix=rigid、coil=flexible)在 Arai 的真实构建体里站得住。
* **超出范围(列出但不计分)**
  * G-CSF–转铁蛋白融合(Bai & Shen 2006 Pharm Res 23:2116,**已核实**):长螺旋 linker H4-2(= A(EAAAK)₄A ×2)
    在 NFS-60 增殖试验里活性最高,EC50 比无 linker 融合体(用 Leu-Glu 二肽相连,Bai 2005 PNAS)低 ~10 倍。
    方向与跨距模型一致,但读出是生物活性,取决于几何以外的太多因素;
  * 转铁蛋白融合插入 ~50 aa 螺旋 linker 后**表达量升 11.2 倍**(Amet 2009 Pharm Res 26:523,**已核实**)——
    库不预测表达,唯一的序列抓手是可开发性(EAAAK 无液相标记);
  * Cys₂His₂ 锌指的 TGEKP 指间 linker 结合 DNA 时变刚(螺旋 capping 贡献结合能),序列 / 长度调节多指 DNA 亲和力
    (点突变降亲和力 10–100 倍)——`partial` 核实,见 §16.5 的更正。要打分需接 DB 查询 + 核定边界(见下)。

### 16.4 局限与下一步

1. 种子集小且偏向"模型搞得定"的充分表征案例;6/6 不能当作准确率。要有统计意义,需扩到几十条带一手引用的受控 swap。
2. 引用已核实,但数值仍是教科书量级。下一步把各论文图里的测量值(逐 n 的 FRET 距离、按构建体的寡聚比例)数字化,
   让 `span_target_A`、`expected` 成为可做亚埃 / 逐残基比对的真值;`verified = partial/no` 的行要补到 `yes`。
3. 目前只对**人工 linker 的物理模型**(1.5 / 3.5 Å/res 和类型→刚性约定)打分,不读数据库。
   `no_geometry_prediction` 里的锌指等案例指向一座该搭的桥:对一个 swap 的几何需求,
   用 `query_canonical` 看库能不能检索出满足该跨距 / 刚性的**天然** linker,把检索能力也纳入回溯验证(roadmap 第 1、2 条)。
4. 跨距模型给的是阈值 / 方向(coil 用上界,螺旋用轴向),不是精确的 end-to-end 分布;`tol` 容差就是为此留的,别拿它做亚埃级比对。

### 16.5 核实记录(2026-11-03,联网核对摘要 / 开放获取全文)

* **已坐实**:Arai 2001 Protein Eng 14(8):529(A(EAAAK)ₙA,n=2–5;FRET 随长度降、CD 螺旋含量随长度升);
  Holliger 1993 PNAS 90(14):6444(同测 5 与 15 残基 linker + 直连);Kortt 1997 Protein Eng 10:423
  (NC10:5/10 残基→二聚、0 残基→三聚、15 残基→单体);Bai & Shen 2006 Pharm Res 23:2116(H4-2 螺旋 linker,EC50 低 ~10 倍);
  Amet 2009 Pharm Res 26:523(螺旋 linker 使 Tf-G-CSF 表达升 11.2 倍)。VH–VL ~32–35 Å 有结构直接测量(4H11:32 Å)佐证。
* **发现并更正的错误**:原种子里锌指那条写的是"**加长** linker 会降低多指 DNA 亲和力"——核实后相反:典型 TGEKP 对理想 B-DNA
  其实**偏短**,结构导向的**更长** linker 反而能显著提升亲和力(有争议);真正稳健的事实是 linker 结合 DNA 时变刚
  (helix capping 贡献结合能)、点突变降亲和力 10–100 倍。已据此改写该行(`verified = partial`,作者归属仍待核。)
* **仍需一手全文**:各论文图里的逐点数值(Arai 每个 n 的距离、各构建体的寡聚比例)未数字化;Bai 2005 PNAS、Huston 1988 PNAS
  按标准引用录入,未逐一开全文。因此数值维持教科书量级,`tol` 容差正为此留。

---

## 17. 刚性 linker 设计系统(`design.py` / `linkerdb design`)

### 17.1 为什么是"检索 + 证据校准 + 外部验证",而不是一个生成网络

有人会问:能不能在这库上训一个从头生成 linker 的模型?我们先做了一个**决策门实验**
(`scripts/discriminator_experiment.py`),因为整个设想的核心假设是"能从序列判断刚性"。结论是**不能**:

* **仅凭序列预测叠合刚性,按 UniProt 分组的 5 折 CV 下 AUC ≈ 0.52**(= 基线 0.71 / 随机),长度单独 0.49;
  连加入结构衍生特征的"泄漏上界"也只有 0.61。随机切分与分组切分同为 0.52——信号弱到没东西可泄漏。
* 更关键的**采样混杂**:rigid 比例随背后结构数**下降**——n=3 时 0.77,n>20 时只剩 **0.47**。
  "rigid"在很大程度上是**构象采样不足**的产物(只抓到一个态就默认判刚),不是 linker 序列的属性。
  这印证了 §13.3 / §14 的发现,现在是全库量化。

所以这个系统**不训练"序列→刚性"分类器**。可信的刚性只有两个来源,系统只用这两个:

1. **你设计进去的 fold**——跨交界的连续 α 螺旋 / β 链本身就是刚性元件(`designed_fold`),这是几何不是预测;
2. **数据库的跨构象证据**——一个 coil 在很多结构(最好是不同功能态)里都没动才算可信的 rigid;低 n 的 rigid 不可信。

### 17.2 三个部件(库当规格 / 来源 / 基准)

* **条件检索 + 受约束多样化(`propose`)**:把设计规格(目标结构域间距、长度、人源、禁糖基化 / 半胱氨酸、
  可开发性上限)变成一个**按证据分层排序**的刚性天然 linker 短名单。复用 `query_canonical` 做规格筛选,
  再按下面的证据层级重排;对 `use = splice` 的 coil linker 可生成**保守单点变体**(同类替换、绝不新增液相标记),
  得到"库里没有、但仍在天然刚性 linker 空间内"的新序列。
* **证据校准(不是学习型判别器)**:每个候选带一个 `trust` 层级——
  `designed_fold`(连续 SSE,构造上刚)> `evidence_strong`(coil,叠合判刚且背后 ≥ 10 个结构)>
  `evidence_weak`(coil,叠合判刚但结构少——**显式标注"可能采样不足"**)> `proxy_only`(仅廉价代理判刚,无跨构象证据)。
  这是**实测证据强度**,天生校准,正是决策门实验告诉我们该用的东西。候选区分 `splice`(可直接嫁接的 coil 序列,
  取自 `coil_sequence`)与 `fold-motif`(零长度 spanning:用一段连续 SSE 直接连接两个结构域,所存序列是侧翼上下文)。
* **外部结构验证(`verify` / `FOLD_PROTOCOL`)**:库不折叠蛋白。既然序列和 DB 先验都预测不了某条新 linker
  在**你这个**融合体里的刚性,真正的裁决交给外部结构预测——而且要出**一组**构象而非单个模型(刚性是跨构象属性)。
  `spans()` 判单个模型里 linker 能否在几何上跨过目标间距;`ensemble_rigidity()` 对预测出的构象组**复用库自身的叠合度量**
  (`superpose.relative_shift` / `classify`,c_rmsd / hinge 阈值完全一致),给出 rigid / flexible 判定。
  `linkerdb design --protocol` 打印完整步骤(AlphaFold-multimer / AF3 / ESMFold 多种子,或 RFdiffusion inpainting 多采样)。

### 17.3 用法

```bash
linkerdb design --distance 25 --human --max-liability 1          # 人源、低液相风险、~25 Å 的刚性连接件
linkerdb design --min-length 5 --max-length 15 --diversify 3     # 可嫁接的短 coil 刚性 linker + 保守变体
linkerdb design --min-trust evidence_strong                       # 只要跨 ≥10 结构验证过的 rigid
linkerdb design --protocol                                        # 打印外部结构验证协议
```

库在这里扮演的角色:**规格与检索底座、刚性 linker 的来源**(4,528 条高置信刚性天然 linker,其中 2,597 条是
构造上刚的连续 SSE 连接)、几何 **oracle**,以及回溯**基准**(`score-golden` / `score-swaps`)。

### 17.4 边界与下一步

1. **这不预测融合体成败**。它缩小并排序候选、给出带实测证据和液相标记的起点;某条 linker 在你具体构建体里好不好用,
   仍要靠外部结构预测 + 实验。排序是**证据层级**,不是成功概率。
2. **"刚性"受采样限制**(§17.1):库说某天然 linker 刚,只代表 PDB 没采到它动;`evidence_weak` 层级就是这个风险的显式标注。
3. **多样化很保守**(只做同类单点替换、不新增液相标记),刻意不离开天然序列空间;要更激进的从头骨架生成,
   应走 RFdiffusion 一类的结构生成,再用本系统的 oracle 和基准来筛与评,而不是指望这 2.7 万条 linker 去训一个更强的生成器。
4. 柔性 linker 没有对应的设计系统:(G4S)ₙ 已够用,且库里柔性样本恰恰缺得最多(§12)。本系统专做**刚性**连接件,
   因为那正是合成工具箱最薄(基本只有 EAAAK)、而库的刚性信号最可信、价值最大的地方。
