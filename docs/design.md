# linkerdb 数据获取层设计

> 范围：本文件只覆盖**数据获取层**：从 PDB 下载结构、解析结构域边界、提取 domain linker、
> 构建非冗余数据集并输出 parquet。下游的特征计算和建模不在本会话范围内。
>
> 版本记录：v0 是初稿（见 git 历史）；v1 是自审后的修订版，第 10 节列出自审发现的问题和对应修正；
> v2 根据 50 结构小样本的人工检查修改了若干规则，正文已更新为最终规则，修改经过记录在第 11 节；
> 第 12 节是已知局限。

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
3. **不连续 domain / 插入型 domain**：`junction_type = insertion`，并记录 junction 序号和 domain 对。
   同一对 domain 可能有两个 junction，下游统计时要注意。
4. **编号不一致**：CDF（作者编号）映射到 SEQRES 后与 seqreschopping 逐个比对 segment 端点：
   * 全部一致 → `numbering_check = agree`；
   * 不一致时，利用"CATH 只按坐标切分，所以正确映射的每个 segment 端点都必须是已建模残基"这一事实裁决
     （小样本中一致的 57 条链 312/312 个端点都已建模）：只有一方的端点全部已建模 → 采用该方
     （`resolved_cdf` / `resolved_seqres`，打 `numbering_resolved` 标志）；两方都满足或都不满足 → `numbering_disagree`，整链排除；
   * 只有其中一个可用 → `cdf_only` / `seqres_only`；
   * 对 CATH，最终采用的映射如果有端点落在未建模残基上，整链按 `mapping_failed` 排除。
5. **domain 长度校验**：映射后每个 domain 内已建模残基数与 CATH 列表（含 unclassified 列表）中的 domain 长度比较
   （允许 ±2）。不一致记为 `domain_length_mismatch`，仅作提示（CATH 的长度统计口径与"已建模"不完全一致）。
6. **几何合理性校验**：junction 两侧残基 `e`、`s` 都已建模时，要求 `d(CA_e, CA_s) ≤ 3.8 Å × (s − e) + 0.5 Å`。
   链上每跨一个残基最多 ~3.8 Å，所以违反说明模型的编号前后不一致；该 linker 打 `geometry_violation`，最终集排除。
   0.5 Å 容差来自小样本：真实结构中相邻残基的 CA–CA 在 3.75–3.95 Å 之间（§11）。
7. **超长 gap**（> 30 个残基）：可能是 CATH 未分类的结构域或大段无序区，打 `long_gap` 标志，最终集排除。
8. **同一 entry 中同一实体的多条链**（A、B、C…）：每个实体只取一条链（§8 的代表选择会自动保证），避免重复。

---

## 8. 数据集构建

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
    linker 级的同源比较只在两侧 CATH 超家族对相同的 linker 之间做，不同超家族标注（含一侧未分类）的同源链不会被比较；
    代表链按分辨率优先选，不一定是最完整的构建体。
13. **CATH 快照与当前 PDB 文件之间的版本差异**：§7 第 4 条的端点裁决只能处理"整体偏移"一类的变化；
    两种编号都自洽但都与当前模型不符的情况（例如重新精修后改了建模范围）无法发现，只能依赖 domain 长度提示和几何校验。
