# linkerdb 数据获取层设计

> 范围：本文件只覆盖**数据获取层**：从 PDB 下载结构、解析结构域边界、提取 domain linker、
> 构建非冗余数据集并输出 parquet。下游的特征计算和建模不在本会话范围内。
>
> 版本记录：v0 是初稿（见 git 历史）；v1 是自审后的修订版，第 10 节列出自审发现的问题和对应修正；
> 第 11 节记录实际运行时采用的方案和参数；第 12 节是已知局限。

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
| 天然性判定 | mmCIF 自身 | `_struct_ref`、`_struct_ref_seq`、`_struct_ref_seq_dif` |

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
5. Junction 分类：
   * `sequential`：两个 domain 都是连续的，并且在序列上不交错；
   * `insertion`：其中一个 domain 是不连续的（插入或回折，如 `1cnsA01 = 1-87,146-243` 中间插着 `1cnsA02 = 88-145`）。
     这种情况下会产生 A→B 和 B→A 两个 junction。

### 4.2 两种 linker 定义（同时输出）

**数据事实**（CATH v4.4.0，全部 358,437 个 junction）：**86.0% 的 junction 两侧 domain 直接相接（gap = 0）**，
gap 为 1–10 的占 9.4%，≥ 30 的占 1.3%。CATH 会把绝大多数残基都分配给某个 domain，
所以按字面意思取"两个 domain 之间的残基"，大多数 linker 都是空的。为此同时输出两种定义：

1. **gap linker（严格定义）**：`e+1 … s-1`，即 CATH 未分配给任何 domain 的残基，长度 ≥ 0。
2. **结构 linker（主定义，`ss_*` 列）**：从 junction 出发，向 N 侧和 C 侧各自延伸，经过的残基
   都是**非规则二级结构**（不是 α 螺旋 H，也不是 β 链 E），遇到第一个规则二级结构残基就停。
   直观理解：结构 linker 是 N 侧 domain 最后一个二级结构元件和 C 侧 domain 第一个二级结构元件之间的那段。
   * 二级结构由内置的 DSSP 算法（Kabsch–Sander 氢键能量）根据第 1 个模型的骨架坐标计算。
     未建模残基没有坐标，按"非规则"处理。
   * 延伸只能在两侧 segment 内部进行，不跨进第三个 domain。每侧最多延伸 `MAX_EXT` 个残基
     （默认 10，经 §8 小样本检查后确认），到达上限就打 `capped_n` / `capped_c` 标志。
   * gap 内部的残基无论二级结构如何都算 linker。
   * 如果 junction 两侧紧邻的残基都处在**同一段连续螺旋**里，结构 linker 长度为 0，
     `junction_ss = 'helix_spanning'`（常见的"螺旋型 linker"，如钙调蛋白中央螺旋）。
     这类情况下没有不依赖人为阈值的 linker 边界，因此不强行划定。

### 4.3 每个 linker 记录的内容

位置（SEQRES 和作者编号）、长度、规范一字母序列、原始三字母残基、两侧各 10 个 flanking 残基及其二级结构、
DSSP 二级结构字符串、缺失残基数、链断裂、非标准残基、UniProt 区间、两侧 domain 的 CATH 分类码等。
完整 schema 见 §9。

---

## 5. 特殊情况

| 情况 | 处理 |
|---|---|
| 缺失残基（SEQRES 中有、坐标中没有） | 序列始终取自 SEQRES，记录 `n_missing` / `frac_observed`。未建模残基的 DSSP 状态记为 `?`（按非规则处理）。linker 全部未建模时打 `fully_disordered` 标志，保留（柔性 linker 正是研究对象），由下游决定是否使用 |
| 链断裂 | 在相邻 SEQRES 位置都已建模的前提下，`C(i)–N(i+1) > 2.0 Å` 记为断裂；缺骨架原子时改用 `CA–CA > 4.2 Å`。SEQRES 上有未建模残基导致的缺口记为 `gap_break`。分别记录 linker 区间（含 ±1 flanking）内的 `n_breaks`。DSSP 中断裂处不计算跨断点的 i−1 氢原子位置 |
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
这类 linker 是人工设计的，必须排除。判定全部基于 mmCIF 中的参考序列注释：

* **链级（不通过则整条链排除）**
  * 该实体的 `_struct_ref` 中 UniProt（`db_name = UNP`）accession **恰好 1 个**。
    0 个记为 `no_uniprot`（包括从头设计蛋白和合成构建体）；≥ 2 个记为 `chimera`（融合体或嵌合体）。
* **linker 级（打标志，最终集排除）**
  * linker（含两侧各 1 个残基）中每个残基都能通过 `_struct_ref_seq` 映射到 UniProt，
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
   * 有任何不一致 → `disagree`，整链排除（无法判断哪个对）；
   * 只有其中一个可用 → `cdf_only` / `seqres_only`，保留并打标志。
5. **domain 长度校验**：映射后每个 domain 内已建模残基数要与 `cath-domain-list.txt` 中的 domain 长度一致
   （允许 ±2）。不一致记为 `domain_length_mismatch`。
6. **几何合理性校验**：junction 两侧残基 `e`、`s` 都已建模时，要求 `d(CA_e, CA_s) ≤ 3.85 Å × (s − e)`。
   违反说明编号或映射有错，整链排除。
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
4. 处理排名第一的候选。如果它在链级 QC 中失败（下载失败、非天然、编号不一致、几何校验失败），
   换同簇下一个候选，每簇最多尝试 3 个。
5. **linker 级去冗余**（弥补实体级聚类的覆盖度漏洞，见 §10 问题 3）：
   * 同一个 UniProt accession 且 UniProt 区间重叠的 linker，只保留排名靠前那条链上的；
   * 两侧 CATH 超家族对相同的 linker，取"linker ± 20 个 flanking 残基"窗口做全局比对，
     一致性 ≥ 40%（分母取较短序列）则视为冗余，保留排名靠前的。
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
| `ext_n`, `ext_c`, `capped_n`, `capped_c`, `junction_ss` | 延伸信息 |
| `flank_n_seq`, `flank_c_seq`, `flank_n_dssp`, `flank_c_dssp` | 各 10 个 flanking 残基 |
| `n_missing`, `frac_observed`, `n_breaks`, `fully_disordered` | 缺失和断裂 |
| `n_nonstandard`, `nonstandard`, `residues_3` | 非标准残基 |
| `ca_span`, `end_to_end` | linker 两端 CA 距离等简单几何量 |
| `uniprot`, `unp_start`, `unp_end`, `unp_ok`, `n_mutations` | UniProt 映射 |
| `numbering_check`, `qc_flags`, `qc_pass` | QC |

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

## 11. 实际采用的方案

（运行后填写，见下文）

## 12. 已知局限

（完成全部流程后填写）
