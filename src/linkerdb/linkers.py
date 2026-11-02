"""Linker extraction: map domain boundaries onto a chain and cut out the
inter-domain segments (docs/design.md §4–§7).

Conventions: positions are SEQRES indices (1-based, inclusive). An empty
interval is written ``start = end + 1``.
"""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np

from . import dssp
from .config import Rules
from .domains.base import ChainDomains, Segment
from .structure import Chain

ENGINEERED_DIF = ("linker", "insertion", "expression tag", "cloning", "deletion", "initiating", "acetylation",
                  "amidation", "chromophore", "leader sequence")
MUTATION_DIF = ("engineered mutation", "conflict", "variant", "microheterogeneity", "see remark")

# linker-level flags that exclude a linker from the final dataset
EXCLUDING_FLAGS = frozenset({"long_gap", "unp_unmapped", "unp_discontinuous", "engineered_linker",
                             "geometry_violation", "overlapping_segments",
                             # chain-level natural-protein failures, copied onto each linker
                             "no_uniprot", "chimera", "not_protein"})


class ChainSkip(Exception):
    """Chain-level QC failure; ``reason`` is recorded in chains.parquet."""

    def __init__(self, reason: str, detail: str = "") -> None:
        super().__init__(f"{reason}: {detail}" if detail else reason)
        self.reason = reason
        self.detail = detail


@dataclass
class MappedDomain:
    domain_id: str
    classification: str | None
    reported_length: int | None
    segments: list[tuple[int, int]]

    @property
    def span(self) -> tuple[int, int]:
        return min(s for s, _ in self.segments), max(e for _, e in self.segments)


@dataclass
class MappedChain:
    domains: list[MappedDomain]
    fragments: list[tuple[int, int]]
    numbering_check: str
    flags: list[str] = field(default_factory=list)


# ---------------------------------------------------------------- mapping
def _map_auth_segments(segs: list[Segment], index: dict[tuple[int, str], int]) -> list[tuple[int, int]] | None:
    out = []
    for seg in segs:
        if seg.start_auth is None or seg.end_auth is None:
            return None
        a = index.get((seg.start_auth.num, seg.start_auth.ins))
        b = index.get((seg.end_auth.num, seg.end_auth.ins))
        if a is None or b is None or b < a:
            return None
        out.append((a, b))
    return out


def _seqres_segments(segs: list[Segment], n: int) -> list[tuple[int, int]] | None:
    out = []
    for seg in segs:
        if not seg.resolved() or not 1 <= seg.start <= seg.end <= n:
            return None
        out.append((seg.start, seg.end))
    return out


def _endpoints_observed(mapping: dict[str, list[tuple[int, int]]], obs: np.ndarray) -> bool:
    return all(obs[s - 1] and obs[e - 1] for segs in mapping.values() for s, e in segs)


def map_domains(chain: Chain, cd: ChainDomains, rules: Rules = Rules()) -> MappedChain:
    """Resolve domain segments to SEQRES indices and cross-check numbering schemes."""
    index = chain.auth_index()
    auth_map: dict[str, list[tuple[int, int]]] | None = {}
    for d in cd.domains:
        if all(s.resolved() for s in d.segments):  # source already in SEQRES numbering (Pfam)
            segs = _seqres_segments(d.segments, chain.n)
        else:
            segs = _map_auth_segments(d.segments, index)
        if segs is None:
            auth_map = None
            break
        auth_map[d.domain_id] = segs
    seq_map: dict[str, list[tuple[int, int]]] | None = None
    if cd.seqres_domains is not None:
        seq_map = {}
        for d in cd.seqres_domains:
            segs = _seqres_segments(d.segments, chain.n)
            if segs is None:
                seq_map = None
                break
            seq_map[d.domain_id] = segs

    obs = chain.observed
    # CATH chops on coordinates, so every segment endpoint of a correct mapping is a
    # modelled residue (verified: 312/312 endpoints in the 50-entry sample).
    two_files = cd.seqres_domains is not None
    if auth_map is not None and seq_map is not None:
        check = "agree" if auth_map == seq_map else "disagree"
    elif auth_map is not None:
        check = "cdf_only" if two_files else "single_source"
    elif seq_map is not None:
        check = "seqres_only"
    else:
        raise ChainSkip("mapping_failed", "domain boundaries could not be located in the structure")
    if check == "disagree":
        ok_auth, ok_seq = _endpoints_observed(auth_map, obs), _endpoints_observed(seq_map, obs)
        if ok_auth and not ok_seq:
            check = "resolved_cdf"
        elif ok_seq and not ok_auth:
            check = "resolved_seqres"
        else:
            diff = sorted(k for k in set(auth_map) | set(seq_map) if auth_map.get(k) != seq_map.get(k))
            raise ChainSkip("numbering_disagree", ",".join(diff))
    chosen = auth_map if check in ("agree", "cdf_only", "single_source", "resolved_cdf") else seq_map
    if two_files and not _endpoints_observed(chosen, obs):
        raise ChainSkip("mapping_failed", "CATH segment endpoint maps to an unmodelled residue")
    info = {d.domain_id: d for d in cd.domains}
    if cd.seqres_domains:
        for d in cd.seqres_domains:
            info.setdefault(d.domain_id, d)
    domains = [MappedDomain(did, info[did].classification, info[did].length, sorted(segs))
               for did, segs in chosen.items()]
    fragments = _map_auth_segments(cd.fragments, index) or [] if cd.fragments else []

    flags = ["numbering_resolved"] if check.startswith("resolved") else []
    if any(d.classification is None for d in domains):
        flags.append("cath_unclassified")
    for d in domains:
        if d.reported_length is None:
            continue
        n_obs = sum(int(obs[s - 1:e].sum()) for s, e in d.segments)
        if abs(n_obs - d.reported_length) > rules.domain_len_tol:
            flags.append("domain_length_mismatch")
            break
    return MappedChain(domains, fragments, check, flags)


# ---------------------------------------------------------------- junctions
@dataclass
class Junction:
    index: int
    n_dom: MappedDomain
    c_dom: MappedDomain
    n_seg: tuple[int, int]
    c_seg: tuple[int, int]
    kind: str


def junctions(mapped: MappedChain) -> tuple[list[Junction], list[str]]:
    pieces = sorted((s, e, k) for k, d in enumerate(mapped.domains) for s, e in d.segments)
    flags = []
    out: list[Junction] = []
    for (s1, e1, k1), (s2, e2, k2) in zip(pieces, pieces[1:]):
        if s2 <= e1:
            flags.append("overlapping_segments")
        if k1 == k2:
            continue  # gap inside one domain
        a, b = mapped.domains[k1], mapped.domains[k2]
        (sa, ea), (sb, eb) = a.span, b.span
        if ea < sb or eb < sa:
            kind = "sequential"      # spans do not overlap (internal gaps do not matter)
        elif (sa <= sb and eb <= ea) or (sb <= sa and ea <= eb):
            kind = "insertion"       # one domain sits inside the other's span
        else:
            kind = "interleaved"     # spans overlap without nesting
        out.append(Junction(len(out) + 1, a, b, (s1, e1), (s2, e2), kind))
    return out, sorted(set(flags))


def drop_minor_segments(mapped: MappedChain, ss: str) -> tuple[MappedChain, list[tuple[str, int, int]]]:
    """Treat domain segments without any regular secondary structure as unassigned.

    CATH sometimes books a short coil excursion (e.g. part of a kinase activation
    loop) to the other domain; its two junctions would produce overlapping
    pseudo-linkers around a domain-internal loop. Such segments carry no SSE that
    could anchor a linker, so they are removed before junction detection.
    """
    dropped = []
    domains = []
    for d in mapped.domains:
        keep = []
        for s, e in d.segments:
            if any(dssp.is_regular(c) for c in ss[s - 1:e]):
                keep.append((s, e))
            else:
                dropped.append((d.domain_id, s, e))
        if keep:
            domains.append(MappedDomain(d.domain_id, d.classification, d.reported_length, keep))
    flags = list(mapped.flags) + (["minor_segment_dropped"] if dropped else [])
    return MappedChain(domains, mapped.fragments, mapped.numbering_check, flags), dropped


# ---------------------------------------------------------------- linkers
def _seq(chain: Chain, a: int, b: int) -> str:
    return "".join(r.one for r in chain.residues[max(a, 1) - 1:max(b, 0)])


def _ss(ss: str, a: int, b: int) -> str:
    return ss[max(a, 1) - 1:max(b, 0)]


def _ca(chain: Chain, i: int) -> np.ndarray | None:
    if not 1 <= i <= chain.n:
        return None
    p = chain.coords[i - 1, 1]
    return None if np.isnan(p).any() else p


def _dist(p, q) -> float | None:
    if p is None or q is None:
        return None
    return float(np.linalg.norm(p - q))


def _unp_check(chain: Chain, a: int, b: int) -> tuple[str | None, list[int | None], str | None]:
    """UniProt positions of residues a..b; returns (accession, positions, failure-flag)."""
    pos: list[int | None] = []
    accs = set()
    for i in range(max(a, 1), min(b, chain.n) + 1):
        hit = chain.unp_position(i)
        if hit is None:
            pos.append(None)
        else:
            accs.add(hit[0])
            pos.append(hit[1])
    if any(p is None for p in pos):
        return (accs.pop() if len(accs) == 1 else None), pos, "unp_unmapped"
    if len(accs) != 1:
        return None, pos, "unp_discontinuous"
    if any(q - p != 1 for p, q in zip(pos, pos[1:])):
        return accs.pop(), pos, "unp_discontinuous"
    return accs.pop(), pos, None


def extract_linkers(chain: Chain, mapped: MappedChain, ss: str, rules: Rules = Rules(),
                    extra_flags: list[str] | None = None) -> list[dict]:
    mapped, _ = drop_minor_segments(mapped, ss)
    jxs, chain_flags = junctions(mapped)
    chain_flags = chain_flags + list(extra_flags or [])
    brk = chain.breaks(rules.peptide_cn_max, rules.ca_ca_max)
    obs = chain.observed
    records = []
    for jx in jxs:
        e, s = jx.n_seg[1], jx.c_seg[0]
        flags = list(chain_flags) + list(mapped.flags)
        gap_len = s - e - 1
        if gap_len < 0:
            gap_len = 0
        if gap_len > rules.long_gap:
            flags.append("long_gap")

        # structural linker: grow over non-regular SS within the flanking segments
        left, ext_n = e + 1, 0
        while left - 1 >= jx.n_seg[0] and ext_n < rules.max_ext and not dssp.is_regular(ss[left - 2]):
            left -= 1
            ext_n += 1
        capped_n = ext_n == rules.max_ext and left - 1 >= jx.n_seg[0] and not dssp.is_regular(ss[left - 2])
        right, ext_c = s - 1, 0
        while right + 1 <= jx.c_seg[1] and ext_c < rules.max_ext and not dssp.is_regular(ss[right]):
            right += 1
            ext_c += 1
        capped_c = ext_c == rules.max_ext and right + 1 <= jx.c_seg[1] and not dssp.is_regular(ss[right])
        if capped_n:
            flags.append("capped_n")
        if capped_c:
            flags.append("capped_c")
        if left == jx.n_seg[0] or right == jx.c_seg[1]:
            flags.append("segment_consumed")
        ss_len = right - left + 1

        if ss_len == 0:
            se, ss_s = ss[e - 1], ss[s - 1]
            no_break_es = not brk[e - 1] if e - 1 < len(brk) else True
            if se in "HI" and ss_s in "HI" and no_break_es:
                jss = "helix_spanning"
            elif se == "E" and ss_s == "E" and no_break_es:
                jss = "strand_spanning"
            else:
                jss = "sse_abutting"
        else:
            jss = "coil"

        # missing residues / breaks over the linker plus its anchor residues
        lo, hi = max(left - 1, 1), min(right + 1, chain.n)
        lk_obs = obs[left - 1:right] if ss_len > 0 else np.zeros(0, dtype=bool)
        n_missing = int((~lk_obs).sum())
        n_breaks = int(brk[lo - 1:hi - 1].sum()) if hi > lo else 0
        if ss_len > 0 and n_missing == ss_len:
            flags.append("fully_disordered")
        if n_missing:
            flags.append("missing_residues")
        if n_breaks:
            flags.append("chain_break")

        # geometric plausibility across the junction (catches numbering errors)
        # A chain cannot span more than ~3.8 Å per residue, whatever is missing in between,
        # so a violation means the SEQRES/author numbering of the model is inconsistent.
        d_es = _dist(_ca(chain, e), _ca(chain, s))
        if d_es is not None and s > e and d_es > rules.ca_step_max * (s - e) + rules.ca_tol:
            flags.append("geometry_violation")

        residues = chain.residues[left - 1:right]
        nonstd = sorted({r.mon_id for r in residues if r.nonstandard})
        if nonstd:
            flags.append("nonstandard_residue")

        acc, unp_pos, unp_flag = _unp_check(chain, lo, hi)
        if unp_flag:
            flags.append(unp_flag)
        difs = [d for d in chain.seq_difs if lo <= d.seq_id <= hi]
        if any(any(k in d.details for k in ENGINEERED_DIF) for d in difs):
            flags.append("engineered_linker")
        n_mut = sum(1 for d in difs if any(k in d.details for k in MUTATION_DIF))

        def unp_at(i):
            k = i - lo
            return unp_pos[k] if 0 <= k < len(unp_pos) else None

        fl = rules.flank
        rec = {
            "linker_id": f"{chain.entry.pdb_id}_{chain.auth_chain}_{jx.index}",
            "pdb_id": chain.entry.pdb_id,
            "auth_asym_id": chain.auth_chain,
            "label_asym_id": chain.label_asym,
            "entity_id": chain.entity_id,
            "method": chain.entry.method,
            "resolution": chain.entry.resolution,
            "n_models": chain.entry.n_models,
            "chain_length": chain.n,
            "n_domains_chain": len(mapped.domains),
            "junction_index": jx.index,
            "junction_type": jx.kind,
            "domain_n": jx.n_dom.domain_id,
            "domain_c": jx.c_dom.domain_id,
            "cath_n": jx.n_dom.classification,
            "cath_c": jx.c_dom.classification,
            "gap_start": e + 1,
            "gap_end": s - 1 if s - 1 >= e + 1 else e,
            "gap_length": gap_len,
            "gap_sequence": _seq(chain, e + 1, s - 1),
            "ss_start": left,
            "ss_end": right,
            "ss_length": ss_len,
            "ss_sequence": _seq(chain, left, right),
            "ss_dssp": _ss(ss, left, right),
            "ss_start_auth": chain.res(left).auth if ss_len else "",
            "ss_end_auth": chain.res(right).auth if ss_len else "",
            "ext_n": ext_n,
            "ext_c": ext_c,
            "junction_ss": jss,
            "flank_n_seq": _seq(chain, left - fl, left - 1),
            "flank_c_seq": _seq(chain, right + 1, right + fl),
            "flank_n_dssp": _ss(ss, left - fl, left - 1),
            "flank_c_dssp": _ss(ss, right + 1, right + fl),
            "n_missing": n_missing,
            "frac_observed": float(lk_obs.mean()) if ss_len else 1.0,
            "n_breaks": n_breaks,
            "n_nonstandard": sum(1 for r in residues if r.nonstandard),
            "nonstandard": nonstd,
            "residues_3": [r.mon_id for r in residues],
            "ca_junction": d_es,
            "ca_span": _dist(_ca(chain, left - 1), _ca(chain, right + 1)),
            "end_to_end": _dist(_ca(chain, left), _ca(chain, right)) if ss_len else None,
            "uniprot": acc,
            "unp_start": unp_at(left) if ss_len else unp_at(s),
            "unp_end": unp_at(right) if ss_len else unp_at(e),
            "unp_lo": unp_pos[0] if unp_pos else None,
            "unp_hi": unp_pos[-1] if unp_pos else None,
            "n_mutations": n_mut,
            "window_seq": _seq(chain, lo - rules.dedup_window, hi + rules.dedup_window),
            "numbering_check": mapped.numbering_check,
            "qc_flags": sorted(set(flags)),
        }
        rec["qc_pass"] = not (set(rec["qc_flags"]) & EXCLUDING_FLAGS)
        records.append(rec)
    return records


def chain_ss(chain: Chain) -> str:
    return dssp.assign(chain.coords, np.array([r.mon_id == "PRO" for r in chain.residues]))


def natural_problem(chain: Chain) -> tuple[str, str] | None:
    """Chain-level natural-protein check (docs/design.md §6): ``(reason, detail)`` or None."""
    if chain.polymer_type and chain.polymer_type != "polypeptide(L)":
        return "not_protein", chain.polymer_type
    accs = chain.uniprot_accessions()
    if not accs:
        return "no_uniprot", ",".join(sorted({r.db_name for r in chain.refs})) or "no reference"
    if len(accs) > 1:
        return "chimera", ",".join(accs)
    return None
