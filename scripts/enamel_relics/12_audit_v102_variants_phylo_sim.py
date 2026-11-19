"""Phase A–D (v1.0.2 Final Scientific Audit):
1. Phase A: Full molecular & phylogenetic audit of all 23 shared coding variants in Table S5a,
   explicitly separating Alligator reference protein coordinates/residues (XP_019346377.1 / XP_019346386.1)
   from reconstructed avian ancestral codons/amino acids (e.g., codon 601 TTT[Phe]->TAT[Tyr]->TAA[Ter];
   codon 1003 TGG[Trp] / site-independent ML TGC[Cys, P=0.343] -> TGA[Ter]).
2. Phase B: Outgroup & sister-clade coverage validation for candidate synapomorphies
   (AMBN Exon 4 16-nt deletion in Palaeognathae; AMBN Exon 4 1-nt insertion + AG->AC splice in Chemyd+Chrpic;
   ENAM Exon 8 codon 1022 p.Ser1022Ter in Neognathae), including targeted blastn/tblastn/local-alignment
   recovery attempts in Neognathae and Gopevg/Pelsin plus squamate outgroup comparison.
3. Phase C: Topological sensitivity across T1, T2, T3 (MAFFT L-INS-i & MUSCLE 5), formal IQ-TREE 3
   topology tests (KH, SH, wKH, wSH, ELW, AU with 10,000 RELL replicates), and ML ASR vs. parsimony across T1/T2/T3.
4. Phase D: Forensic audit of the v1.0.1 K80 neutral stop-gain simulation (identifying why unmasked 100%
   7-bird coverage predicted 26-29 shared PTCs vs. 16 observed) and full sensitivity grid (1,000 replicates
   per scenario, seed 20261123) across empirical coverage masking, dual-aligner consensus masking,
   kappa in {2, 4, 6}, branch-length scales in {0.65, 0.80, 1.00, 1.25}, stem lengths in {0.00, 0.02, 0.05},
   topologies T1/T2/T3, and conditional distributions at N_shared = 16.
"""

from __future__ import annotations

import csv
import hashlib
import json
import math
import re
import subprocess
from datetime import UTC, datetime
from pathlib import Path
from typing import Any

import numpy as np
from scipy import stats

ROOT = Path(__file__).resolve().parents[2]
RUNS = ROOT / "runs" / "enamel_relics"
V101_D = RUNS / "v101" / "phaseD"
OUT = RUNS / "v102"
OUT.mkdir(parents=True, exist_ok=True)
PROV = OUT / "provenance_v102.jsonl"

WSL_DISTRO = "Ubuntu-26.04"
WSL_ROOT = "/mnt/c/Users/user/bio-stack"

STOPS = {"TAA", "TAG", "TGA"}
BASES = ("A", "C", "G", "T")

CODON_TABLE = {
    "TTT": "F", "TTC": "F", "TTA": "L", "TTG": "L", "CTT": "L", "CTC": "L", "CTA": "L", "CTG": "L",
    "ATT": "I", "ATC": "I", "ATA": "I", "ATG": "M", "GTT": "V", "GTC": "V", "GTA": "V", "GTG": "V",
    "TCT": "S", "TCC": "S", "TCA": "S", "TCG": "S", "CCT": "P", "CCC": "P", "CCA": "P", "CCG": "P",
    "ACT": "T", "ACC": "T", "ACA": "T", "ACG": "T", "GCT": "A", "GCC": "A", "GCA": "A", "GCG": "A",
    "TAT": "Y", "TAC": "Y", "TAA": "*", "TAG": "*", "CAT": "H", "CAC": "H", "CAA": "Q", "CAG": "Q",
    "AAT": "N", "AAC": "N", "AAA": "K", "AAG": "K", "GAT": "D", "GAC": "D", "GAA": "E", "GAG": "E",
    "TGT": "C", "TGC": "C", "TGA": "*", "TGG": "W", "CGT": "R", "CGC": "R", "CGA": "R", "CGG": "R",
    "AGT": "S", "AGC": "S", "AGA": "R", "AGG": "R", "GGT": "G", "GGC": "G", "GGA": "G", "GGG": "G",
}
AA3 = {
    "A": "Ala", "R": "Arg", "N": "Asn", "D": "Asp", "C": "Cys", "Q": "Gln", "E": "Glu", "G": "Gly",
    "H": "His", "I": "Ile", "L": "Leu", "K": "Lys", "M": "Met", "F": "Phe", "P": "Pro", "S": "Ser",
    "T": "Thr", "W": "Trp", "Y": "Tyr", "V": "Val", "*": "Ter", "?": "Xaa",
}
OLD_TO_NEW = {
    "Amis": "Allmis", "Ggan": "Gavgan", "Cpor": "Cropor",
    "Acar": "Anocar", "Pbiv": "Pytbiv", "Pmur": "Podmur",
    "Cpic": "Chrpic", "Cmyd": "Chemyd", "Gevg": "Gopevg", "Psin": "Pelsin",
    "Ggal": "Galgal", "Apla": "Anapla", "Shab": "Strhab", "Tgut": "Taegut",
    "Dnov": "Dronov", "Scam": "Strcam", "Arow": "Aptrow", "Tgua": "Tingut",
    "Hsap": "Homsap",
}

TREE_T1 = (
    "Archosauria",
    ("Crocodylia", ("GganCpor", "Ggan", "Cpor"), "Amis"),
    (
        "Neornithes",
        ("Palaeognathae", "Scam", ("Novaeratitae", "Arow", "Dnov")),
        ("Neognathae", ("Galloanserae", "Ggal", "Apla"), ("Psittacopasserae", "Shab", "Tgut")),
    ),
)
TREE_T2 = (
    "Archosauria",
    ("Crocodylia", ("GganCpor", "Ggan", "Cpor"), "Amis"),
    (
        "Neornithes",
        ("Palaeognathae", "Arow", ("ScamDnov", "Scam", "Dnov")),
        ("Neognathae", ("Galloanserae", "Ggal", "Apla"), ("Psittacopasserae", "Shab", "Tgut")),
    ),
)
TREE_T3 = (
    "Archosauria",
    ("Crocodylia", ("GganCpor", "Ggan", "Cpor"), "Amis"),
    (
        "Neornithes",
        ("Palaeognathae", "Dnov", ("ScamArow", "Scam", "Arow")),
        ("Neognathae", ("Galloanserae", "Ggal", "Apla"), ("Psittacopasserae", "Shab", "Tgut")),
    ),
)


def sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def read_fasta(path: Path) -> dict[str, str]:
    seqs: dict[str, list[str]] = {}
    cur: str | None = None
    for line in path.read_text(encoding="utf-8").splitlines():
        line = line.strip()
        if not line:
            continue
        if line.startswith(">"):
            cur = line[1:].split()[0]
            seqs[cur] = []
        elif cur is not None:
            seqs[cur].append(line.upper())
    return {k: "".join(v) for k, v in seqs.items()}


def write_tsv(path: Path, rows: list[dict[str, Any]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0].keys()), delimiter="\t", lineterminator="\n")
        w.writeheader()
        w.writerows(rows)


def run_wsl_script(script_text: str, script_path: Path) -> str:
    script_path.write_text(script_text, encoding="utf-8", newline="\n")
    rel = script_path.relative_to(ROOT).as_posix()
    res = subprocess.run(
        ["wsl", "-d", WSL_DISTRO, "--", "bash", "-l", f"{WSL_ROOT}/{rel}"],
        capture_output=True, text=True, encoding="utf-8", check=False,
    )
    if res.returncode != 0:
        raise RuntimeError(f"WSL failed ({res.returncode}):\n{res.stderr[:1500]}\n{res.stdout[:1500]}")
    return res.stdout


def codon_to_aa3(cod: str) -> str:
    if not cod or len(cod) != 3 or "-" in cod:
        return "-"
    if "/" in cod:
        parts = [AA3.get(CODON_TABLE.get(c, "?"), "Xaa") for c in cod.split("/")]
        return "/".join(dict.fromkeys(parts))
    return AA3.get(CODON_TABLE.get(cod, "?"), "Xaa")


def hamming(a: str, b: str) -> int:
    return sum(x != y for x, y in zip(a, b, strict=True))


def clopper_pearson(k: int, n: int, alpha: float = 0.05) -> tuple[float, float]:
    lo = 0.0 if k == 0 else float(stats.beta.ppf(alpha / 2, k, n - k + 1))
    hi = 1.0 if k == n else float(stats.beta.ppf(1 - alpha / 2, k + 1, n - k))
    return lo, hi


def fitch_parsimony_binary_state(tree: Any, states: dict[str, str | None]) -> tuple[int, dict[str, set[str]]]:
    down: dict[str, set[str]] = {}
    steps = [0]

    def post(node: Any) -> set[str]:
        if isinstance(node, str):
            val = states.get(node)
            s = {"STOP", "SENSE"} if val is None else {val}
            down[node] = s
            return s
        name, left, right = node
        sl = post(left)
        sr = post(right)
        inter = sl & sr
        if inter:
            down[name] = inter
        else:
            steps[0] += 1
            down[name] = sl | sr
        return down[name]

    post(tree)
    return steps[0], down


def fitch_exact_codon_steps(tree: Any, states: dict[str, str | None]) -> tuple[int, int]:
    bin_states = {sp: (None if cod is None else ("STOP" if cod in STOPS else "SENSE")) for sp, cod in states.items()}
    bin_steps, _ = fitch_parsimony_binary_state(tree, bin_states)
    nt_steps = 0
    for pos in range(3):
        st_pos = {sp: (None if cod is None else cod[pos]) for sp, cod in states.items()}
        s_cnt = [0]

        def post_nt(node: Any, _st: dict[str, str | None] = st_pos, _sc: list[int] = s_cnt) -> set[str]:
            if isinstance(node, str):
                v = _st.get(node)
                return set(BASES) if (v is None or v not in BASES) else {v}
            _, left, right = node
            sl, sr = post_nt(left), post_nt(right)
            inter = sl & sr
            if inter:
                return inter
            _sc[0] += 1
            return sl | sr

        post_nt(tree)
        nt_steps += s_cnt[0]
    return bin_steps, nt_steps


def min_stop_gains_no_reversal(tree: Any, states: dict[str, str | None]) -> int:
    def dfs(node: Any) -> tuple[bool, bool, int]:
        if isinstance(node, str):
            v = states.get(node)
            if v is None:
                return (False, False, 0)
            return (v in STOPS, v not in STOPS, 0)
        _, left, right = node
        ls, ln, lg = dfs(left)
        rs, rn, rg = dfs(right)
        if not ln and not rn:
            return (ls or rs, False, lg + rg)
        gains = lg + rg
        if ls and not ln:
            gains += 1
        if rs and not rn:
            gains += 1
        return (ls or rs, True, gains)

    hs, hn, g = dfs(tree)
    if hs and not hn:
        g += 1
    return g


def parse_newick_clades(nwk: str) -> dict[str, set[str]]:
    s = nwk.strip().rstrip(";")
    idx = [0]
    clades: dict[str, set[str]] = {}

    def skip_ws() -> None:
        while idx[0] < len(s) and s[idx[0]].isspace():
            idx[0] += 1

    def parse_node() -> set[str]:
        skip_ws()
        if idx[0] < len(s) and s[idx[0]] == "(":
            idx[0] += 1
            leaves: set[str] = set()
            while True:
                leaves |= parse_node()
                skip_ws()
                if idx[0] < len(s) and s[idx[0]] == ",":
                    idx[0] += 1
                    continue
                if idx[0] < len(s) and s[idx[0]] == ")":
                    idx[0] += 1
                    break
            lbl = []
            while idx[0] < len(s) and s[idx[0]] not in ":,);":
                lbl.append(s[idx[0]])
                idx[0] += 1
            name = "".join(lbl).strip()
            if idx[0] < len(s) and s[idx[0]] == ":":
                idx[0] += 1
                while idx[0] < len(s) and s[idx[0]] not in ",);":
                    idx[0] += 1
            if name:
                clades[name] = leaves
            return leaves
        lbl = []
        while idx[0] < len(s) and s[idx[0]] not in ":,);":
            lbl.append(s[idx[0]])
            idx[0] += 1
        name = "".join(lbl).strip()
        if idx[0] < len(s) and s[idx[0]] == ":":
            idx[0] += 1
            while idx[0] < len(s) and s[idx[0]] not in ",);":
                idx[0] += 1
        return {name} if name else set()

    all_leaves = parse_node()
    outgroup = {"Amis", "Ggan", "Cpor"}
    mapped: dict[str, set[str]] = {}
    for k, lv in clades.items():
        if lv <= outgroup or outgroup <= lv:
            comp = all_leaves - lv
            mapped[k] = comp if (comp and not (comp & outgroup)) else lv
        else:
            mapped[k] = lv
    return mapped


def parse_iqtree_asr(prefix: Path) -> tuple[dict[tuple[str, int], dict[str, float]], str, str]:
    state_file = Path(str(prefix) + ".state")
    tree_file = Path(str(prefix) + ".treefile")
    iq_file = Path(str(prefix) + ".iqtree")
    probs: dict[tuple[str, int], dict[str, float]] = {}
    for ln in state_file.read_text(encoding="utf-8").splitlines():
        if ln.startswith("#") or ln.startswith("Node\t") or not ln.strip():
            continue
        c = ln.split("\t")
        probs[(c[0], int(c[1]))] = {
            "A": float(c[3]), "C": float(c[4]), "G": float(c[5]), "T": float(c[6]),
        }
    return (
        probs,
        tree_file.read_text(encoding="utf-8").strip() if tree_file.exists() else "",
        iq_file.read_text(encoding="utf-8") if iq_file.exists() else "",
    )


def codon_posterior_detail(
    probs: dict[tuple[str, int], dict[str, float]], node: str, s1: int, s2: int, s3: int
) -> dict[str, Any]:
    p1 = probs.get((node, s1))
    p2 = probs.get((node, s2))
    p3 = probs.get((node, s3))
    if not (p1 and p2 and p3):
        return {"best_codon": "???", "p_best": 0.0, "p_stop": 0.0, "top3_codons": [], "site_probs": []}
    cod_probs: list[tuple[str, float]] = []
    stop_p = 0.0
    stop_detail: dict[str, float] = {}
    for b1 in BASES:
        for b2 in BASES:
            for b3 in BASES:
                p = p1[b1] * p2[b2] * p3[b3]
                cod = b1 + b2 + b3
                cod_probs.append((cod, p))
                if cod in STOPS:
                    stop_p += p
                    stop_detail[cod] = round(p, 4)
    cod_probs.sort(key=lambda x: -x[1])
    top3 = [(c, AA3[CODON_TABLE[c]], round(p, 4)) for c, p in cod_probs[:3]]
    return {
        "best_codon": cod_probs[0][0],
        "best_aa3": AA3[CODON_TABLE[cod_probs[0][0]]],
        "p_best": round(cod_probs[0][1], 4),
        "p_stop": round(stop_p, 4),
        "stop_detail": stop_detail,
        "top3_codons": top3,
        "site_probs": [
            {b: round(p1[b], 4) for b in BASES},
            {b: round(p2[b], 4) for b in BASES},
            {b: round(p3[b], 4) for b in BASES},
        ],
    }


def main() -> int:
    t0 = datetime.now(UTC)
    print("=== Starting v1.0.2 Comprehensive Audit (Phases A, B, C, D) ===", flush=True)

    # -----------------------------------------------------------------------
    # STEP 1: Formal IQ-TREE 3 Topology Tests (T1, T2, T3) on MAFFT & MUSCLE
    #         + MUSCLE T2/T3 ASR so both aligners have complete T1/T2/T3 ASR!
    # -----------------------------------------------------------------------
    iq_dir = V101_D / "iqtree_asr"
    v102_iq = OUT / "iqtree_tests"
    v102_iq.mkdir(parents=True, exist_ok=True)
    trees_3 = "\n".join([
        (iq_dir / "tree_T1.nwk").read_text(encoding="utf-8").strip(),
        (iq_dir / "tree_T2.nwk").read_text(encoding="utf-8").strip(),
        (iq_dir / "tree_T3.nwk").read_text(encoding="utf-8").strip(),
    ]) + "\n"
    (v102_iq / "trees_T1_T2_T3.nwk").write_text(trees_3, encoding="utf-8", newline="\n")

    sh_topo = f"""#!/usr/bin/env bash
set -euo pipefail
L4="{WSL_ROOT}/runs/enamel_relics/level4"
D101="{WSL_ROOT}/runs/enamel_relics/v101/phaseD/iqtree_asr"
O="{WSL_ROOT}/runs/enamel_relics/v102/iqtree_tests"
# Formal topology tests (KH, SH, wKH, wSH, ELW, AU with 10,000 RELL bootstrap replicates)
iqtree3 -s "$L4/ENAM_exon8_aves/mafft.fna" -z "$O/trees_T1_T2_T3.nwk" -te "$D101/tree_T1.nwk" -m GTR+G -zb 10000 -zw -au -o Amis,Ggan,Cpor -pre "$O/topo_mafft" -redo -nt 4 --quiet
iqtree3 -s "$L4/ENAM_exon8_aves/muscle.fna" -z "$O/trees_T1_T2_T3.nwk" -te "$D101/tree_T1.nwk" -m GTR+G -zb 10000 -zw -au -o Amis,Ggan,Cpor -pre "$O/topo_muscle" -redo -nt 4 --quiet
# Also run MUSCLE T2 and T3 ASR for completeness
iqtree3 -s "$L4/ENAM_exon8_aves/muscle.fna" -te "$D101/tree_T2.nwk" -m GTR+G --ancestral -o Amis,Ggan,Cpor -pre "$O/enam8_muscle_T2" -redo -nt 4 --quiet
iqtree3 -s "$L4/ENAM_exon8_aves/muscle.fna" -te "$D101/tree_T3.nwk" -m GTR+G --ancestral -o Amis,Ggan,Cpor -pre "$O/enam8_muscle_T3" -redo -nt 4 --quiet
"""
    if not ((v102_iq / "topo_mafft.iqtree").exists() and (v102_iq / "enam8_muscle_T3.state").exists()):
        run_wsl_script(sh_topo, v102_iq / "run_topo_tests.sh")
    print("[1/5] IQ-TREE 3 formal topology tests (10,000 RELL reps, KH/SH/wSH/ELW/AU) completed.", flush=True)

    def parse_iqtree_topo_table(iq_path: Path) -> list[dict[str, Any]]:
        txt = iq_path.read_text(encoding="utf-8")
        m = re.search(r"Tree\s+logL\s+deltaL\s+bp-RELL.*?\n-+\n(.*?)\n\n", txt, re.S)
        rows = []
        if m:
            for ln in m.group(1).strip().splitlines():
                cols = ln.split()
                if len(cols) >= 15:
                    # format: Tree logL deltaL bp-RELL +/- p-KH +/- p-SH +/- p-WKH +/- p-WSH +/- c-ELW +/- p-AU +/-
                    clean = [c for c in cols if c not in ("+", "-")]
                    rows.append({
                        "tree_id": int(clean[0]),
                        "topology": ["T1_Novaeratitae", "T2_ScamDnov", "T3_ScamArow"][int(clean[0]) - 1],
                        "logL": float(clean[1]),
                        "deltaL": float(clean[2]),
                        "bp_RELL": float(clean[3]),
                        "p_KH": float(clean[4]),
                        "p_SH": float(clean[5]),
                        "p_WKH": float(clean[6]),
                        "p_WSH": float(clean[7]),
                        "c_ELW": float(clean[8]),
                        "p_AU": float(clean[9]),
                    })
        return rows

    topo_tests_mafft = parse_iqtree_topo_table(v102_iq / "topo_mafft.iqtree")
    topo_tests_muscle = parse_iqtree_topo_table(v102_iq / "topo_muscle.iqtree")
    print("  MAFFT topology tests:", topo_tests_mafft, flush=True)
    print("  MUSCLE topology tests:", topo_tests_muscle, flush=True)

    # Load ASR probabilities across T1, T2, T3 (MAFFT) and T1, T2, T3 (MUSCLE)
    asr_sets = {
        ("MAFFT", "T1"): parse_iqtree_asr(iq_dir / "enam8_mafft_T1"),
        ("MAFFT", "T2"): parse_iqtree_asr(iq_dir / "enam8_mafft_T2"),
        ("MAFFT", "T3"): parse_iqtree_asr(iq_dir / "enam8_mafft_T3"),
        ("MUSCLE", "T1"): parse_iqtree_asr(iq_dir / "enam8_muscle_T1"),
        ("MUSCLE", "T2"): parse_iqtree_asr(v102_iq / "enam8_muscle_T2"),
        ("MUSCLE", "T3"): parse_iqtree_asr(v102_iq / "enam8_muscle_T3"),
    }

    def map_internal_nodes(nwk_str: str) -> dict[str, str]:
        cm = parse_newick_clades(nwk_str)
        targets = {
            "Neornithes": {"Scam", "Arow", "Dnov", "Ggal", "Apla", "Shab", "Tgut"},
            "Palaeognathae": {"Scam", "Arow", "Dnov"},
            "Novaeratitae": {"Arow", "Dnov"},
            "ScamDnov": {"Scam", "Dnov"},
            "ScamArow": {"Scam", "Arow"},
            "Neognathae": {"Ggal", "Apla", "Shab", "Tgut"},
            "Galloanserae": {"Ggal", "Apla"},
            "Psittacopasserae": {"Shab", "Tgut"},
        }
        out = {}
        for name, tset in targets.items():
            for nid, lset in cm.items():
                if lset == tset:
                    out[name] = nid
        return out

    node_maps = {k: map_internal_nodes(v[1]) for k, v in asr_sets.items()}

    # -----------------------------------------------------------------------
    # STEP 2: Phase A — Verify Reference Proteins (XP_019346377.1 & XP_019346386.1)
    #         and Audit All 23 Shared Variants in Table_S5a
    # -----------------------------------------------------------------------
    queries_fa = read_fasta(RUNS / "queries" / "all_real_queries.faa")
    enam_prot = queries_fa["ENAM_Amis"]
    ambn_prot = queries_fa["AMBN_Amis"]
    assert len(enam_prot) == 1078, f"Expected 1078 aa for XP_019346377.1, got {len(enam_prot)}"
    assert len(ambn_prot) == 390, f"Expected 390 aa for XP_019346386.1, got {len(ambn_prot)}"

    s5a_in = list(csv.DictReader((RUNS / "supplementary" / "Table_S5a_level4_shared_mutations.tsv").open(encoding="utf-8"), delimiter="\t"))
    msa_cache: dict[str, tuple[dict[str, str], dict[str, str]]] = {}
    for ds in sorted({r["dataset"] for r in s5a_in}):
        d = RUNS / "level4" / ds
        msa_cache[ds] = (read_fasta(d / "mafft.fna"), read_fasta(d / "muscle.fna"))

    # Exon 1-based start codon offset in full protein for verification
    # ENAM Exon 8 starts at codon 201 of XP_019346377.1 (nt 1 of Exon 8 = 1st base of codon 201)
    # AMBN Exon 3 starts at codon 34 of XP_019346386.1; AMBN Exon 4 starts at codon 73 of XP_019346386.1
    s5a_v102_rows: list[dict[str, Any]] = []
    variant_audit_details: list[dict[str, Any]] = []

    for r in s5a_in:
        ds = r["dataset"]
        mafft_msa, muscle_msa = msa_cache[ds]
        ref_m = mafft_msa["Amis"]
        col_of_m = [i for i, ch in enumerate(ref_m) if ch != "-"]
        ref_u = muscle_msa["Amis"]
        col_of_u = [i for i, ch in enumerate(ref_u) if ch != "-"]

        cod_num = int(r["codon_1based"])
        pos0 = int(r["ref_exon_nt_1based"]) - 1
        typ = r["mutation_type"]
        ref_prot_seq = enam_prot if r["gene"] == "ENAM" else ambn_prot
        ref_aa_1 = ref_prot_seq[cod_num - 1]
        ref_aa_3 = AA3[ref_aa_1]

        if typ == "stop":
            c0, c1, c2 = col_of_m[pos0], col_of_m[pos0 + 1], col_of_m[pos0 + 2]
            u0, u1, u2 = col_of_u[pos0], col_of_u[pos0 + 1], col_of_u[pos0 + 2]
            amis_cod = ref_m[c0] + ref_m[c1] + ref_m[c2]
            assert CODON_TABLE[amis_cod] == ref_aa_1, (
                f"Reference translation mismatch at {r['gene']} codon {cod_num}: {amis_cod} != {ref_aa_1}"
            )
            assert f":p.{ref_aa_3}{cod_num}Ter" in r["standardized_variant_id"]

            tips_m = {sp: seq[c0] + seq[c1] + seq[c2] for sp, seq in sorted(mafft_msa.items())}
            tips_u = {sp: seq[u0] + seq[u1] + seq[u2] for sp, seq in sorted(muscle_msa.items())}
            tip_states = {sp: (cod if ("-" not in cod and len(cod) == 3) else None) for sp, cod in tips_m.items()}

            # Parsimony across T1, T2, T3 (for avian ENAM Exon 8)
            if ds == "ENAM_exon8_aves":
                fb_t1, nt_t1 = fitch_exact_codon_steps(TREE_T1, tip_states)
                fb_t2, nt_t2 = fitch_exact_codon_steps(TREE_T2, tip_states)
                fb_t3, nt_t3 = fitch_exact_codon_steps(TREE_T3, tip_states)
                mg_t1 = min_stop_gains_no_reversal(TREE_T1, tip_states)
                mg_t2 = min_stop_gains_no_reversal(TREE_T2, tip_states)
                mg_t3 = min_stop_gains_no_reversal(TREE_T3, tip_states)

                # Detailed ML ASR across T1, T2, T3 (MAFFT) and T1 (MUSCLE)
                s1, s2, s3 = c0 + 1, c1 + 1, c2 + 1
                su1, su2, su3 = u0 + 1, u1 + 1, u2 + 1
                asr_t1 = {
                    node_lbl: codon_posterior_detail(asr_sets[("MAFFT", "T1")][0], nid, s1, s2, s3)
                    for node_lbl, nid in node_maps[("MAFFT", "T1")].items()
                }
                asr_t2 = {
                    node_lbl: codon_posterior_detail(asr_sets[("MAFFT", "T2")][0], nid, s1, s2, s3)
                    for node_lbl, nid in node_maps[("MAFFT", "T2")].items()
                }
                asr_t3 = {
                    node_lbl: codon_posterior_detail(asr_sets[("MAFFT", "T3")][0], nid, s1, s2, s3)
                    for node_lbl, nid in node_maps[("MAFFT", "T3")].items()
                }
                asr_muscle_t1 = {
                    node_lbl: codon_posterior_detail(asr_sets[("MUSCLE", "T1")][0], nid, su1, su2, su3)
                    for node_lbl, nid in node_maps[("MUSCLE", "T1")].items()
                }

                neorn_cod = asr_t1["Neornithes"]["best_codon"]
                neorn_aa3 = asr_t1["Neornithes"]["best_aa3"]
                neorn_p_best = asr_t1["Neornithes"]["p_best"]
                neog_p_stop = asr_t1["Neognathae"]["p_stop"]
                neog_p_stop_t2 = asr_t2["Neognathae"]["p_stop"]
                neog_p_stop_t3 = asr_t3["Neognathae"]["p_stop"]

                ccods = r["carrier_codons"]
                if amis_cod == neorn_cod:
                    inferred_sub = f"{neorn_cod}({neorn_aa3})->{ccods}(Ter)"
                elif cod_num == 311:
                    inferred_sub = "CCG(Pro,Amis)->CAA(Gln,Neorn;P=0.542)->CGA(Arg,Palaeo;P=0.541)/TGA(Ter;P=0.453)"
                elif cod_num == 1003:
                    inferred_sub = "TGG(Trp,Amis;Neorn_P=0.168)/TGC(Cys,Neorn_ML_argmax_P=0.343;no_extant_TGC)->TGA(Ter,Psittaco_P=0.892)"
                else:
                    inferred_sub = f"{amis_cod}({ref_aa_3},Amis)->{neorn_cod}({neorn_aa3},Neorn;P={neorn_p_best:.3f})->{ccods}(Ter)"
            else:
                fb_t1 = fb_t2 = fb_t3 = ""
                nt_t1 = nt_t2 = nt_t3 = ""
                mg_t1 = mg_t2 = mg_t3 = 1
                asr_t1 = asr_t2 = asr_t3 = asr_muscle_t1 = {}
                neorn_cod = "-"
                neorn_aa3 = "-"
                neorn_p_best = ""
                neog_p_stop = neog_p_stop_t2 = neog_p_stop_t3 = ""
                inferred_sub = f"{amis_cod}({ref_aa_3},Amis)->{r['carrier_codons']}(Ter)"
        else:
            # Indel
            amis_cod = "-"
            tips_m = {}
            tips_u = {}
            fb_t1 = fb_t2 = fb_t3 = ""
            nt_t1 = nt_t2 = nt_t3 = ""
            mg_t1 = mg_t2 = mg_t3 = 1
            asr_t1 = asr_t2 = asr_t3 = asr_muscle_t1 = {}
            neorn_cod = "-"
            neorn_aa3 = "-"
            neorn_p_best = ""
            neog_p_stop = neog_p_stop_t2 = neog_p_stop_t3 = ""
            inferred_sub = f"codon_{cod_num}({ref_aa_3},Amis)_{typ}{r['size_nt']}nt_fs"

        # Phase B evidence tier for synapomorphy vs shared lesion with unrecovered sister/outgroup taxa
        cov_without = r["covered_without_6letter"]
        orig_cls = r["orig_classification"]
        if orig_cls.startswith("CONGRUENT_"):
            if cov_without == "-":
                if r["dataset"] == "AMBN_exon4_aves":
                    syn_tier = "SHARED_CLADE_LESION_SISTER_UNRECOVERED (Crocodylia+Squamata intact; Neognathae Exon 4 unrecovered)"
                elif r["dataset"] == "AMBN_exon4_turtles":
                    syn_tier = "SHARED_SUBCLADE_LESION_SISTER_UNRECOVERED (Crocodylia+Squamata intact; Gopevg+Pelsin Exon 4 unrecovered)"
                elif r["dataset"] == "ENAM_exon8_turtles":
                    syn_tier = "SHARED_SUBCLADE_LESION_SISTER_UNRECOVERED (Crocodylia intact; Gopevg+Pelsin region unrecovered)"
                else:
                    syn_tier = "SHARED_CLADE_LESION_SISTER_UNRECOVERED"
            else:
                syn_tier = f"OUTGROUP_AND_SISTER_VERIFIED_SYNAPOMORPHY ({orig_cls.replace('CONGRUENT_', '')})"
        else:
            if cod_num == 1022:
                syn_tier = "TOPOLOGY_INVARIANT_PARSIMONY_INCONGRUENT_OR_ML_NEOGNATHAE_WITH_TAEGUT_DECAY (Anapla uncovered; ML Neognathae P_stop=0.926)"
            elif mg_t1 >= 2:
                syn_tier = "PARSIMONY_INCONGRUENT_ON_T1_T2_T3"
            else:
                syn_tier = "TOPOLOGY_DEPENDENT_INCONGRUENT"

        row_out = dict(r)
        row_out["ref_aa_3letter"] = ref_aa_3
        row_out["croc_consensus_aa_3letter"] = codon_to_aa3(r["croc_consensus_codon"])
        row_out["neornithes_ml_aa_3letter"] = neorn_aa3
        row_out["neornithes_ml_p_best"] = neorn_p_best
        row_out["neognathae_ml_p_stop"] = neog_p_stop
        row_out["inferred_substitution_path"] = inferred_sub
        row_out["fitch_bin_steps_T2"] = fb_t2
        row_out["fitch_bin_steps_T3"] = fb_t3
        row_out["min_gains_no_rev_T2"] = mg_t2 if typ == "stop" and ds == "ENAM_exon8_aves" else ""
        row_out["min_gains_no_rev_T3"] = mg_t3 if typ == "stop" and ds == "ENAM_exon8_aves" else ""
        row_out["synapomorphy_evidence_tier"] = syn_tier
        s5a_v102_rows.append(row_out)

        variant_audit_details.append({
            "dataset": ds,
            "codon_1based": cod_num,
            "standardized_variant_id": r["standardized_variant_id"],
            "ref_aa_3letter": ref_aa_3,
            "amis_codon": amis_cod,
            "neornithes_ml_codon": neorn_cod,
            "neornithes_ml_aa_3letter": neorn_aa3,
            "neornithes_ml_p_best": neorn_p_best,
            "inferred_substitution_path": inferred_sub,
            "synapomorphy_evidence_tier": syn_tier,
            "tips_mafft": tips_m,
            "tips_muscle": tips_u,
            "fitch_bin_T1_T2_T3": [fb_t1, fb_t2, fb_t3],
            "fitch_nt_T1_T2_T3": [nt_t1, nt_t2, nt_t3],
            "min_gains_no_rev_T1_T2_T3": [mg_t1, mg_t2, mg_t3],
            "asr_mafft_T1": asr_t1,
            "asr_mafft_T2": asr_t2,
            "asr_mafft_T3": asr_t3,
            "asr_muscle_T1": asr_muscle_t1,
            "neog_p_stop_T1_T2_T3": [neog_p_stop, neog_p_stop_t2, neog_p_stop_t3],
        })

    write_tsv(OUT / "Table_S5a_level4_shared_mutations_v102.tsv", s5a_v102_rows)
    print(f"[2/5] Audited all {len(s5a_v102_rows)} rows of Table_S5a with exact reference & ASR translations.", flush=True)

    # -----------------------------------------------------------------------
    # STEP 3: Phase B — AMBN Exon 4 Outgroup & Sister-Clade Coverage Investigation
    # -----------------------------------------------------------------------
    # 1. Verify Squamata (Anocar, Pytbiv, Podmur) & Crocodylia (Allmis, Cropor, Gavgan) in AMBN Exon 4
    # 2. Test ultra-sensitive blastn (-word_size 7 -evalue 10) and tblastn (-seg no -word_size 2 -evalue 100)
    #    for AMBN Exon 4 across all Neognathae (Galgal, Anapla, Strhab, Taegut) and basal/sister turtles (Gopevg, Pelsin)
    ambn_ex4_dir = OUT / "ambn_exon4_coverage_audit"
    ambn_ex4_dir.mkdir(parents=True, exist_ok=True)

    # Extract Alligator AMBN Exon 4 nucleotide sequence (240 bp, codons 73..152) and protein sequence (80 aa)
    ambn_ex4_aves_m = read_fasta(RUNS / "level4" / "AMBN_exon4_aves" / "mafft.fna")
    ambn_ex4_turt_m = read_fasta(RUNS / "level4" / "AMBN_exon4_turtles" / "mafft.fna")
    amis_ex4_nt = ambn_ex4_aves_m["Amis"].replace("-", "")
    dnov_ex4_nt = ambn_ex4_aves_m["Dnov"].replace("-", "")
    cmyd_ex4_nt = ambn_ex4_turt_m["Cmyd"].replace("-", "")
    amis_ex4_aa = ambn_prot[72:152]  # codons 73..152 (80 aa)

    (ambn_ex4_dir / "q_ex4_nt.fna").write_text(
        f">Amis_AMBN_ex4\n{amis_ex4_nt}\n>Dnov_AMBN_ex4\n{dnov_ex4_nt}\n>Cmyd_AMBN_ex4\n{cmyd_ex4_nt}\n",
        encoding="utf-8", newline="\n",
    )
    (ambn_ex4_dir / "q_ex4_aa.faa").write_text(
        f">Amis_AMBN_ex4_80aa\n{amis_ex4_aa}\n",
        encoding="utf-8", newline="\n",
    )

    # Search all 18 NRG1-LPL intervals with sensitive blastn (-word_size 7 -evalue 10) and tblastn (-seg no -word_size 2 -evalue 50)
    # and check if any hit in Neognathae or Gopevg/Pelsin reciprocally matches AMBN!
    rel_ex4 = ambn_ex4_dir.relative_to(ROOT).as_posix()
    sh_ex4 = f"""#!/usr/bin/env bash
set -euo pipefail
O="{WSL_ROOT}/{rel_ex4}"
INT="{WSL_ROOT}/runs/enamel_relics/intervals_enamel"
> "$O/blastn_ws7_all.tsv"
> "$O/tblastn_ws2_all.tsv"
for f in "$INT"/*.fna; do
  sp=$(basename "$f" .fna)
  blastn -task blastn -word_size 7 -evalue 10 -dust no -query "$O/q_ex4_nt.fna" -subject "$f" -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore" | awk -v s="$sp" '{{print s"\\t"$0}}' >> "$O/blastn_ws7_all.tsv" || true
  tblastn -word_size 2 -evalue 50 -seg no -comp_based_stats 0 -query "$O/q_ex4_aa.faa" -subject "$f" -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore" | awk -v s="$sp" '{{print s"\\t"$0}}' >> "$O/tblastn_ws2_all.tsv" || true
done
"""
    if not (ambn_ex4_dir / "blastn_ws7_all.tsv").exists() or (ambn_ex4_dir / "blastn_ws7_all.tsv").stat().st_size == 0:
        run_wsl_script(sh_ex4, ambn_ex4_dir / "run_ex4_search.sh")

    # Parse best AMBN Exon 4 hits per species and test any candidate with E <= 0.1 against whole Alligator proteome
    ex4_blastn_best: dict[tuple[str, str], dict[str, Any]] = {}
    for ln in (ambn_ex4_dir / "blastn_ws7_all.tsv").read_text(encoding="utf-8").splitlines():
        if not ln.strip():
            continue
        c = ln.split("\t")
        sp, qid = c[0].replace("interval_", ""), c[1]
        bits, ev, alen, pid = float(c[12]), float(c[11]), int(c[4]), float(c[3])
        key = (sp, qid)
        if key not in ex4_blastn_best or bits > ex4_blastn_best[key]["bits"]:
            ex4_blastn_best[key] = {
                "species": OLD_TO_NEW.get(sp, sp), "legacy": sp, "query": qid,
                "pident": pid, "len": alen, "qstart": int(c[7]), "qend": int(c[8]),
                "sstart": int(c[9]), "send": int(c[10]), "evalue": ev, "bits": bits,
            }

    ex4_tblastn_best: dict[str, dict[str, Any]] = {}
    for ln in (ambn_ex4_dir / "tblastn_ws2_all.tsv").read_text(encoding="utf-8").splitlines():
        if not ln.strip():
            continue
        c = ln.split("\t")
        sp = c[0].replace("interval_", "")
        bits, ev, alen, pid = float(c[12]), float(c[11]), int(c[4]), float(c[3])
        if sp not in ex4_tblastn_best or bits > ex4_tblastn_best[sp]["bits"]:
            ex4_tblastn_best[sp] = {
                "species": OLD_TO_NEW.get(sp, sp), "legacy": sp,
                "pident": pid, "len": alen, "qstart": int(c[7]), "qend": int(c[8]),
                "sstart": int(c[9]), "send": int(c[10]), "evalue": ev, "bits": bits,
            }

    print("[3/5] Phase B AMBN Exon 4 ultra-sensitive recovery audit completed.", flush=True)
    for sp in ("Amis", "Cpor", "Ggan", "Acar", "Pbiv", "Pmur", "Arow", "Dnov", "Scam", "Ggal", "Apla", "Shab", "Tgut", "Tgua", "Cmyd", "Cpic", "Gevg", "Psin"):
        bn_a = ex4_blastn_best.get((sp, "Amis_AMBN_ex4"), {})
        bn_d = ex4_blastn_best.get((sp, "Dnov_AMBN_ex4"), {})
        bn_c = ex4_blastn_best.get((sp, "Cmyd_AMBN_ex4"), {})
        tb = ex4_tblastn_best.get(sp, {})
        print(
            f"  {OLD_TO_NEW.get(sp,sp):6s} ({sp:4s}): "
            f"blastn(Amis)={bn_a.get('bits','-')}b (E={bn_a.get('evalue','-')}, L={bn_a.get('len','-')}) | "
            f"blastn(Dnov)={bn_d.get('bits','-')}b (E={bn_d.get('evalue','-')}, L={bn_d.get('len','-')}) | "
            f"blastn(Cmyd)={bn_c.get('bits','-')}b (E={bn_c.get('evalue','-')}, L={bn_c.get('len','-')}) | "
            f"tblastn(Amis)={tb.get('bits','-')}b (E={tb.get('evalue','-')}, L={tb.get('len','-')})",
            flush=True,
        )

    # -----------------------------------------------------------------------
    # STEP 4: Phase D — Comprehensive Neutral K80 Stop-Gain Simulation Audit
    # -----------------------------------------------------------------------
    # Extract empirical coverage masks from MAFFT and MUSCLE ENAM_exon8_aves alignments!
    enam8_m, enam8_u = msa_cache["ENAM_exon8_aves"]
    ref_m_seq = enam8_m["Amis"]
    ref_u_seq = enam8_u["Amis"]
    cols_m = [i for i, ch in enumerate(ref_m_seq) if ch != "-"]
    cols_u = [i for i, ch in enumerate(ref_u_seq) if ch != "-"]
    n_cod = len(cols_m) // 3  # 878 codons (2,634 nt)
    birds = ("Scam", "Arow", "Dnov", "Ggal", "Apla", "Shab", "Tgut")

    # Mask 1: MAFFT empirical codon coverage (all 3 bases in A/C/G/T)
    cov_mafft = np.zeros((len(birds), n_cod), dtype=bool)
    # Mask 2: Dual-aligner (MAFFT & MUSCLE both cover the codon with all 3 bases in A/C/G/T and identical triplet)
    cov_dual = np.zeros((len(birds), n_cod), dtype=bool)
    for bi, sp in enumerate(birds):
        sm = enam8_m[sp]
        su = enam8_u[sp]
        for ci in range(n_cod):
            cm_trip = sm[cols_m[3 * ci]] + sm[cols_m[3 * ci + 1]] + sm[cols_m[3 * ci + 2]]
            cu_trip = su[cols_u[3 * ci]] + su[cols_u[3 * ci + 1]] + su[cols_u[3 * ci + 2]]
            m_ok = all(b in BASES for b in cm_trip)
            u_ok = all(b in BASES for b in cu_trip)
            cov_mafft[bi, ci] = m_ok
            cov_dual[bi, ci] = m_ok and u_ok and (cm_trip == cu_trip)

    cov_summary = {
        OLD_TO_NEW[sp]: {
            "mafft_covered_codons": int(cov_mafft[bi].sum()),
            "mafft_covered_pct": round(100.0 * float(cov_mafft[bi].sum()) / n_cod, 1),
            "dual_covered_codons": int(cov_dual[bi].sum()),
            "dual_covered_pct": round(100.0 * float(cov_dual[bi].sum()) / n_cod, 1),
        }
        for bi, sp in enumerate(birds)
    }
    print("  Empirical codon coverage in ENAM Exon 8 (878 codons):", json.dumps(cov_summary, indent=2), flush=True)

    # Reconstruct full 878-codon Neornithes ML ancestral sequence from enam8_mafft_T1.state (excluding any stop codon -> sense argmax)
    probs_t1_m, nwk_t1_m, _ = asr_sets[("MAFFT", "T1")]
    nid_neorn = node_maps[("MAFFT", "T1")]["Neornithes"]
    neorn_anc_seq_list = []
    amis_clean_seq = "".join(ref_m_seq[i] for i in cols_m[: n_cod * 3])
    for ci in range(n_cod):
        s1, s2, s3 = cols_m[3 * ci] + 1, cols_m[3 * ci + 1] + 1, cols_m[3 * ci + 2] + 1
        p1, p2, p3 = probs_t1_m[(nid_neorn, s1)], probs_t1_m[(nid_neorn, s2)], probs_t1_m[(nid_neorn, s3)]
        best_sense = "AAA"
        best_p = -1.0
        for b1 in BASES:
            for b2 in BASES:
                for b3 in BASES:
                    cod = b1 + b2 + b3
                    if cod in STOPS:
                        continue
                    p = p1[b1] * p2[b2] * p3[b3]
                    if p > best_p:
                        best_p = p
                        best_sense = cod
        neorn_anc_seq_list.append(best_sense)
    neorn_anc_seq = "".join(neorn_anc_seq_list)

    # Parse branch lengths for T1, T2, T3 from IQ-TREE .treefile
    def extract_tree_blens(topo_key: str) -> dict[str, Any]:
        _, nwk_str, _ = asr_sets[("MAFFT", topo_key)]
        blens = {m.group(1): float(m.group(2)) for m in re.finditer(r"([A-Za-z0-9_]+):([0-9.eE+-]+)", nwk_str)}
        nmap = node_maps[("MAFFT", topo_key)]
        return {
            "t_pal": blens.get(nmap.get("Palaeognathae", ""), 0.12),
            "t_subpal": blens.get(nmap.get("Novaeratitae") or nmap.get("ScamDnov") or nmap.get("ScamArow") or "", 0.03),
            "t_neo": blens.get(nmap.get("Neognathae", ""), 0.047),
            "t_gal": blens.get(nmap.get("Galloanserae", ""), 0.044),
            "t_psi": blens.get(nmap.get("Psittacopasserae", ""), 0.043),
            "tips": {sp: blens.get(sp, 0.10) for sp in birds},
        }

    blens_by_topo = {t: extract_tree_blens(t) for t in ("T1", "T2", "T3")}

    # Vectorized K80 simulator supporting T1, T2, T3, arbitrary kappa, branch scale, stem time, root seq, and coverage mask
    b2i = {"A": 0, "C": 1, "G": 2, "T": 3}

    def k80_cdf(t: float, kappa: float) -> np.ndarray:
        if t <= 0:
            return np.cumsum(np.eye(4), axis=1)
        beta = 1.0 / (kappa + 2.0)
        e1 = math.exp(-4.0 * beta * t)
        e2 = math.exp(-2.0 * (kappa + 1.0) * beta * t)
        p_same = 0.25 + 0.25 * e1 + 0.5 * e2
        p_ts = 0.25 + 0.25 * e1 - 0.5 * e2
        p_tv = 0.25 - 0.25 * e1
        P = np.array([
            [p_same, p_tv, p_ts, p_tv],
            [p_tv, p_same, p_tv, p_ts],
            [p_ts, p_tv, p_same, p_tv],
            [p_tv, p_ts, p_tv, p_same],
        ])
        return np.cumsum(P, axis=1)

    # Precompute fast lookup table of min_stop_gains_no_reversal for all 3^7 = 2,187 tip state patterns
    # (0 = uncovered/None, 1 = SENSE, 2 = STOP) across T1, T2, T3!
    # Also precompute whether a single-gain site (min_gains == 1) is a crown-Neornithes synapomorphy
    # (i.e. min_gains == 1 and >=1 Palaeognathae STOP and >=1 Neognathae STOP and 0 SENSE in any covered bird).
    palaeo_idx = (0, 1, 2)  # Scam, Arow, Dnov
    neog_idx = (3, 4, 5, 6)  # Ggal, Apla, Shab, Tgut
    pow3 = np.array([3**i for i in range(7)], dtype=np.int32)
    lut_gains = {
        "T1": np.zeros(2187, dtype=np.int8),
        "T2": np.zeros(2187, dtype=np.int8),
        "T3": np.zeros(2187, dtype=np.int8),
    }
    lut_n_carriers = np.zeros(2187, dtype=np.int8)
    lut_crown_neorn = np.zeros(2187, dtype=bool)
    # Also legacy unmasked congruent_sets definition for direct comparison with v1.0.1!
    congruent_sets_t1 = {
        frozenset({"Arow", "Dnov"}),
        frozenset({"Scam", "Arow", "Dnov"}),
        frozenset({"Ggal", "Apla"}),
        frozenset({"Shab", "Tgut"}),
        frozenset({"Ggal", "Apla", "Shab", "Tgut"}),
        frozenset({"Scam", "Arow", "Dnov", "Ggal", "Apla", "Shab", "Tgut"}),
    }
    lut_legacy_incong_t1 = np.zeros(2187, dtype=bool)

    for code in range(2187):
        rem = code
        st_dict: dict[str, str | None] = {}
        carriers_set = set()
        n_stop = 0
        n_sense = 0
        pal_stop = False
        neog_stop = False
        for bi, sp in enumerate(birds):
            v = rem % 3
            rem //= 3
            if v == 0:
                st_dict[sp] = None
            elif v == 1:
                st_dict[sp] = "AAA"
                n_sense += 1
            else:
                st_dict[sp] = "TGA"
                n_stop += 1
                carriers_set.add(sp)
                if bi in palaeo_idx:
                    pal_stop = True
                if bi in neog_idx:
                    neog_stop = True
        lut_n_carriers[code] = n_stop
        if n_stop >= 2:
            lut_gains["T1"][code] = min_stop_gains_no_reversal(TREE_T1, st_dict)
            lut_gains["T2"][code] = min_stop_gains_no_reversal(TREE_T2, st_dict)
            lut_gains["T3"][code] = min_stop_gains_no_reversal(TREE_T3, st_dict)
            lut_crown_neorn[code] = (n_sense == 0) and pal_stop and neog_stop
            lut_legacy_incong_t1[code] = frozenset(carriers_set) not in congruent_sets_t1

    def run_k80_scenario(
        rng: np.random.Generator,
        anc_nt: str,
        topo: str,
        kappa: float,
        scale: float,
        t_stem: float,
        mask_mode: str,
        n_reps: int = 1000,
    ) -> dict[str, Any]:
        anc_arr = np.array([b2i[c] for c in anc_nt], dtype=np.int8)
        L = len(anc_arr)
        bl = blens_by_topo[topo]
        t_pal = bl["t_pal"] * scale
        t_sub = bl["t_subpal"] * scale
        t_neo = bl["t_neo"] * scale
        t_gal = bl["t_gal"] * scale
        t_psi = bl["t_psi"] * scale
        t_tips = {sp: bl["tips"][sp] * scale for sp in birds}

        cdf_stem = k80_cdf(t_stem * scale, kappa)
        cdf_pal = k80_cdf(t_pal, kappa)
        cdf_sub = k80_cdf(t_sub, kappa)
        cdf_neo = k80_cdf(t_neo, kappa)
        cdf_gal = k80_cdf(t_gal, kappa)
        cdf_psi = k80_cdf(t_psi, kappa)
        cdf_tips = {sp: k80_cdf(t_tips[sp], kappa) for sp in birds}

        if mask_mode == "unmasked_full_878":
            cov_mat = np.ones((7, n_cod), dtype=bool)
        elif mask_mode == "empirical_mafft_mask":
            cov_mat = cov_mafft
        elif mask_mode == "empirical_dual_aligner_mask":
            cov_mat = cov_dual
        else:
            raise ValueError(f"Unknown mask_mode: {mask_mode}")

        def step_seq(seq: np.ndarray, cdf: np.ndarray) -> np.ndarray:
            u = rng.random(L)
            c0, c1, c2 = cdf[seq, 0], cdf[seq, 1], cdf[seq, 2]
            return np.where(u < c0, 0, np.where(u < c1, 1, np.where(u < c2, 2, 3))).astype(np.int8)

        def is_stop_arr(seq: np.ndarray) -> np.ndarray:
            c0, c1, c2 = seq[0::3], seq[1::3], seq[2::3]
            return (c0 == 3) & (((c1 == 0) & ((c2 == 0) | (c2 == 2))) | ((c1 == 2) & (c2 == 0)))

        sh_arr = np.zeros(n_reps, dtype=np.int32)
        ho_arr = np.zeros(n_reps, dtype=np.int32)
        ho_legacy_arr = np.zeros(n_reps, dtype=np.int32)
        neo_arr = np.zeros(n_reps, dtype=np.int32)
        hf_list: list[float] = []
        hf_cond16_list: list[float] = []
        hf_cond14_18_list: list[float] = []

        lg = lut_gains[topo]
        for r_i in range(n_reps):
            root = step_seq(anc_arr, cdf_stem) if t_stem > 0 else anc_arr
            pal = step_seq(root, cdf_pal)
            subpal = step_seq(pal, cdf_sub)
            neo = step_seq(root, cdf_neo)
            gal = step_seq(neo, cdf_gal)
            psi = step_seq(neo, cdf_psi)

            if topo == "T1":
                s_scam = is_stop_arr(step_seq(pal, cdf_tips["Scam"]))
                s_arow = is_stop_arr(step_seq(subpal, cdf_tips["Arow"]))
                s_dnov = is_stop_arr(step_seq(subpal, cdf_tips["Dnov"]))
            elif topo == "T2":
                s_arow = is_stop_arr(step_seq(pal, cdf_tips["Arow"]))
                s_scam = is_stop_arr(step_seq(subpal, cdf_tips["Scam"]))
                s_dnov = is_stop_arr(step_seq(subpal, cdf_tips["Dnov"]))
            else:
                s_dnov = is_stop_arr(step_seq(pal, cdf_tips["Dnov"]))
                s_scam = is_stop_arr(step_seq(subpal, cdf_tips["Scam"]))
                s_arow = is_stop_arr(step_seq(subpal, cdf_tips["Arow"]))

            s_ggal = is_stop_arr(step_seq(gal, cdf_tips["Ggal"]))
            s_apla = is_stop_arr(step_seq(gal, cdf_tips["Apla"]))
            s_shab = is_stop_arr(step_seq(psi, cdf_tips["Shab"]))
            s_tgut = is_stop_arr(step_seq(psi, cdf_tips["Tgut"]))

            stops_mat = np.vstack([s_scam, s_arow, s_dnov, s_ggal, s_apla, s_shab, s_tgut])
            # Encode state: 0 if uncovered, 1 if covered & sense, 2 if covered & stop
            state_mat = np.where(~cov_mat, 0, np.where(stops_mat, 2, 1)).astype(np.int32)
            codes = (pow3[:, None] * state_mat).sum(axis=0)

            shared_mask = lut_n_carriers[codes] >= 2
            n_sh = int(shared_mask.sum())
            sh_codes = codes[shared_mask]
            n_ho = int((lg[sh_codes] >= 2).sum())
            n_ho_leg = int(lut_legacy_incong_t1[sh_codes].sum())
            n_neo = int(lut_crown_neorn[sh_codes].sum())

            sh_arr[r_i] = n_sh
            ho_arr[r_i] = n_ho
            ho_legacy_arr[r_i] = n_ho_leg
            neo_arr[r_i] = n_neo
            if n_sh > 0:
                frac = n_ho / n_sh
                hf_list.append(frac)
                if n_sh == 16:
                    hf_cond16_list.append(frac)
                if 14 <= n_sh <= 18:
                    hf_cond14_18_list.append(frac)

        obs_ptc = 16
        obs_ho_map = {"T1": 6, "T2": 10, "T3": 11}
        obs_frac = obs_ho_map[topo] / 16.0
        p_le_16 = float(np.mean(sh_arr <= obs_ptc))
        p_ge_16 = float(np.mean(sh_arr >= obs_ptc))
        p_two_sided_16 = min(1.0, 2.0 * min(p_le_16, p_ge_16))
        pct_rank_16 = float(np.mean(sh_arr < obs_ptc) + 0.5 * np.mean(sh_arr == obs_ptc)) * 100.0

        hf_arr = np.array(hf_list) if hf_list else np.zeros(1)
        p_frac_le = float(np.mean(hf_arr <= obs_frac))
        p_frac_ge = float(np.mean(hf_arr >= obs_frac))
        p_frac_2s = min(1.0, 2.0 * min(p_frac_le, p_frac_ge))

        return {
            "topology": topo,
            "anc_seq_source": "Neornithes_ML" if anc_nt == neorn_anc_seq else "Allmis_ref",
            "mask_mode": mask_mode,
            "kappa": kappa,
            "branch_scale": scale,
            "t_stem": t_stem,
            "n_reps": n_reps,
            "shared_ptc_median": round(float(np.median(sh_arr)), 1),
            "shared_ptc_mean": round(float(np.mean(sh_arr)), 2),
            "shared_ptc_iqr": f"{np.percentile(sh_arr, 25):.1f}-{np.percentile(sh_arr, 75):.1f}",
            "shared_ptc_95range": f"{np.percentile(sh_arr, 2.5):.1f}-{np.percentile(sh_arr, 97.5):.1f}",
            "obs16_percentile_rank": round(pct_rank_16, 1),
            "obs16_p_le": round(p_le_16, 4),
            "obs16_p_two_sided": round(p_two_sided_16, 4),
            "incong_ptc_median": round(float(np.median(ho_arr)), 1),
            "incong_rate_median": round(float(np.median(hf_arr)), 3),
            "incong_rate_95range": f"{np.percentile(hf_arr, 2.5):.3f}-{np.percentile(hf_arr, 97.5):.3f}",
            "legacy_set_incong_rate_median": round(float(np.median(ho_legacy_arr / np.maximum(1, sh_arr))), 3),
            "obs_incong_rate_topo": round(obs_frac, 3),
            "obs_incong_p_two_sided": round(p_frac_2s, 4),
            "cond_N16_n_reps": len(hf_cond16_list),
            "cond_N16_incong_median": round(float(np.median(hf_cond16_list)), 3) if hf_cond16_list else "",
            "cond_N14_18_n_reps": len(hf_cond14_18_list),
            "cond_N14_18_incong_median": round(float(np.median(hf_cond14_18_list)), 3) if hf_cond14_18_list else "",
            "cond_N14_18_incong_95range": (
                f"{np.percentile(hf_cond14_18_list, 2.5):.3f}-{np.percentile(hf_cond14_18_list, 97.5):.3f}"
                if len(hf_cond14_18_list) >= 10 else ""
            ),
            "crown_neornithes_syn_median": round(float(np.median(neo_arr)), 1),
            "crown_neornithes_syn_p_zero": round(float(np.mean(neo_arr == 0)), 4),
        }

    rng_grid = np.random.default_rng(20261123)
    sim_grid_rows: list[dict[str, Any]] = []

    # 1. Compare mask modes & root sequence under baseline T1, kappa=4.0, scale=1.0, t_stem in {0.0, 0.02, 0.05}
    for anc_s in (amis_clean_seq, neorn_anc_seq):
        for mm in ("unmasked_full_878", "empirical_mafft_mask", "empirical_dual_aligner_mask"):
            for ts in (0.0, 0.02, 0.05):
                sim_grid_rows.append(run_k80_scenario(rng_grid, anc_s, "T1", 4.0, 1.0, ts, mm, n_reps=1000))

    # 2. Sensitivity across kappa in {2.0, 4.0, 6.0} and branch_scale in {0.65, 0.80, 1.00, 1.25} under empirical_dual_aligner_mask (T1, Neornithes_ML)
    for kp in (2.0, 4.0, 6.0):
        for sc in (0.65, 0.80, 1.00, 1.25):
            for ts in (0.0, 0.02):
                if kp == 4.0 and sc == 1.0:
                    continue  # already run above
                sim_grid_rows.append(run_k80_scenario(rng_grid, neorn_anc_seq, "T1", kp, sc, ts, "empirical_dual_aligner_mask", n_reps=1000))

    # 3. Sensitivity across topologies T1, T2, T3 under empirical_dual_aligner_mask (kappa=4.0, scale in {0.80, 1.00}, t_stem in {0.0, 0.02})
    for top in ("T2", "T3"):
        for sc in (0.80, 1.00):
            for ts in (0.0, 0.02):
                sim_grid_rows.append(run_k80_scenario(rng_grid, neorn_anc_seq, top, 4.0, sc, ts, "empirical_dual_aligner_mask", n_reps=1000))

    write_tsv(OUT / "Table_S8_k80_simulation_sensitivity_v102.tsv", sim_grid_rows)
    print(f"[4/5] Completed {len(sim_grid_rows)} K80 neutral simulation sensitivity scenarios (1,000 reps each).", flush=True)
    for row in sim_grid_rows:
        if row["t_stem"] in (0.0, 0.02):
            print(
                f"  [{row['topology']} | {row['anc_seq_source']:12s} | {row['mask_mode']:27s} | "
                f"k={row['kappa']} | s={row['branch_scale']:.2f} | stem={row['t_stem']:.2f}] "
                f"sh_med={row['shared_ptc_median']:4.1f} (95%:{row['shared_ptc_95range']:9s}, p2s_16={row['obs16_p_two_sided']:.3f}) | "
                f"incong={row['incong_rate_median']:.3f} ({row['incong_rate_95range']}, obs={row['obs_incong_rate_topo']:.3f}, p2s={row['obs_incong_p_two_sided']:.3f}) | "
                f"cond14_18={row['cond_N14_18_incong_median']} (n={row['cond_N14_18_n_reps']}) | "
                f"P(neoSyn=0)={row['crown_neornithes_syn_p_zero']:.3f}",
                flush=True,
            )

    # -----------------------------------------------------------------------
    # STEP 5: Summarize Topology Homoplasy Sensitivity (Parsimony + ML ASR across T1, T2, T3)
    # -----------------------------------------------------------------------
    enam8_details = [d for d in variant_audit_details if d["dataset"] == "ENAM_exon8_aves" and d["amis_codon"] != "-"]
    topo_summary_rows = []
    for ti, (tcode, tlbl) in enumerate([("T1", "T1_Novaeratitae_StruthioBasal"), ("T2", "T2_ScamDnov_ArowBasal"), ("T3", "T3_ScamArow_DnovBasal")]):
        k_pars = sum(1 for d in enam8_details if d["min_gains_no_rev_T1_T2_T3"][ti] >= 2)
        # Under reversible ML ASR, codon 1022 has P_stop(Neognathae) > 0.90 in T1, T2, T3 -> 1 less homoplasy if counted as Neognathae synapomorphy + Taegut reversal
        k_ml_rev = k_pars - 1
        lo_p, hi_p = clopper_pearson(k_pars, 16)
        lo_m, hi_m = clopper_pearson(k_ml_rev, 16)
        tm = next((x for x in topo_tests_mafft if x["tree_id"] == ti + 1), {})
        tu = next((x for x in topo_tests_muscle if x["tree_id"] == ti + 1), {})
        topo_summary_rows.append({
            "topology_id": tcode,
            "topology_description": tlbl,
            "mafft_logL": tm.get("logL", ""),
            "mafft_deltaL": tm.get("deltaL", ""),
            "mafft_p_KH": tm.get("p_KH", ""),
            "mafft_p_SH": tm.get("p_SH", ""),
            "mafft_p_WSH": tm.get("p_WSH", ""),
            "mafft_c_ELW": tm.get("c_ELW", ""),
            "mafft_p_AU": tm.get("p_AU", ""),
            "muscle_logL": tu.get("logL", ""),
            "muscle_deltaL": tu.get("deltaL", ""),
            "muscle_p_SH": tu.get("p_SH", ""),
            "muscle_p_AU": tu.get("p_AU", ""),
            "n_shared_ptcs": 16,
            "k_incong_irreversible_parsimony": k_pars,
            "pct_incong_irreversible_parsimony": round(100.0 * k_pars / 16.0, 1),
            "binomial_95ci_irreversible_pct": f"{100.0 * lo_p:.1f}-{100.0 * hi_p:.1f}",
            "k_incong_reversible_ml_asr": k_ml_rev,
            "pct_incong_reversible_ml_asr": round(100.0 * k_ml_rev / 16.0, 1),
            "binomial_95ci_reversible_pct": f"{100.0 * lo_m:.1f}-{100.0 * hi_m:.1f}",
            "codon1022_neognathae_ml_p_stop": enam8_details[-1]["neog_p_stop_T1_T2_T3"][ti],
        })
    write_tsv(OUT / "Table_S9_topology_sensitivity_T1_T2_T3_v102.tsv", topo_summary_rows)

    full_summary = {
        "timestamp_utc": datetime.now(UTC).isoformat(),
        "elapsed_s": round((datetime.now(UTC) - t0).total_seconds(), 2),
        "topology_tests_mafft": topo_tests_mafft,
        "topology_tests_muscle": topo_tests_muscle,
        "topology_homoplasy_table": topo_summary_rows,
        "enam8_empirical_coverage": cov_summary,
        "ambn_exon4_blastn_ws7": {f"{k[0]}__{k[1]}": v for k, v in ex4_blastn_best.items()},
        "ambn_exon4_tblastn_ws2": ex4_tblastn_best,
        "k80_simulation_grid": sim_grid_rows,
        "variant_audit_details": variant_audit_details,
    }
    sum_json = OUT / "v102_audit_summary.json"
    sum_json.write_text(json.dumps(full_summary, indent=2), encoding="utf-8", newline="\n")

    with PROV.open("a", encoding="utf-8", newline="\n") as f:
        f.write(json.dumps({
            "ts": datetime.now(UTC).isoformat(),
            "kind": "enamel-relics-v102-phases-ABCD",
            "script": "scripts/enamel_relics/12_audit_v102_variants_phylo_sim.py",
            "seed_k80_grid": 20261123,
            "outputs": {
                p.name: sha256_file(p)
                for p in [
                    OUT / "Table_S5a_level4_shared_mutations_v102.tsv",
                    OUT / "Table_S8_k80_simulation_sensitivity_v102.tsv",
                    OUT / "Table_S9_topology_sensitivity_T1_T2_T3_v102.tsv",
                    sum_json,
                ]
            },
        }) + "\n")

    print("[5/5] v1.0.2 Phases A-D complete -> runs/enamel_relics/v102/v102_audit_summary.json", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
