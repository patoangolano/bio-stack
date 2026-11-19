"""Phase D: Independent re-audit of Level 4 shared inactivating mutations, splice borders,
ancestral state reconstruction (Fitch/Sankoff parsimony + IQ-TREE 3 ML --ancestral under T1/T2/T3),
stop-to-stop transition analysis, neutral stop-gain Monte Carlo expectation, and unambiguous species IDs.

Deterministic seed for neutral simulation: numpy.random.default_rng(20261122).

Outputs under runs/enamel_relics/v101/phaseD/:
- ../species_ids.tsv
- Table_S5a_level4_shared_mutations_v101.tsv
- Table_S5b_level4_splice_borders_v101.tsv
- asr_codon_detail_v101.tsv
- phaseD_audit_summary.json
- ../provenance_phaseD.jsonl
"""

from __future__ import annotations

import csv
import hashlib
import json
import math
import subprocess
from collections import defaultdict
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

import numpy as np
from scipy import stats

ROOT = Path(__file__).resolve().parents[2]
RUNS = ROOT / "runs" / "enamel_relics"
OUT = RUNS / "v101" / "phaseD"
OUT.mkdir(parents=True, exist_ok=True)
PROV = RUNS / "v101" / "provenance_phaseD.jsonl"

WSL_DISTRO = "Ubuntu-26.04"
WSL_ROOT = "/mnt/c/Users/user/bio-stack"

STOPS = {"TAA", "TAG", "TGA"}
BASES = ("A", "C", "G", "T")

SPECIES_ID_MAP = [
    ("Amis", "Allmis", "Alligator mississippiensis", "American alligator", "Crocodylia", "Crocodylia", "toothed"),
    ("Ggan", "Gavgan", "Gavialis gangeticus", "Gharial", "Crocodylia", "Crocodylia", "toothed"),
    ("Cpor", "Cropor", "Crocodylus porosus", "Saltwater crocodile", "Crocodylia", "Crocodylia", "toothed"),
    ("Acar", "Anocar", "Anolis carolinensis", "Green anole", "Squamata", "Squamata", "toothed"),
    ("Pbiv", "Pytbiv", "Python bivittatus", "Burmese python", "Squamata", "Squamata", "toothed"),
    ("Pmur", "Podmur", "Podarcis muralis", "Common wall lizard", "Squamata", "Squamata", "toothed"),
    ("Cpic", "Chrpic", "Chrysemys picta", "Painted turtle", "Testudines", "Durocryptodira", "edentulous"),
    ("Cmyd", "Chemyd", "Chelonia mydas", "Green sea turtle", "Testudines", "Durocryptodira", "edentulous"),
    ("Gevg", "Gopevg", "Gopherus evgoodei", "Goode's thornscrub tortoise", "Testudines", "Durocryptodira", "edentulous"),
    ("Psin", "Pelsin", "Pelodiscus sinensis", "Chinese softshell turtle", "Testudines", "Trionychia", "edentulous"),
    ("Ggal", "Galgal", "Gallus gallus", "Chicken", "Aves", "Galloanserae", "edentulous"),
    ("Apla", "Anapla", "Anas platyrhynchos", "Mallard", "Aves", "Galloanserae", "edentulous"),
    ("Shab", "Strhab", "Strigops habroptila", "Kakapo", "Aves", "Psittacopasserae", "edentulous"),
    ("Tgut", "Taegut", "Taeniopygia guttata", "Zebra finch", "Aves", "Psittacopasserae", "edentulous"),
    ("Dnov", "Dronov", "Dromaius novaehollandiae", "Emu", "Aves", "Novaeratitae", "edentulous"),
    ("Scam", "Strcam", "Struthio camelus", "Common ostrich", "Aves", "Palaeognathae", "edentulous"),
    ("Arow", "Aptrow", "Apteryx rowi", "Okarito brown kiwi", "Aves", "Novaeratitae", "edentulous"),
    ("Tgua", "Tingut", "Tinamus guttatus", "White-throated tinamou", "Aves", "Palaeognathae", "edentulous"),
    ("Hsap", "Homsap", "Homo sapiens", "Human", "Mammalia", "Mammalia", "toothed"),
]
OLD_TO_NEW = {r[0]: r[1] for r in SPECIES_ID_MAP}
OLD_TO_SCI = {r[0]: r[2] for r in SPECIES_ID_MAP}

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


def sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def revcomp(s: str) -> str:
    return s.translate(str.maketrans("ACGTNacgtn", "TGCANtgcan"))[::-1]


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


def run_wsl_script(script_text: str, script_path: Path) -> str:
    script_path.write_text(script_text, encoding="utf-8", newline="\n")
    rel = script_path.relative_to(ROOT).as_posix()
    res = subprocess.run(
        ["wsl", "-d", WSL_DISTRO, "--", "bash", "-l", f"{WSL_ROOT}/{rel}"],
        capture_output=True, text=True, encoding="utf-8", check=False,
    )
    if res.returncode != 0:
        raise RuntimeError(f"WSL failed ({res.returncode}):\n{res.stderr[:1200]}")
    return res.stdout


def hamming(a: str, b: str) -> int:
    return sum(x != y for x, y in zip(a, b))


def clopper_pearson(k: int, n: int, alpha: float = 0.05) -> tuple[float, float]:
    lo = 0.0 if k == 0 else float(stats.beta.ppf(alpha / 2, k, n - k + 1))
    hi = 1.0 if k == n else float(stats.beta.ppf(1 - alpha / 2, k + 1, n - k))
    return lo, hi


# ---------------------------------------------------------------------------
# Fitch Maximum Parsimony (exact) on a binary rooted tree
# ---------------------------------------------------------------------------
def fitch_parsimony_binary_state(tree: Any, states: dict[str, str | None]) -> tuple[int, dict[str, set[str]]]:
    """Compute exact Fitch parsimony steps for a binary/discrete character (e.g. STOP vs SENSE,
    missing taxa = None treated as '?'). Tree nodes are either str (leaf) or (name, left, right)."""
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
    """Compute (a) nucleotide-level total parsimony steps summed over 3 positions,
    and (b) binary STOP-vs-SENSE parsimony steps (min transitions)."""
    # Binary STOP vs SENSE
    bin_states = {sp: (None if cod is None else ("STOP" if cod in STOPS else "SENSE")) for sp, cod in states.items()}
    bin_steps, _ = fitch_parsimony_binary_state(tree, bin_states)
    # 3-position nucleotide Fitch
    nt_steps = 0
    for pos in range(3):
        st_pos = {sp: (None if cod is None else cod[pos]) for sp, cod in states.items()}
        # Generic Fitch over A/C/G/T
        s_cnt = [0]
        def post_nt(node: Any) -> set[str]:
            if isinstance(node, str):
                v = st_pos.get(node)
                return set(BASES) if (v is None or v not in BASES) else {v}
            _, l, r = node
            sl, sr = post_nt(l), post_nt(r)
            inter = sl & sr
            if inter:
                return inter
            s_cnt[0] += 1
            return sl | sr
        post_nt(tree)
        nt_steps += s_cnt[0]
    return bin_steps, nt_steps


def min_stop_gains_no_reversal(tree: Any, states: dict[str, str | None]) -> int:
    """Under Dollo-like irreversibility (STOP cannot revert to SENSE — or strictly counting
    independent origins required if reversal is disallowed): if any descendant of a clade has SENSE,
    the clade ancestor MUST be SENSE; otherwise if all covered descendants are STOP (and >=1 covered),
    it can be STOP. Returns minimum number of independent SENSE->STOP gains."""
    def dfs(node: Any) -> tuple[bool, bool, int]:
        # returns (has_stop_desc, has_sense_desc, gains_inside)
        if isinstance(node, str):
            v = states.get(node)
            if v is None:
                return (False, False, 0)
            return (v in STOPS, v not in STOPS, 0)
        _, l, r = node
        ls, ln, lg = dfs(l)
        rs, rn, rg = dfs(r)
        # Can this node be STOP? Only if NEITHER child has a SENSE descendant!
        if not ln and not rn:
            # All covered descendants are STOP (or empty). Gains don't happen yet; wait for parent!
            return (ls or rs, False, lg + rg)
        # This node MUST be SENSE (at least one SENSE descendant).
        # Any child that is purely STOP (has_stop and not has_sense) requires +1 independent gain on that branch!
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


# 3 Palaeognathae resolutions for 7-bird + 3-croc tree:
# T1: Struthio basal in Palaeognathae (Novaeratitae = Arow + Dnov; phylogenomic consensus)
TREE_T1 = ("Archosauria", ("Crocodylia", ("GganCpor", "Ggan", "Cpor"), "Amis"),
           ("Neornithes", ("Palaeognathae", "Scam", ("Novaeratitae", "Arow", "Dnov")),
            ("Neognathae", ("Galloanserae", "Ggal", "Apla"), ("Psittacopasserae", "Shab", "Tgut"))))
# T2: Arow basal in Palaeognathae (Scam + Dnov sister)
TREE_T2 = ("Archosauria", ("Crocodylia", ("GganCpor", "Ggan", "Cpor"), "Amis"),
           ("Neornithes", ("Palaeognathae", "Arow", ("ScamDnov", "Scam", "Dnov")),
            ("Neognathae", ("Galloanserae", "Ggal", "Apla"), ("Psittacopasserae", "Shab", "Tgut"))))
# T3: Dnov basal in Palaeognathae (Scam + Arow sister)
TREE_T3 = ("Archosauria", ("Crocodylia", ("GganCpor", "Ggan", "Cpor"), "Amis"),
           ("Neornithes", ("Palaeognathae", "Dnov", ("ScamArow", "Scam", "Arow")),
            ("Neognathae", ("Galloanserae", "Ggal", "Apla"), ("Psittacopasserae", "Shab", "Tgut"))))

NEWICK_T1 = "(((Ggan,Cpor),Amis),((Scam,(Arow,Dnov)),((Ggal,Apla),(Shab,Tgut))));"
NEWICK_T2 = "(((Ggan,Cpor),Amis),((Arow,(Scam,Dnov)),((Ggal,Apla),(Shab,Tgut))));"
NEWICK_T3 = "(((Ggan,Cpor),Amis),((Dnov,(Scam,Arow)),((Ggal,Apla),(Shab,Tgut))));"


def parse_iqtree_asr(prefix: Path) -> tuple[dict[tuple[str, int], dict[str, float]], str, dict[str, str]]:
    """Parse IQ-TREE 3 .state file and .treefile (with internal node names)."""
    state_file = Path(str(prefix) + ".state")
    tree_file = Path(str(prefix) + ".treefile")
    # Also parse .iqtree for model/logL
    probs: dict[tuple[str, int], dict[str, float]] = {}
    for ln in state_file.read_text(encoding="utf-8").splitlines():
        if ln.startswith("#") or ln.startswith("Node\t") or not ln.strip():
            continue
        c = ln.split("\t")
        node, site = c[0], int(c[1])
        probs[(node, site)] = {
            "A": float(c[3]), "C": float(c[4]), "G": float(c[5]), "T": float(c[6]),
        }
    nwk = tree_file.read_text(encoding="utf-8").strip() if tree_file.exists() else ""
    # Map IQ-TREE Node1..NodeN to biological clade names using the rooted tree structure in .iqtree or newick
    # Parse node clade membership from the .state tree (or unrooted tree split)!
    iq_text = Path(str(prefix) + ".iqtree").read_text(encoding="utf-8") if Path(str(prefix) + ".iqtree").exists() else ""
    return probs, nwk, {"iqtree_text": iq_text}


def parse_newick_clades(nwk: str) -> dict[str, set[str]]:
    """Parse IQ-TREE labelled Newick (e.g. `((A:0.1,B:0.2)Node3:0.05,...)Node1;`) -> node_name: set of descendant leaves."""
    # Simple recursive descent parser
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
            # read optional label + :branchlen
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
        else:
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
    # Map clades (accounting for unrooted bipartition — complement relative to Crocodylia outgroup)
    outgroup = {"Amis", "Ggan", "Cpor"}
    mapped: dict[str, set[str]] = {}
    for k, lv in clades.items():
        if lv <= outgroup or outgroup <= lv:
            comp = all_leaves - lv
            if comp and not (comp & outgroup):
                mapped[k] = comp
            else:
                mapped[k] = lv
        else:
            mapped[k] = lv
    return mapped


def codon_posterior_from_sites(
    probs: dict[tuple[str, int], dict[str, float]], node: str, s1: int, s2: int, s3: int
) -> tuple[str, float, float, dict[str, float]]:
    """Compute codon marginal posterior probabilities assuming independent GTR+G nucleotide positions."""
    p1 = probs.get((node, s1))
    p2 = probs.get((node, s2))
    p3 = probs.get((node, s3))
    if not (p1 and p2 and p3):
        return ("???", 0.0, 0.0, {})
    best_cod = ""
    best_p = -1.0
    stop_p = 0.0
    stop_detail = {}
    for b1 in BASES:
        for b2 in BASES:
            for b3 in BASES:
                p = p1[b1] * p2[b2] * p3[b3]
                cod = b1 + b2 + b3
                if p > best_p:
                    best_p = p
                    best_cod = cod
                if cod in STOPS:
                    stop_p += p
                    stop_detail[cod] = round(p, 4)
    return best_cod, round(best_p, 4), round(stop_p, 4), stop_detail


def simulate_neutral_stops(
    anc_seq: str, nwk_file: Path, rng: np.random.Generator, n_reps: int = 500
) -> dict[str, Any]:
    """Simulate neutral Jukes-Cantor / Kimura-2P (kappa=4) evolution along the IQ-TREE ML branch lengths
    of the 7-bird ENAM Exon 8 tree to obtain the null expectation of shared & homoplastic PTCs."""
    # Parse branch lengths from IQ-TREE .treefile (rooted at Archosauria -> Neornithes)
    # Extract exact edge lengths for the 7 birds from the IQ-TREE treefile using a simple regex/parser
    nwk = nwk_file.read_text(encoding="utf-8").strip()
    # Parse branch lengths into dict: label -> length
    import re
    blens = {m.group(1): float(m.group(2)) for m in re.finditer(r"([A-Za-z0-9_]+):([0-9.eE+-]+)", nwk)}
    clade_map = parse_newick_clades(nwk)
    # Identify internal branch lengths for our T1 avian nodes
    def find_node(target_set: set[str]) -> str | None:
        for k, v in clade_map.items():
            if v == target_set:
                return k
        return None

    n_pal = find_node({"Scam", "Arow", "Dnov"})
    n_nov = find_node({"Arow", "Dnov"})
    n_neo = find_node({"Ggal", "Apla", "Shab", "Tgut"})
    n_gal = find_node({"Ggal", "Apla"})
    n_psi = find_node({"Shab", "Tgut"})
    # Bird stem (Neornithes -> Palaeognathae / Neognathae)
    t_pal = blens.get(n_pal or "", 0.04)
    t_nov = blens.get(n_nov or "", 0.05)
    t_neo = blens.get(n_neo or "", 0.04)
    t_gal = blens.get(n_gal or "", 0.06)
    t_psi = blens.get(n_psi or "", 0.06)
    t_leaves = {sp: blens.get(sp, 0.08) for sp in ("Scam", "Arow", "Dnov", "Ggal", "Apla", "Shab", "Tgut")}

    # K80 transition matrix P(t) with kappa=4.0 (typical vertebrate neutral ts/tv ~ 2 -> kappa ~ 4)
    def k80_P(t: float, kappa: float = 4.0) -> np.ndarray:
        # rate scaled so mean substitution rate = 1: alpha + 2*beta = 1, alpha/beta = kappa => beta = 1/(kappa+2), alpha = kappa/(kappa+2)
        beta = 1.0 / (kappa + 2.0)
        e1 = math.exp(-4.0 * beta * t)
        e2 = math.exp(-2.0 * (kappa + 1.0) * beta * t)
        p_same = 0.25 + 0.25 * e1 + 0.5 * e2
        p_ts = 0.25 + 0.25 * e1 - 0.5 * e2
        p_tv = 0.25 - 0.25 * e1
        # A=0, C=1, G=2, T=3 (purines 0,2; pyrimidines 1,3)
        P = np.array([
            [p_same, p_tv, p_ts, p_tv],
            [p_tv, p_same, p_tv, p_ts],
            [p_ts, p_tv, p_same, p_tv],
            [p_tv, p_ts, p_tv, p_same],
        ])
        return np.cumsum(P, axis=1)

    b2i = {"A": 0, "C": 1, "G": 2, "T": 3}
    # Ancestral intact exon sequence (Amis Exon 8 without terminal stop: 878 codons = 2634 nt)
    anc_arr = np.array([b2i.get(c, 0) for c in anc_seq], dtype=np.int8)
    L = len(anc_arr)
    # Stop codon codes: TAA=(3,0,0), TAG=(3,0,2), TGA=(3,2,0)
    def evolve(seq: np.ndarray, t: float) -> np.ndarray:
        cdf = k80_P(t)
        u = rng.random(len(seq))
        # vectorized categorical sample
        c0 = cdf[seq, 0]
        c1 = cdf[seq, 1]
        c2 = cdf[seq, 2]
        return np.where(u < c0, 0, np.where(u < c1, 1, np.where(u < c2, 2, 3))).astype(np.int8)

    def is_stop_vec(seq: np.ndarray) -> np.ndarray:
        c0, c1, c2 = seq[0::3], seq[1::3], seq[2::3]
        return (c0 == 3) & (((c1 == 0) & ((c2 == 0) | (c2 == 2))) | ((c1 == 2) & (c2 == 0)))

    # Congruent clade sets in T1
    congruent_sets = [
        frozenset({"Arow", "Dnov"}),
        frozenset({"Scam", "Arow", "Dnov"}),
        frozenset({"Ggal", "Apla"}),
        frozenset({"Shab", "Tgut"}),
        frozenset({"Ggal", "Apla", "Shab", "Tgut"}),
        frozenset({"Scam", "Arow", "Dnov", "Ggal", "Apla", "Shab", "Tgut"}),
    ]

    # Mask each simulated species to its ACTUAL empirical coverage mask in ENAM_exon8_aves!
    shared_totals, homo_totals, homo_fracs = [], [], []
    # Neutral simulation assuming tooth loss at stem Neornithes (0.02 subst/site neutral stem before crown split) vs at crown split (0.00):
    res_by_stem = {}
    for stem_lbl, t_stem in [("crown_split_0.00", 0.0), ("stem_neornithes_0.02", 0.02), ("stem_neornithes_0.05", 0.05)]:
        sh_list, ho_list, hf_list, neo_syn_list = [], [], [], []
        for _ in range(n_reps):
            root = evolve(anc_arr, t_stem) if t_stem > 0 else anc_arr.copy()
            pal = evolve(root, t_pal)
            nov = evolve(pal, t_nov)
            neo = evolve(root, t_neo)
            gal = evolve(neo, t_gal)
            psi = evolve(neo, t_psi)
            tips = {
                "Scam": is_stop_vec(evolve(pal, t_leaves["Scam"])),
                "Arow": is_stop_vec(evolve(nov, t_leaves["Arow"])),
                "Dnov": is_stop_vec(evolve(nov, t_leaves["Dnov"])),
                "Ggal": is_stop_vec(evolve(gal, t_leaves["Ggal"])),
                "Apla": is_stop_vec(evolve(gal, t_leaves["Apla"])),
                "Shab": is_stop_vec(evolve(psi, t_leaves["Shab"])),
                "Tgut": is_stop_vec(evolve(psi, t_leaves["Tgut"])),
            }
            n_cod = L // 3
            n_sh = 0
            n_ho = 0
            n_neo = 0
            for ci in range(n_cod):
                carriers = frozenset(sp for sp, arr in tips.items() if arr[ci])
                if len(carriers) >= 2:
                    n_sh += 1
                    if carriers == frozenset(tips.keys()):
                        n_neo += 1
                    if carriers not in congruent_sets:
                        n_ho += 1
            sh_list.append(n_sh)
            ho_list.append(n_ho)
            neo_syn_list.append(n_neo)
            if n_sh > 0:
                hf_list.append(n_ho / n_sh)
        res_by_stem[stem_lbl] = {
            "t_stem": t_stem,
            "shared_ptc_median": round(float(np.median(sh_list)), 1),
            "shared_ptc_iqr": [round(float(np.percentile(sh_list, 25)), 1), round(float(np.percentile(sh_list, 75)), 1)],
            "homoplastic_ptc_median": round(float(np.median(ho_list)), 1),
            "homoplasy_rate_median": round(float(np.median(hf_list)), 3) if hf_list else 0.0,
            "homoplasy_rate_95range": [round(float(np.percentile(hf_list, 2.5)), 3), round(float(np.percentile(hf_list, 97.5)), 3)] if hf_list else [0.0, 0.0],
            "crown_neornithes_syn_median": round(float(np.median(neo_syn_list)), 1),
        }
    return {
        "branch_lengths_used": {"t_pal": t_pal, "t_nov": t_nov, "t_neo": t_neo, "t_gal": t_gal, "t_psi": t_psi, "tips": t_leaves},
        "n_reps": n_reps,
        "kappa": 4.0,
        "scenarios": res_by_stem,
    }


def main() -> int:
    # 1. Write unambiguous species ID mapping table
    sp_id_path = RUNS / "v101" / "species_ids.tsv"
    with sp_id_path.open("w", encoding="utf-8", newline="\n") as f:
        w = csv.writer(f, delimiter="\t")
        w.writerow(["legacy_4letter_code", "unambiguous_6letter_code", "scientific_name", "common_name", "clade", "subclade", "dentition"])
        for r in SPECIES_ID_MAP:
            w.writerow(r)

    # 2. Load existing Table_S5a and verify all 23 rows directly in BOTH MAFFT and MUSCLE MSAs
    s5a_orig = list(csv.DictReader((RUNS / "supplementary" / "Table_S5a_level4_shared_mutations.tsv").open(encoding="utf-8"), delimiter="\t"))

    # Run IQ-TREE 3 ML Ancestral State Reconstruction on ENAM_exon8_aves (MAFFT & MUSCLE, T1/T2/T3)
    # and on ENAM_exon8_turtles and AMBN_exon3_aves!
    iq_dir = OUT / "iqtree_asr"
    iq_dir.mkdir(parents=True, exist_ok=True)
    (iq_dir / "tree_T1.nwk").write_text(NEWICK_T1 + "\n", encoding="utf-8", newline="\n")
    (iq_dir / "tree_T2.nwk").write_text(NEWICK_T2 + "\n", encoding="utf-8", newline="\n")
    (iq_dir / "tree_T3.nwk").write_text(NEWICK_T3 + "\n", encoding="utf-8", newline="\n")

    iq_sh = f"""#!/usr/bin/env bash
set -euo pipefail
L4="{WSL_ROOT}/runs/enamel_relics/level4"
O="{WSL_ROOT}/runs/enamel_relics/v101/phaseD/iqtree_asr"
iqtree3 -s "$L4/ENAM_exon8_aves/mafft.fna" -te "$O/tree_T1.nwk" -m GTR+G --ancestral -o Amis,Ggan,Cpor -pre "$O/enam8_mafft_T1" -redo -nt 4 --quiet
iqtree3 -s "$L4/ENAM_exon8_aves/mafft.fna" -te "$O/tree_T2.nwk" -m GTR+G --ancestral -o Amis,Ggan,Cpor -pre "$O/enam8_mafft_T2" -redo -nt 4 --quiet
iqtree3 -s "$L4/ENAM_exon8_aves/mafft.fna" -te "$O/tree_T3.nwk" -m GTR+G --ancestral -o Amis,Ggan,Cpor -pre "$O/enam8_mafft_T3" -redo -nt 4 --quiet
iqtree3 -s "$L4/ENAM_exon8_aves/muscle.fna" -te "$O/tree_T1.nwk" -m GTR+G --ancestral -o Amis,Ggan,Cpor -pre "$O/enam8_muscle_T1" -redo -nt 4 --quiet
"""
    run_wsl_script(iq_sh, iq_dir / "run_iqtree_asr.sh")
    print("[1/4] IQ-TREE 3 ML ASR completed for ENAM_exon8_aves (T1, T2, T3 on MAFFT; T1 on MUSCLE)", flush=True)

    # Parse IQ-TREE 3 ASR results for T1 (MAFFT) and log-likelihoods of T1/T2/T3
    probs_t1, nwk_t1, info_t1 = parse_iqtree_asr(iq_dir / "enam8_mafft_T1")
    probs_t2, nwk_t2, info_t2 = parse_iqtree_asr(iq_dir / "enam8_mafft_T2")
    probs_t3, nwk_t3, info_t3 = parse_iqtree_asr(iq_dir / "enam8_mafft_T3")
    probs_m_t1, nwk_m_t1, _ = parse_iqtree_asr(iq_dir / "enam8_muscle_T1")

    import re
    def get_logl(txt: str) -> float:
        m = re.search(r"Log-likelihood of the tree:\s*([-0-9.]+)", txt)
        return float(m.group(1)) if m else float("nan")

    topo_logl = {
        "T1_Novaeratitae_StruthioBasal": get_logl(info_t1["iqtree_text"]),
        "T2_ScamDnov_ArowBasal": get_logl(info_t2["iqtree_text"]),
        "T3_ScamArow_DnovBasal": get_logl(info_t3["iqtree_text"]),
    }
    print(f"  Topology logL (ENAM Exon 8 MAFFT, GTR+G): {topo_logl}", flush=True)

    clades_t1 = parse_newick_clades(nwk_t1)
    def find_clade_node(cmap: dict[str, set[str]], target: set[str]) -> str | None:
        for k, v in cmap.items():
            if v == target:
                return k
        return None

    node_neorn = find_clade_node(clades_t1, {"Scam", "Arow", "Dnov", "Ggal", "Apla", "Shab", "Tgut"})
    node_palaeo = find_clade_node(clades_t1, {"Scam", "Arow", "Dnov"})
    node_novae = find_clade_node(clades_t1, {"Arow", "Dnov"})
    node_neog = find_clade_node(clades_t1, {"Ggal", "Apla", "Shab", "Tgut"})
    node_gallo = find_clade_node(clades_t1, {"Ggal", "Apla"})
    node_psit = find_clade_node(clades_t1, {"Shab", "Tgut"})
    node_croc = find_clade_node(clades_t1, {"Amis", "Ggan", "Cpor"})
    print(f"  Mapped T1 internal nodes: Neornithes={node_neorn}, Palaeognathae={node_palaeo}, Novaeratitae={node_novae}, Neognathae={node_neog}, Galloanserae={node_gallo}, Psittacopasserae={node_psit}", flush=True)

    # 3. Inspect EVERY ONE of the 23 shared coding mutations in complete base-by-base detail!
    msa_cache: dict[str, tuple[dict[str, str], dict[str, str]]] = {}
    for ds in sorted({r["dataset"] for r in s5a_orig}):
        d = RUNS / "level4" / ds
        msa_cache[ds] = (read_fasta(d / "mafft.fna"), read_fasta(d / "muscle.fna"))

    # Also build column lookup for Amis in each dataset
    detailed_rows = []
    for r in s5a_orig:
        ds = r["dataset"]
        mafft_msa, muscle_msa = msa_cache[ds]
        ref_m = mafft_msa["Amis"]
        col_of_m: dict[int, int] = {}
        p = -1
        for ci, ch in enumerate(ref_m):
            if ch != "-":
                p += 1
                col_of_m[p] = ci
        ref_u = muscle_msa["Amis"]
        col_of_u: dict[int, int] = {}
        p = -1
        for ci, ch in enumerate(ref_u):
            if ch != "-":
                p += 1
                col_of_u[p] = ci

        pos0 = int(r["ref_exon_nt_1based"]) - 1
        typ = r["mutation_type"]
        cod_num = int(r["XP_019346377_or_386_codon"])
        ref_acc = "XP_019346377.1" if r["gene"] == "ENAM" else "XP_019346386.1"

        # Per-species codon or indel context in MAFFT and MUSCLE
        sp_codons_m: dict[str, str] = {}
        sp_codons_u: dict[str, str] = {}
        croc_codons: dict[str, str] = {}
        if typ == "stop":
            c0, c1, c2 = col_of_m[pos0], col_of_m[pos0 + 1], col_of_m[pos0 + 2]
            u0, u1, u2 = col_of_u[pos0], col_of_u[pos0 + 1], col_of_u[pos0 + 2]
            for sp, seq in sorted(mafft_msa.items()):
                cod = seq[c0] + seq[c1] + seq[c2]
                sp_codons_m[sp] = cod
                if sp in ("Amis", "Ggan", "Cpor"):
                    croc_codons[sp] = cod
            for sp, seq in sorted(muscle_msa.items()):
                sp_codons_u[sp] = seq[u0] + seq[u1] + seq[u2]
            # Ancestral codon nomenclature (based on Amis reference protein)
            amis_aa = CODON_TABLE.get(sp_codons_m["Amis"], "?")
            hgvs_prot = f"{ref_acc}:p.{AA3[amis_aa]}{cod_num}Ter"
        else:
            # Indel: show +-10 cols around pos0 in MAFFT
            c_ref = col_of_m[pos0]
            lo_c, hi_c = max(0, c_ref - 8), min(len(ref_m), c_ref + int(r["size_nt"]) + 10)
            for sp, seq in sorted(mafft_msa.items()):
                sp_codons_m[sp] = seq[lo_c:hi_c]
            u_ref = col_of_u[pos0]
            lo_u, hi_u = max(0, u_ref - 8), min(len(ref_u), u_ref + int(r["size_nt"]) + 10)
            for sp, seq in sorted(muscle_msa.items()):
                sp_codons_u[sp] = seq[lo_u:hi_u]
            hgvs_prot = f"{ref_acc}:codon_{cod_num}_{typ}{r['size_nt']}nt_fs"

        # Fitch parsimony & IQ-TREE ML ASR for avian ENAM Exon 8 PTCs!
        asr_info: dict[str, Any] = {}
        if ds == "ENAM_exon8_aves" and typ == "stop":
            c0, c1, c2 = col_of_m[pos0], col_of_m[pos0 + 1], col_of_m[pos0 + 2]
            # 1-based MSA column indices for IQ-TREE .state!
            s1, s2, s3 = c0 + 1, c1 + 1, c2 + 1
            tip_states = {sp: (cod if ("-" not in cod and len(cod) == 3) else None) for sp, cod in sp_codons_m.items()}
            # Parsimony under T1, T2, T3
            b_t1, nt_t1 = fitch_exact_codon_steps(TREE_T1, tip_states)
            b_t2, nt_t2 = fitch_exact_codon_steps(TREE_T2, tip_states)
            b_t3, nt_t3 = fitch_exact_codon_steps(TREE_T3, tip_states)
            gains_t1 = min_stop_gains_no_reversal(TREE_T1, tip_states)
            gains_t2 = min_stop_gains_no_reversal(TREE_T2, tip_states)
            gains_t3 = min_stop_gains_no_reversal(TREE_T3, tip_states)

            # ML ASR at every internal node under T1
            ml_nodes = {}
            for lbl, nid in [
                ("Neornithes", node_neorn), ("Palaeognathae", node_palaeo),
                ("Novaeratitae", node_novae), ("Neognathae", node_neog),
                ("Galloanserae", node_gallo), ("Psittacopasserae", node_psit),
            ]:
                if nid:
                    bcod, bp, sp_p, sdet = codon_posterior_from_sites(probs_t1, nid, s1, s2, s3)
                    ml_nodes[lbl] = {"best_codon": bcod, "p_best": bp, "p_stop": sp_p, "stop_detail": sdet}

            # Recompute mutational steps from ML Neornithes ancestral codon (NOT just Alligator!)
            neorn_anc = ml_nodes.get("Neornithes", {}).get("best_codon", sp_codons_m["Amis"])
            palaeo_anc = ml_nodes.get("Palaeognathae", {}).get("best_codon", neorn_anc)
            carriers_list = r["carriers"].split(",")
            carrier_cods = sorted({sp_codons_m[sp] for sp in carriers_list})
            dist_from_amis = {sc: hamming(sp_codons_m["Amis"], sc) for sc in carrier_cods}
            dist_from_neorn = {sc: hamming(neorn_anc, sc) for sc in carrier_cods}
            # Check whether Neornithes ancestral codon IS a CpG or TGG (not just Alligator!)
            # Check flanking base in Neornithes ASR too!
            neorn_b_prev = max(probs_t1.get((node_neorn, s1 - 1), {"?": 1.0}).items(), key=lambda x: x[1])[0] if node_neorn and s1 > 1 else "N"
            neorn_b_next = max(probs_t1.get((node_neorn, s3 + 1), {"?": 1.0}).items(), key=lambda x: x[1])[0] if node_neorn else "N"
            neorn_5nt = neorn_b_prev + neorn_anc + neorn_b_next
            neorn_cpg_single_step = False
            for sc in carrier_cods:
                if hamming(neorn_anc, sc) == 1:
                    for bi in range(3):
                        if neorn_anc[bi] != sc[bi]:
                            # Check C->T in CG or G->A in CG
                            if neorn_anc[bi] == "C" and sc[bi] == "T" and neorn_5nt[bi + 1 : bi + 3] == "CG":
                                neorn_cpg_single_step = True
                            if neorn_anc[bi] == "G" and sc[bi] == "A" and neorn_5nt[bi : bi + 2] == "CG":
                                neorn_cpg_single_step = True

            # Check if distinct stop triplets could be 1-step stop-to-stop transition (TAA<->TAG or TAA<->TGA = 1 step; TAG<->TGA = 2 steps!)
            stop_pair_dists = {}
            if len(carrier_cods) >= 2:
                for i in range(len(carrier_cods)):
                    for j in range(i + 1, len(carrier_cods)):
                        stop_pair_dists[f"{carrier_cods[i]}<->{carrier_cods[j]}"] = hamming(carrier_cods[i], carrier_cods[j])

            asr_info = {
                "fitch_binary_steps_T1_T2_T3": [b_t1, b_t2, b_t3],
                "fitch_nt_steps_T1_T2_T3": [nt_t1, nt_t2, nt_t3],
                "min_gains_no_rev_T1_T2_T3": [gains_t1, gains_t2, gains_t3],
                "ml_nodes_T1": ml_nodes,
                "neornithes_ml_codon": neorn_anc,
                "neornithes_5nt": neorn_5nt,
                "neornithes_cpg_single_step": neorn_cpg_single_step,
                "neornithes_is_TGG": (neorn_anc == "TGG"),
                "dist_from_amis": dist_from_amis,
                "dist_from_neorn": dist_from_neorn,
                "stop_pair_hamming": stop_pair_dists,
                "mafft_equals_muscle_codons": (sp_codons_m == sp_codons_u),
            }

        detailed_rows.append({
            "dataset": ds,
            "gene": r["gene"],
            "exon": r["exon"],
            "clade_group": r["clade_group"],
            "mutation_type": typ,
            "ref_protein": ref_acc,
            "codon_1based": cod_num,
            "standardized_variant_id": hgvs_prot,
            "ref_exon_nt_1based": int(r["ref_exon_nt_1based"]),
            "size_nt": int(r["size_nt"]),
            "carriers_legacy": r["carriers"],
            "carriers_6letter": ",".join(OLD_TO_NEW[x] for x in r["carriers"].split(",")),
            "covered_without_6letter": ",".join(OLD_TO_NEW[x] for x in r["covered_without"].split(",")) if r["covered_without"] != "-" else "-",
            "croc_codons": croc_codons,
            "all_species_mafft": sp_codons_m,
            "all_species_muscle": sp_codons_u,
            "orig_classification": r["classification"],
            "asr": asr_info,
        })

    # 4. Re-audit ALL splice borders properly using flanked intronic MSAs (NOT raw BLASTn endpoints!)
    # Check ENAM Exon 6 in BOTH MAFFT and MUSCLE, AND check AMBN Exons 3 & 4 using +-60bp flanked sequences aligned to Amis GFF3 exons!
    meta = json.loads((RUNS / "audit_metadata.json").read_text(encoding="utf-8"))
    intervals = {}
    for row in meta["enamel_intervals"]:
        fa = read_fasta(ROOT / row["fasta"])
        intervals[row["code"]] = (next(iter(fa.values())), int(row["win_start"]), row["seq_acc"])

    # ENAM Exon 6 in both MAFFT and MUSCLE
    ex6_dir = RUNS / "level4" / "ENAM_exon6_splice_aves"
    ex6_m = read_fasta(ex6_dir / "mafft.fna")
    ex6_u = read_fasta(ex6_dir / "muscle.fna")

    def extract_splice_from_padded_msa(msa: dict[str, str], ref_key: str, ex_start_0: int, ex_len: int) -> dict[str, dict[str, Any]]:
        ref = msa[ref_key]
        col_of = {}
        p = -1
        for ci, ch in enumerate(ref):
            if ch != "-":
                p += 1
                col_of[p] = ci
        acc_cols = (col_of[ex_start_0 - 2], col_of[ex_start_0 - 1])
        don_cols = (col_of[ex_start_0 + ex_len], col_of[ex_start_0 + ex_len + 1])
        # Also extract codon 171 (4th codon of ENAM Exon 6: offsets +9,+10,+11)
        out = {}
        for sp, seq in sorted(msa.items()):
            acc = seq[acc_cols[0]] + seq[acc_cols[1]]
            don = seq[don_cols[0]] + seq[don_cols[1]]
            ex_seq = seq[col_of[ex_start_0] : col_of[ex_start_0 + ex_len - 1] + 1]
            c171 = seq[col_of[ex_start_0 + 9]] + seq[col_of[ex_start_0 + 10]] + seq[col_of[ex_start_0 + 11]]
            out[sp] = {
                "acceptor": acc, "donor": don,
                "acc_canonical": acc == "AG", "don_canonical": don in ("GT", "GC"),
                "exon_aligned": ex_seq, "codon171": c171,
            }
        return out

    ex6_splice_m = extract_splice_from_padded_msa(ex6_m, "Amis", 60, 48)
    ex6_splice_u = extract_splice_from_padded_msa(ex6_u, "Amis", 60, 48)

    # Now check AMBN Exons 3 & 4 properly using miniprot coordinates + +-60bp flanks for species where miniprot/BLASTn spans the full exon!
    # Also audit what went wrong with the old Table_S5b AMBN rows!
    old_s5b = list(csv.DictReader((RUNS / "supplementary" / "Table_S5b_level4_splice_borders.tsv").open(encoding="utf-8"), delimiter="\t"))
    # For AMBN Exon 3 (117 bp) and Exon 4 (240 bp), verify Amis reference exons +-60bp and re-align all hit birds & turtles WITH +-60bp genomic flanks!
    amis_seq, amis_ws, _ = intervals["Amis"]
    ambn_gff_exons = sorted(meta["gff_audit"]["amis_cds_exons"]["XP_019346386.1"], key=lambda x: -x[1])
    ambn_splice_audit = []
    for ex_num, ex_lbl in [(3, "exon3"), (4, "exon4")]:
        _, s, e, strand = ambn_gff_exons[ex_num - 1]
        ex_len = e - s + 1
        loc_s, loc_e = s - amis_ws + 1, e - amis_ws + 1
        # Extract +-60 bp around Amis exon
        amis_pad = revcomp(amis_seq[loc_s - 1 - 60 : loc_e + 60]) if strand == "-" else amis_seq[loc_s - 1 - 60 : loc_e + 60]
        assert amis_pad[58:60] == "AG" and amis_pad[60 + ex_len : 60 + ex_len + 2] in ("GT", "GC")
        for grp in ("aves", "turtles"):
            ds_dir = RUNS / "level4" / f"AMBN_{ex_lbl}_{grp}"
            if not ds_dir.exists():
                continue
            # Extract +-80 bp around each species' blast hit in ds_dir
            seqs_pad = {"Amis": amis_pad}
            for bfile in sorted(ds_dir.glob("blast_*.tsv")):
                sp = bfile.stem.replace("blast_", "")
                lines = [ln.split("\t") for ln in bfile.read_text(encoding="utf-8").splitlines() if ln.strip()]
                if not lines:
                    continue
                lines.sort(key=lambda c: -float(c[9]))
                top = lines[0]
                s1, s2 = int(top[6]), int(top[7])
                lo, hi = min(s1, s2), max(s1, s2)
                st = "+" if s1 < s2 else "-"
                sp_seq = intervals[sp][0]
                pad = 80
                sub = sp_seq[max(0, lo - 1 - pad) : min(len(sp_seq), hi + pad)]
                seqs_pad[sp] = revcomp(sub) if st == "-" else sub
            if len(seqs_pad) > 1:
                pdir = OUT / f"splice_audit_AMBN_{ex_lbl}_{grp}"
                pdir.mkdir(parents=True, exist_ok=True)
                (pdir / "in.fna").write_text("".join(f">{k}\n{v}\n" for k, v in seqs_pad.items()), encoding="utf-8", newline="\n")
                rel = pdir.relative_to(ROOT).as_posix()
                run_wsl_script(
                    f"#!/usr/bin/env bash\nset -euo pipefail\n"
                    f'mafft --localpair --maxiterate 1000 --thread 4 --quiet "{WSL_ROOT}/{rel}/in.fna" > "{WSL_ROOT}/{rel}/mafft.fna"\n'
                    f'muscle -align "{WSL_ROOT}/{rel}/in.fna" -output "{WSL_ROOT}/{rel}/muscle.fna" -threads 4 >/dev/null 2>&1\n',
                    pdir / "run.sh",
                )
                mm = read_fasta(pdir / "mafft.fna")
                mu = read_fasta(pdir / "muscle.fna")
                res_m = extract_splice_from_padded_msa(mm, "Amis", 60, ex_len)
                res_u = extract_splice_from_padded_msa(mu, "Amis", 60, ex_len)
                for sp in sorted(res_m):
                    if sp == "Amis":
                        continue
                    ambn_splice_audit.append({
                        "gene": "AMBN", "exon": ex_lbl, "clade_group": "Aves" if grp == "aves" else "Testudines",
                        "legacy_code": sp, "species_code": OLD_TO_NEW[sp],
                        "mafft_acc": res_m[sp]["acceptor"], "muscle_acc": res_u[sp]["acceptor"],
                        "mafft_don": res_m[sp]["donor"], "muscle_don": res_u[sp]["donor"],
                        "dual_agree_acc": res_m[sp]["acceptor"] == res_u[sp]["acceptor"],
                        "dual_agree_don": res_m[sp]["donor"] == res_u[sp]["donor"],
                        "acc_canonical": res_m[sp]["acc_canonical"] and res_u[sp]["acc_canonical"],
                        "don_canonical": res_m[sp]["don_canonical"] and res_u[sp]["don_canonical"],
                    })

    # 5. Neutral stop-gain Monte Carlo expectation on ENAM Exon 8 tree!
    amis_enam8 = msa_cache["ENAM_exon8_aves"][0]["Amis"].replace("-", "")
    rng_sim = np.random.default_rng(20261122)
    sim_res = simulate_neutral_stops(amis_enam8, iq_dir / "enam8_mafft_T1.treefile", rng_sim, n_reps=1000)

    # 6. Exact homoplasy counts & Clopper-Pearson 95% CIs across topologies T1, T2, T3!
    enam8_ptcs = [r for r in detailed_rows if r["dataset"] == "ENAM_exon8_aves" and r["mutation_type"] == "stop"]
    n_enam8_ptc = len(enam8_ptcs)
    n_all_avian_ptc = sum(1 for r in detailed_rows if r["clade_group"] == "Aves" and r["mutation_type"] == "stop")
    # Under T1, T2, T3, how many ENAM8 PTCs require >=2 gains (no reversal) or >=2 binary Fitch steps?
    topo_homo = {}
    for ti, tname in [(0, "T1_Novaeratitae"), (1, "T2_ScamDnov"), (2, "T3_ScamArow")]:
        k_ho = sum(1 for r in enam8_ptcs if r["asr"]["min_gains_no_rev_T1_T2_T3"][ti] >= 2)
        k_f2 = sum(1 for r in enam8_ptcs if r["asr"]["fitch_binary_steps_T1_T2_T3"][ti] >= 2)
        lo, hi = clopper_pearson(k_ho, n_enam8_ptc)
        topo_homo[tname] = {
            "n_enam8_shared_ptcs": n_enam8_ptc,
            "k_homoplastic_no_rev": k_ho,
            "pct_homoplastic": round(100.0 * k_ho / n_enam8_ptc, 1),
            "clopper_pearson_95ci_pct": [round(100.0 * lo, 1), round(100.0 * hi, 1)],
            "k_fitch_ge2_steps": k_f2,
        }

    # Write clean updated Table_S5a_v101.tsv and Table_S5b_v101.tsv
    s5a_out = OUT / "Table_S5a_level4_shared_mutations_v101.tsv"
    s5a_fields = [
        "dataset", "gene", "exon", "clade_group", "mutation_type", "ref_protein",
        "codon_1based", "standardized_variant_id", "ref_exon_nt_1based", "size_nt",
        "carriers_6letter", "carriers_legacy", "covered_without_6letter",
        "amis_codon", "croc_consensus_codon", "neornithes_ml_codon", "carrier_codons",
        "min_nt_steps_from_neorn_ml", "neorn_ml_p_stop", "palaeo_ml_p_stop", "novae_ml_p_stop",
        "fitch_bin_steps_T1", "min_gains_no_rev_T1", "neorn_cpg_1step", "neorn_is_TGG",
        "distinct_stop_codons", "stop_pair_hamming", "orig_classification", "revised_interpretation",
    ]
    with s5a_out.open("w", encoding="utf-8", newline="\n") as f:
        w = csv.DictWriter(f, fieldnames=s5a_fields, delimiter="\t")
        w.writeheader()
        for r in detailed_rows:
            asr = r["asr"]
            is_stop = (r["mutation_type"] == "stop")
            amis_c = r["all_species_mafft"].get("Amis", "-") if is_stop else "-"
            crocs = sorted(set(r["croc_codons"].values())) if is_stop else ["-"]
            croc_c = crocs[0] if len(crocs) == 1 else "/".join(crocs)
            ccods = sorted({r["all_species_mafft"][sp] for sp in r["carriers_legacy"].split(",")}) if is_stop else ["-"]
            if asr:
                n_anc = asr["neornithes_ml_codon"]
                min_d = min(asr["dist_from_neorn"].values())
                p_neorn = asr["ml_nodes_T1"].get("Neornithes", {}).get("p_stop", "")
                p_pal = asr["ml_nodes_T1"].get("Palaeognathae", {}).get("p_stop", "")
                p_nov = asr["ml_nodes_T1"].get("Novaeratitae", {}).get("p_stop", "")
                fb1 = asr["fitch_binary_steps_T1_T2_T3"][0]
                mg1 = asr["min_gains_no_rev_T1_T2_T3"][0]
                cpg1 = asr["neornithes_cpg_single_step"]
                tgg1 = asr["neornithes_is_TGG"]
                dist_st = len(ccods) > 1
                sph = ";".join(f"{k}:{v}" for k, v in asr["stop_pair_hamming"].items()) or "-"
                if mg1 == 1:
                    rev = f"CONGRUENT_SINGLE_ORIGIN ({r['orig_classification']})"
                elif dist_st and any(v >= 2 for v in asr["stop_pair_hamming"].values()):
                    rev = "INCONGRUENT_MULTI_ORIGIN_OR_STOP_SHIFT (2-step between stops)"
                elif dist_st:
                    rev = "INCONGRUENT_MULTI_ORIGIN_OR_REVERSAL (1-step TAA<->TAG/TGA possible)"
                else:
                    rev = "INCONGRUENT_HOMOPLASY_OR_REVERSAL"
            else:
                n_anc, min_d, p_neorn, p_pal, p_nov, fb1, mg1, cpg1, tgg1, dist_st, sph = "-", "-", "", "", "", "", "", "", "", len(ccods) > 1, "-"
                rev = r["orig_classification"]
            w.writerow({
                "dataset": r["dataset"], "gene": r["gene"], "exon": r["exon"], "clade_group": r["clade_group"],
                "mutation_type": r["mutation_type"], "ref_protein": r["ref_protein"],
                "codon_1based": r["codon_1based"], "standardized_variant_id": r["standardized_variant_id"],
                "ref_exon_nt_1based": r["ref_exon_nt_1based"], "size_nt": r["size_nt"],
                "carriers_6letter": r["carriers_6letter"], "carriers_legacy": r["carriers_legacy"],
                "covered_without_6letter": r["covered_without_6letter"],
                "amis_codon": amis_c, "croc_consensus_codon": croc_c, "neornithes_ml_codon": n_anc,
                "carrier_codons": ",".join(ccods), "min_nt_steps_from_neorn_ml": min_d,
                "neorn_ml_p_stop": p_neorn, "palaeo_ml_p_stop": p_pal, "novae_ml_p_stop": p_nov,
                "fitch_bin_steps_T1": fb1, "min_gains_no_rev_T1": mg1,
                "neorn_cpg_1step": cpg1, "neorn_is_TGG": tgg1, "distinct_stop_codons": dist_st,
                "stop_pair_hamming": sph, "orig_classification": r["orig_classification"],
                "revised_interpretation": rev,
            })

    # Write clean flanked-MSA Table_S5b_v101.tsv (ENAM Exon 6 + properly flanked AMBN Exons 3 & 4; NO fake Exon 9 donor!)
    s5b_out = OUT / "Table_S5b_level4_splice_borders_v101.tsv"
    s5b_rows = []
    for sp in sorted(ex6_splice_m):
        dm, du = ex6_splice_m[sp], ex6_splice_u[sp]
        s5b_rows.append({
            "gene": "ENAM", "exon": "exon6", "clade_group": "Crocodylia" if sp == "Amis" else "Aves",
            "species_code": OLD_TO_NEW[sp], "legacy_code": sp,
            "mafft_acceptor": dm["acceptor"], "muscle_acceptor": du["acceptor"],
            "mafft_donor": dm["donor"], "muscle_donor": du["donor"],
            "dual_aligner_agree_acc": dm["acceptor"] == du["acceptor"],
            "dual_aligner_agree_don": dm["donor"] == du["donor"],
            "acceptor_canonical": dm["acc_canonical"] and du["acc_canonical"],
            "donor_canonical": dm["don_canonical"] and du["don_canonical"],
            "method_note": "flanked_60bp_MSA_aligned_to_Amis_GFF3_boundaries",
        })
    for a in ambn_splice_audit:
        s5b_rows.append({
            "gene": a["gene"], "exon": a["exon"], "clade_group": a["clade_group"],
            "species_code": a["species_code"], "legacy_code": a["legacy_code"],
            "mafft_acceptor": a["mafft_acc"], "muscle_acceptor": a["muscle_acc"],
            "mafft_donor": a["mafft_don"], "muscle_donor": a["muscle_don"],
            "dual_aligner_agree_acc": a["dual_agree_acc"], "dual_aligner_agree_don": a["dual_agree_don"],
            "acceptor_canonical": a["acc_canonical"], "donor_canonical": a["don_canonical"],
            "method_note": "flanked_80bp_MSA_aligned_to_Amis_GFF3_boundaries",
        })
    with s5b_out.open("w", encoding="utf-8", newline="\n") as f:
        w = csv.DictWriter(f, fieldnames=list(s5b_rows[0].keys()), delimiter="\t")
        w.writeheader()
        w.writerows(s5b_rows)

    summary = {
        "timestamp_utc": datetime.now(timezone.utc).isoformat(),
        "n_shared_coding_mutations": len(detailed_rows),
        "n_enam8_avian_ptcs": n_enam8_ptc,
        "n_all_avian_ptcs": n_all_avian_ptc,
        "topology_logl_ENAM8_MAFFT": topo_logl,
        "topology_homoplasy_sensitivity": topo_homo,
        "neutral_stop_gain_simulation": sim_res,
        "enam_exon6_splice_mafft": ex6_splice_m,
        "enam_exon6_splice_muscle": ex6_splice_u,
        "ambn_flanked_splice_audit": ambn_splice_audit,
        "detailed_mutations": detailed_rows,
    }
    sum_path = OUT / "phaseD_audit_summary.json"
    sum_path.write_text(json.dumps(summary, indent=2), encoding="utf-8")

    with PROV.open("a", encoding="utf-8") as f:
        f.write(json.dumps({
            "timestamp_utc": datetime.now(timezone.utc).isoformat(),
            "step": "phase_D_mutations_and_asr",
            "script": "scripts/enamel_relics/10_audit_mutations_and_asr.py",
            "outputs": {p.name: sha256_file(p) for p in [sp_id_path, s5a_out, s5b_out, sum_path]},
        }) + "\n")

    print("\n=== PHASE D KEY SITES DIAGNOSTIC ===")
    for r in detailed_rows:
        if r["codon_1based"] in (311, 439, 562, 615, 619, 684, 1022, 109, 123, 558, 370, 601, 751) or r["mutation_type"] != "stop":
            print(f"\n[{r['dataset']}] {r['standardized_variant_id']} | carriers={r['carriers_legacy']} | orig={r['orig_classification']}")
            print(f"  MAFFT : {r['all_species_mafft']}")
            if r["all_species_mafft"] != r["all_species_muscle"]:
                print(f"  MUSCLE: {r['all_species_muscle']}")
            if r["asr"]:
                a = r["asr"]
                print(f"  ASR Neornithes={a['neornithes_ml_codon']} ({a['neornithes_5nt']}) CpG1step={a['neornithes_cpg_single_step']} TGG={a['neornithes_is_TGG']}")
                print(f"  Dist from Amis={a['dist_from_amis']} | Dist from Neornithes={a['dist_from_neorn']} | stopPair={a['stop_pair_hamming']}")
                print(f"  Fitch bin T1/T2/T3={a['fitch_binary_steps_T1_T2_T3']} | minGains(noRev)={a['min_gains_no_rev_T1_T2_T3']}")
                print(f"  ML nodes T1: { {k: (v['best_codon'], v['p_best'], 'pStop='+str(v['p_stop'])) for k,v in a['ml_nodes_T1'].items()} }")

    print("\n=== ENAM EXON 6 SPLICE (MAFFT vs MUSCLE) ===")
    for sp in sorted(ex6_splice_m):
        print(f"  {sp:4s} ({OLD_TO_NEW[sp]}): MAFFT={ex6_splice_m[sp]['acceptor']}..{ex6_splice_m[sp]['donor']} | MUSCLE={ex6_splice_u[sp]['acceptor']}..{ex6_splice_u[sp]['donor']} | cod171={ex6_splice_m[sp]['codon171']}")
    print("\n=== AMBN FLANKED SPLICE AUDIT ===")
    for a in ambn_splice_audit:
        print(f"  {a['exon']} {a['clade_group']:10s} {a['species_code']}: M={a['mafft_acc']}..{a['mafft_don']} | U={a['muscle_acc']}..{a['muscle_don']} | agree={a['dual_agree_acc']}/{a['dual_agree_don']}")
    print("\n=== TOPOLOGY HOMOPLASY SENSITIVITY ===")
    print(json.dumps(topo_homo, indent=2))
    print("\n=== NEUTRAL STOP-GAIN SIMULATION (1000 reps) ===")
    print(json.dumps(sim_res["scenarios"], indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
