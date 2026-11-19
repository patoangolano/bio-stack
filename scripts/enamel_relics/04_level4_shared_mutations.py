"""Phase 2.4 — Exon-aware Level 4 Shared Inactivating Mutation Analysis
(dual-aligner concordance: MAFFT L-INS-i + MUSCLE 5, homopolymer filter >=4 nt,
CpG/TGG context annotation, canonical GT-AG splice-site audit, and phylogenetic
congruence vs homoplasy classification) across:
  1. Avian ENAM Exon 8 (full terminal coding exon, codons 201-1078 of XP_019346377.1)
  2. Avian ENAM N-terminal Exons 6-7 + Intron 6 GT->AT splice-site check (from Arow)
  3. Testudine ENAM Exon 8 (4 turtles + 3 crocodylians)
  4. Avian AMBN major coding exons (Exons 1, 3, 4, 5, 6, 9 + splice borders)
  5. Testudine AMBN coding exons (Exons 1, 3, 4, 5, 6, 9 + Intron GT->AT splice audit)

Also audits Cpic ENAM's 58-kb miniprot intron jump to confirm whether the N-terminal
block is bona fide Exon 4 or a spurious distal alignment.

Usage:
    uv run python scripts/enamel_relics/04_level4_shared_mutations.py
"""

from __future__ import annotations

import csv
import hashlib
import json
import re
import subprocess
from collections import defaultdict
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / "runs" / "enamel_relics"
L4 = OUT / "level4"
TABLES = OUT / "tables"

CROCS = ["Amis", "Ggan", "Cpor"]
AVES = ["Ggal", "Apla", "Shab", "Tgut", "Dnov", "Scam", "Arow", "Tgua"]
TURTLES = ["Cpic", "Cmyd", "Gevg", "Psin"]
STOPS = {"TAA", "TAG", "TGA"}

# Accepted topologies for strict MRCA congruence checks
# Aves:
#   Palaeognathae = {Scam, Tgua, Arow, Dnov}; Notopalaeognathae = {Tgua, Arow, Dnov}; Novaeratitae = {Arow, Dnov}
#   Neognathae = {Ggal, Apla, Shab, Tgut}; Galloanserae = {Ggal, Apla}; Psittacopasserae = {Shab, Tgut}
AVIAN_CLADES: list[tuple[str, set[str]]] = [
    ("Neornithes", set(AVES)),
    ("Palaeognathae", {"Scam", "Tgua", "Arow", "Dnov"}),
    ("Notopalaeognathae", {"Tgua", "Arow", "Dnov"}),
    ("Novaeratitae", {"Arow", "Dnov"}),
    ("Neognathae", {"Ggal", "Apla", "Shab", "Tgut"}),
    ("Galloanserae", {"Ggal", "Apla"}),
    ("Psittacopasserae", {"Shab", "Tgut"}),
]

# Testudines (Cryptodira):
#   Cryptodira = {Psin, Cmyd, Cpic, Gevg}
#   Durocryptodira = {Cmyd, Cpic, Gevg}
#   Testudinoidea = {Cpic, Gevg}
TURTLE_CLADES: list[tuple[str, set[str]]] = [
    ("Cryptodira", set(TURTLES)),
    ("Durocryptodira", {"Cmyd", "Cpic", "Gevg"}),
    ("Testudinoidea", {"Cpic", "Gevg"}),
]

COMP = str.maketrans("ACGTNacgtn", "TGCANtgcan")
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


def revcomp(s: str) -> str:
    return s.translate(COMP)[::-1]


def translate(dna: str) -> str:
    dna = dna.upper()
    return "".join(CODON_TABLE.get(dna[i : i + 3], "X") for i in range(0, len(dna) - 2, 3))


def read_fasta(path: Path) -> dict[str, str]:
    seqs: dict[str, list[str]] = {}
    cur = ""
    for ln in path.read_text(encoding="utf-8").splitlines():
        if ln.startswith(">"):
            cur = ln[1:].split()[0]
            seqs[cur] = []
        elif cur:
            seqs[cur].append(ln.strip().upper())
    return {k: "".join(v) for k, v in seqs.items()}


def sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def load_intervals(meta: dict[str, Any]) -> dict[str, tuple[str, int, str]]:
    """Return sp -> (seq_string, win_start_1based, seq_acc)."""
    out: dict[str, tuple[str, int, str]] = {}
    for row in meta["enamel_intervals"]:
        seqs = read_fasta(ROOT / row["fasta"])
        out[row["code"]] = (next(iter(seqs.values())), int(row["win_start"]), row["seq_acc"])
    return out


def extract_subseq(seq: str, start_1based: int, end_1based: int, strand: str) -> str:
    lo = max(1, min(start_1based, end_1based))
    hi = min(len(seq), max(start_1based, end_1based))
    sub = seq[lo - 1 : hi]
    return revcomp(sub) if strand == "-" else sub


def parse_miniprot_exons(paf_path: Path) -> list[dict[str, Any]]:
    """Parse best miniprot alignment into exact exon blocks (0-based target coordinates on local interval)
    and query aa ranges + flanking intron lengths."""
    if not paf_path.exists():
        return []
    lines = [l for l in paf_path.read_text(encoding="utf-8").splitlines() if l.strip()]
    if not lines:
        return []
    best = max(lines, key=lambda l: int(next(t[5:] for t in l.split("\t")[12:] if t.startswith("AS:i:"))))
    c = best.split("\t")
    qini, qfim = int(c[2]), int(c[3])
    strand = c[4]
    tini, tfim = int(c[7]), int(c[8])
    cg = next(t[5:] for t in c[12:] if t.startswith("cg:Z:"))

    # Walk CIGAR in coding direction (on '-' strand, target position descends from tfim down to tini!)
    exons: list[dict[str, Any]] = []
    cur_q_start = qini
    cur_t_start = tini if strand == "+" else tfim
    qpos = qini
    tpos = cur_t_start

    def close_exon(t_end: int, q_end: int, split_phase: int = 0) -> None:
        nonlocal cur_q_start, cur_t_start
        lo = min(cur_t_start, t_end)
        hi = max(cur_t_start, t_end)
        if hi > lo:
            exons.append({
                "strand": strand,
                "q_start_aa": cur_q_start,
                "q_end_aa": q_end,
                "t_lo_0": lo,
                "t_hi_0": hi,
                "bp": hi - lo,
                "split_phase_after": split_phase,
            })

    for n_str, op in re.findall(r"(\d+)([A-Z])", cg):
        n = int(n_str)
        sgn = 1 if strand == "+" else -1
        if op == "M":
            qpos += n
            tpos += sgn * 3 * n
        elif op == "D":
            tpos += sgn * 3 * n
        elif op == "I":
            qpos += n
        elif op == "F":
            # n nt insertion in target relative to query frame
            tpos += sgn * n
        elif op == "G":
            # n nt present in target for a deleted/frameshifted codon (consumes 1 query aa and n target nt)
            qpos += 1
            tpos += sgn * n
        elif op == "N":
            close_exon(tpos, qpos, split_phase=0)
            tpos += sgn * n
            cur_t_start = tpos
            cur_q_start = qpos
        elif op == "U":
            # 1 nt before intron, (n-3) nt intron, 2 nt after intron (total target advance = n, query +1)
            close_exon(tpos + sgn * 1, qpos + 1, split_phase=1)
            tpos += sgn * (n - 2)
            cur_t_start = tpos
            tpos += sgn * 2
            qpos += 1
            cur_q_start = qpos - 1
        elif op == "V":
            # 2 nt before intron, (n-3) nt intron, 1 nt after intron
            close_exon(tpos + sgn * 2, qpos + 1, split_phase=2)
            tpos += sgn * (n - 1)
            cur_t_start = tpos
            tpos += sgn * 1
            qpos += 1
            cur_q_start = qpos - 1
    close_exon(tpos, qpos, split_phase=0)
    return exons


def run_wsl_cmd(script_text: str, script_path: Path) -> str:
    script_path.write_text(script_text, encoding="utf-8", newline="\n")
    wsl_p = "/mnt/c/" + str(script_path.resolve())[3:].replace("\\", "/")
    res = subprocess.run(["wsl", "-d", "Ubuntu-26.04", "--", "bash", "-l", wsl_p], capture_output=True, text=True, check=False)
    if res.returncode != 0:
        raise RuntimeError(f"WSL script {script_path.name} failed:\nSTDERR:\n{res.stderr}\nSTDOUT:\n{res.stdout}")
    return res.stdout


def homopolymer_near(seq_aligned: str, i: int, j: int) -> bool:
    """True if >= 4 identical consecutive bases occur within +-4 columns around [i, j) (gaps removed)."""
    window = seq_aligned[max(0, i - 5) : min(len(seq_aligned), j + 5)].replace("-", "")
    return any(window[k : k + 4] == window[k] * 4 for k in range(max(0, len(window) - 3)))


def extract_msa_lesions(
    msa: dict[str, str],
    ref_key: str,
    croc_keys: list[str],
    codon_offset_1based: int,
    reading_frame_shift: int = 0,
) -> dict[str, set[tuple[str, int, int, bool]]]:
    """Identify deletions, insertions, and in-frame stops relative to Amis reference CDS exon.
    pos_ref is 0-based nt offset within the reference exon sequence."""
    ref = msa[ref_key]
    ncol = len(ref)
    pos_ref: list[int] = []
    p = -1
    for ch in ref:
        if ch != "-":
            p += 1
        pos_ref.append(p)
    croc_base = [sum(msa[cr][c] != "-" for cr in croc_keys if cr in msa) >= max(2, len(croc_keys) - 1) for c in range(ncol)]
    col_of: dict[int, int] = {pos_ref[c]: c for c in range(ncol) if ref[c] != "-"}
    ref_len = len(ref.replace("-", ""))

    out: dict[str, set[tuple[str, int, int, bool]]] = {}
    for sp, s in msa.items():
        if sp in croc_keys:
            continue
        cov = [i for i, b in enumerate(s) if b != "-"]
        if not cov:
            continue
        ini, fim = cov[0], cov[-1]
        ev: set[tuple[str, int, int, bool]] = set()
        c = ini
        while c <= fim:
            if s[c] == "-" and croc_base[c]:
                d = c
                while d <= fim and not (s[d] != "-" and croc_base[d]) and (s[d] == "-" or not croc_base[d]):
                    d += 1
                n = sum(1 for k in range(c, d) if croc_base[k] and s[k] == "-")
                if n > 0:
                    ev.add(("del", pos_ref[c], n, homopolymer_near(ref, c, d)))
                c = d
            elif s[c] != "-" and not croc_base[c]:
                d = c
                while d <= fim and not croc_base[d]:
                    d += 1
                n = sum(1 for k in range(c, d) if s[k] != "-")
                if n > 0:
                    ev.add(("ins", pos_ref[c], n, homopolymer_near(s, c, d)))
                c = d
            else:
                c += 1
        # In-frame stops on Amis coding frame
        for nt0 in range(reading_frame_shift, ref_len - 2, 3):
            cols = [col_of.get(nt0 + x) for x in range(3)]
            if any(x is None for x in cols):
                continue
            c0, c1, c2 = cols[0], cols[1], cols[2]
            assert c0 is not None and c1 is not None and c2 is not None
            if not (ini <= c0 and c2 <= fim):
                continue
            # Exclude final native stop codon of exon if Amis itself has a stop there
            ref_cod = ref[c0] + ref[c1] + ref[c2]
            if ref_cod in STOPS:
                continue
            cod = s[c0] + s[c1] + s[c2]
            if "-" not in cod and cod in STOPS:
                ev.add(("stop", nt0, 3, False))
        out[sp] = ev
    return out


def score_dual_aligner_dataset(
    dataset_id: str,
    gene: str,
    exon_label: str,
    clade_group: str,
    mafft_msa: dict[str, str],
    muscle_msa: dict[str, str],
    clade_defs: list[tuple[str, set[str]]],
    codon_offset_1based: int,
    reading_frame_shift: int = 0,
) -> tuple[list[dict[str, Any]], dict[str, Any]]:
    """Score shared lesions concordant between MAFFT L-INS-i and MUSCLE 5,
    plus crocodylian negative-control lesion count and full CpG/TGG context."""
    ev_m = extract_msa_lesions(mafft_msa, "Amis", CROCS, codon_offset_1based, reading_frame_shift)
    ev_u = extract_msa_lesions(muscle_msa, "Amis", CROCS, codon_offset_1based, reading_frame_shift)

    # Negative control: lesions in Ggan and Cpor vs Amis
    croc_ctrl = {"mafft": 0, "muscle": 0}
    for aname, msa in (("mafft", mafft_msa), ("muscle", muscle_msa)):
        for cr in ("Ggan", "Cpor"):
            if cr not in msa:
                continue
            sub_ev = extract_msa_lesions(msa, "Amis", ["Amis"], codon_offset_1based, reading_frame_shift)
            croc_ctrl[aname] += sum(1 for e in sub_ev.get(cr, set()) if (e[0] == "stop" or e[2] % 3 != 0) and not e[3])

    def shared_map(ev: dict[str, set[tuple[str, int, int, bool]]]) -> dict[tuple[str, int, int], set[str]]:
        m: dict[tuple[str, int, int], set[str]] = defaultdict(set)
        for sp, events in ev.items():
            for typ, pos, n, homo in events:
                if homo:
                    continue
                if typ == "stop" or (n % 3 != 0):
                    m[(typ, pos, n)].add(sp)
        return {k: v for k, v in m.items() if len(v) >= 2}

    s1, s2 = shared_map(ev_m), shared_map(ev_u)
    robust = {k: s1[k] & s2[k] for k in (s1.keys() & s2.keys()) if len(s1[k] & s2[k]) >= 2}

    # Reference indexing for coverage + ancestral codon + CpG dinucleotide context
    ref = mafft_msa["Amis"]
    ref_ungapped = ref.replace("-", "")
    col_of: dict[int, int] = {}
    p = -1
    for col, ch in enumerate(ref):
        if ch != "-":
            p += 1
            col_of[p] = col

    target_species = [sp for sp in mafft_msa if sp not in CROCS]

    def covers(sp: str, typ: str, pos: int, n: int) -> bool:
        s = mafft_msa[sp]
        if typ == "stop":
            # Species covers a stop codon only if all 3 codon positions are non-gap bases
            cols_cod = [col_of.get(pos + x) for x in range(3)]
            if any(c is None for c in cols_cod):
                return False
            return all(s[c] != "-" for c in cols_cod if c is not None)
        cols = [col_of[x] for x in range(max(0, pos - 3), min(len(ref_ungapped), pos + n + 3)) if x in col_of]
        return len(cols) >= 2 and s[cols[0]] != "-" and s[cols[-1]] != "-"

    rows: list[dict[str, Any]] = []
    for (typ, pos, n), carriers in sorted(robust.items(), key=lambda x: x[0][1]):
        # Ignore leading/trailing overhangs at position 0 or end of exon (intronic flanks)
        if typ in ("ins", "del") and (pos <= 2 or pos >= len(ref_ungapped) - 3):
            continue
        covered = {sp for sp in target_species if covers(sp, typ, pos, n)}
        covered_without = covered - carriers
        # Smallest clade containing all carriers
        valid_clades = [x for x in clade_defs if carriers <= x[1]]
        mrca_name, mrca_set = min(valid_clades, key=lambda x: len(x[1])) if valid_clades else ("Unresolved", set())
        congruent = bool(mrca_set) and ((covered & mrca_set) == carriers)

        codon_in_exon = (pos - reading_frame_shift) // 3
        full_protein_codon = codon_offset_1based + codon_in_exon
        anc_codon = ref_ungapped[pos : pos + 3] if typ == "stop" else "-"
        carrier_stops = sorted({
            "".join(mafft_msa[sp][col_of[pos + x]] for x in range(3))
            for sp in carriers
        }) if typ == "stop" else []

        # Full CpG context check (including cross-codon boundaries: pos-1..pos+3)
        ctx_5 = ref_ungapped[max(0, pos - 1) : min(len(ref_ungapped), pos + 4)] if typ == "stop" else "-"
        cpg_site = False
        if typ == "stop":
            # Check if any mutated base (differing from anc_codon in any carrier stop) is part of a CG in ref
            for st in carrier_stops:
                for b_i in range(3):
                    if b_i < len(anc_codon) and b_i < len(st) and st[b_i] != anc_codon[b_i]:
                        g_i = pos + b_i
                        if (g_i > 0 and ref_ungapped[g_i - 1 : g_i + 1] == "CG") or (
                            g_i + 1 < len(ref_ungapped) and ref_ungapped[g_i : g_i + 2] == "CG"
                        ):
                            cpg_site = True
        tgg_hotspot = (typ == "stop" and anc_codon == "TGG")
        multi_stop = (typ == "stop" and len(carrier_stops) > 1)

        if not congruent:
            classification = "INCONGRUENT_HOMOPLASY"
        elif multi_stop:
            classification = f"CONGRUENT_{mrca_name}_DISTINCT_STOPS"
        else:
            classification = f"CONGRUENT_{mrca_name}"

        rows.append({
            "dataset": dataset_id,
            "gene": gene,
            "exon": exon_label,
            "clade_group": clade_group,
            "mutation_type": typ,
            "ref_exon_nt_1based": pos + 1,
            "XP_019346377_or_386_codon": full_protein_codon,
            "size_nt": n,
            "carriers": ",".join(sorted(carriers)),
            "n_carriers": len(carriers),
            "covered_without": ",".join(sorted(covered_without)) or "-",
            "ancestral_croc_codon": anc_codon,
            "carrier_stop_codons": ",".join(carrier_stops) if carrier_stops else "-",
            "ref_5nt_context": ctx_5,
            "at_CpG_dinucleotide": cpg_site,
            "at_TGG_codon": tgg_hotspot,
            "distinct_stop_codons": multi_stop,
            "mrca_clade": mrca_name,
            "tree_congruent": congruent,
            "classification": classification,
        })

    # Per-species total inactivating lesions (concordant or per aligner) for diagnostics
    per_sp_counts = {
        sp: {
            "mafft": sum(1 for e in ev_m.get(sp, set()) if (e[0] == "stop" or e[2] % 3 != 0) and not e[3]),
            "muscle": sum(1 for e in ev_u.get(sp, set()) if (e[0] == "stop" or e[2] % 3 != 0) and not e[3]),
            "aligned_nt": len(mafft_msa[sp].replace("-", "")),
        }
        for sp in target_species
    }

    summary = {
        "dataset": dataset_id,
        "gene": gene,
        "exon": exon_label,
        "n_target_species_aligned": len(target_species),
        "target_species": sorted(target_species),
        "croc_control_lesions": croc_ctrl,
        "shared_mafft": len(s1),
        "shared_muscle": len(s2),
        "shared_robust_both": len(robust),
        "congruent_count": sum(1 for r in rows if r["tree_congruent"]),
        "incongruent_count": sum(1 for r in rows if not r["tree_congruent"]),
        "per_species_counts": per_sp_counts,
    }
    return rows, summary


def main() -> int:
    t0 = datetime.now(timezone.utc)
    L4.mkdir(parents=True, exist_ok=True)
    TABLES.mkdir(parents=True, exist_ok=True)

    meta = json.loads((OUT / "audit_metadata.json").read_text(encoding="utf-8"))
    intervals = load_intervals(meta)
    amis_seq, amis_win_start, _ = intervals["Amis"]

    print("=== [1/6] Extracting & verifying exact Alligator RefSeq CDS exons (ENAM & AMBN) ===")
    # ENAM XP_019346377.1 (8 CDS exons on '-' strand -> sort descending for 5'->3' transcription order)
    enam_gff_exons = sorted(meta["gff_audit"]["amis_cds_exons"]["XP_019346377.1"], key=lambda x: -x[1])
    ambn_gff_exons = sorted(meta["gff_audit"]["amis_cds_exons"]["XP_019346386.1"], key=lambda x: -x[1])

    enam_ref_exons: list[dict[str, Any]] = []
    cum_nt = 0
    for idx, (seqid, s, e, strand) in enumerate(enam_gff_exons, 1):
        loc_s = s - amis_win_start + 1
        loc_e = e - amis_win_start + 1
        dna = extract_subseq(amis_seq, loc_s, loc_e, strand)
        # Flanking intronic dinucleotides (in coding orientation!)
        # Upstream acceptor is right before 5' end (on '-' strand: loc_e+1..loc_e+2 revcomp'd)
        acc_di = revcomp(amis_seq[loc_e : loc_e + 2]) if strand == "-" else amis_seq[loc_s - 3 : loc_s - 1]
        don_di = revcomp(amis_seq[loc_s - 3 : loc_s - 1]) if strand == "-" else amis_seq[loc_e : loc_e + 2]
        enam_ref_exons.append({
            "exon_num": idx,
            "abs_start": s,
            "abs_end": e,
            "bp": len(dna),
            "cds_nt_start_0": cum_nt,
            "codon_start_1": (cum_nt // 3) + 1,
            "frame_shift": (3 - (cum_nt % 3)) % 3,
            "acceptor_2bp": acc_di,
            "donor_2bp": don_di,
            "dna": dna,
        })
        cum_nt += len(dna)
    enam_full_cds = "".join(x["dna"] for x in enam_ref_exons)
    enam_ref_prot = read_fasta(OUT / "queries" / "mp_ENAM_Amis.faa")["ENAM_Amis"]
    assert translate(enam_full_cds[:-3]) == enam_ref_prot, "ENAM CDS translation mismatch!"
    print(f"  ENAM verified: {len(enam_ref_exons)} exons, {len(enam_full_cds)} bp -> 100% match to XP_019346377.1 ({len(enam_ref_prot)} aa)")
    for ex in enam_ref_exons:
        print(f"    Exon {ex['exon_num']}: {ex['bp']:4d} bp (codons {ex['codon_start_1']:4d}..{(ex['cds_nt_start_0']+ex['bp'])//3:4d}, shift={ex['frame_shift']}) [{ex['acceptor_2bp']}..{ex['donor_2bp']}]")

    ambn_ref_exons: list[dict[str, Any]] = []
    cum_nt = 0
    for idx, (seqid, s, e, strand) in enumerate(ambn_gff_exons, 1):
        loc_s = s - amis_win_start + 1
        loc_e = e - amis_win_start + 1
        dna = extract_subseq(amis_seq, loc_s, loc_e, strand)
        acc_di = revcomp(amis_seq[loc_e : loc_e + 2]) if strand == "-" else amis_seq[loc_s - 3 : loc_s - 1]
        don_di = revcomp(amis_seq[loc_s - 3 : loc_s - 1]) if strand == "-" else amis_seq[loc_e : loc_e + 2]
        ambn_ref_exons.append({
            "exon_num": idx,
            "abs_start": s,
            "abs_end": e,
            "bp": len(dna),
            "cds_nt_start_0": cum_nt,
            "codon_start_1": (cum_nt // 3) + 1,
            "frame_shift": (3 - (cum_nt % 3)) % 3,
            "acceptor_2bp": acc_di,
            "donor_2bp": don_di,
            "dna": dna,
        })
        cum_nt += len(dna)
    ambn_full_cds = "".join(x["dna"] for x in ambn_ref_exons)
    ambn_ref_prot = read_fasta(OUT / "queries" / "mp_AMBN_Amis.faa")["AMBN_Amis"]
    assert translate(ambn_full_cds[:-3]) == ambn_ref_prot, "AMBN CDS translation mismatch!"
    print(f"  AMBN verified: {len(ambn_ref_exons)} exons, {len(ambn_full_cds)} bp -> 100% match to XP_019346386.1 ({len(ambn_ref_prot)} aa)")
    for ex in ambn_ref_exons:
        print(f"    Exon {ex['exon_num']}: {ex['bp']:4d} bp (codons {ex['codon_start_1']:4d}..{(ex['cds_nt_start_0']+ex['bp'])//3:4d}, shift={ex['frame_shift']}) [{ex['acceptor_2bp']}..{ex['donor_2bp']}]")

    # Audit Cpic ENAM's 58 kb miniprot jump
    print("\n=== [2/6] Auditing Cpic ENAM 58-kb miniprot intron jump & Arow ENAM Exon 6-7 splice sites ===")
    cpic_exons = parse_miniprot_exons(OUT / "miniprot" / "mp_ENAM_Amis_Cpic.paf")
    for idx, ex in enumerate(cpic_exons, 1):
        wstart = intervals["Cpic"][1]
        print(f"  Cpic ENAM mp_exon {idx}: q={ex['q_start_aa']}-{ex['q_end_aa']} bp={ex['bp']} abs={wstart+ex['t_lo_0']}-{wstart+ex['t_hi_0']-1} ({ex['strand']})")

    arow_exons = parse_miniprot_exons(OUT / "miniprot" / "mp_ENAM_Amis_Arow.paf")
    for idx, ex in enumerate(arow_exons, 1):
        aseq, wstart, _ = intervals["Arow"]
        # flanking 2bp splice dinucleotides on '-' strand
        lo, hi = ex["t_lo_0"], ex["t_hi_0"]
        if ex["strand"] == "-":
            acc_2 = revcomp(aseq[hi : hi + 2])
            don_2 = revcomp(aseq[lo - 2 : lo])
        else:
            acc_2 = aseq[lo - 2 : lo]
            don_2 = aseq[hi : hi + 2]
        print(f"  Arow ENAM mp_exon {idx}: q={ex['q_start_aa']}-{ex['q_end_aa']} bp={ex['bp']} abs={wstart+lo}-{wstart+hi-1} ({ex['strand']}) splice=[{acc_2}..{don_2}]")

    # Extract Crocodylian (Ggan, Cpor) coding exons directly from their miniprot PAFs
    croc_enam_ex8: dict[str, str] = {"Amis": enam_ref_exons[7]["dna"][:-3]}  # 2634 bp = 878 codons
    croc_ambn_exons: dict[int, dict[str, str]] = {ex["exon_num"]: {"Amis": ex["dna"][:-3] if ex["exon_num"] == 9 else ex["dna"]} for ex in ambn_ref_exons}

    for cr in ("Ggan", "Cpor"):
        cr_seq = intervals[cr][0]
        # ENAM exon 8 is the final (largest) exon in miniprot (>2400 bp)
        ex_list = parse_miniprot_exons(OUT / "miniprot" / f"mp_ENAM_Amis_{cr}.paf")
        big = max(ex_list, key=lambda x: x["bp"])
        croc_enam_ex8[cr] = extract_subseq(cr_seq, big["t_lo_0"] + 1, big["t_hi_0"], big["strand"])
        # AMBN exons (9 exons in crocodylians)
        a_list = parse_miniprot_exons(OUT / "miniprot" / f"mp_AMBN_Amis_{cr}.paf")
        if len(a_list) == 9:
            for idx, ex in enumerate(a_list, 1):
                croc_ambn_exons[idx][cr] = extract_subseq(cr_seq, ex["t_lo_0"] + 1, ex["t_hi_0"], ex["strand"])
        else:
            # Match by query start aa
            for ex in a_list:
                best_ref = min(ambn_ref_exons, key=lambda r: abs((r["cds_nt_start_0"] // 3) - ex["q_start_aa"]))
                croc_ambn_exons[best_ref["exon_num"]][cr] = extract_subseq(cr_seq, ex["t_lo_0"] + 1, ex["t_hi_0"], ex["strand"])

    # Helper to run dc-megablast / blastn of a query sequence against a list of species intervals, extract homologous spans, and run MAFFT + MUSCLE
    def build_and_align_dataset(
        ds_name: str,
        queries: dict[str, str],  # {query_label: dna_coding_strand}
        target_sps: list[str],
        croc_refs: dict[str, str],
        task: str = "dc-megablast",
        extra_flags: str = "",
        include_flank_bp: int = 0,
    ) -> tuple[dict[str, str], dict[str, str], dict[str, dict[str, Any]]]:
        ds_dir = L4 / ds_name
        ds_dir.mkdir(parents=True, exist_ok=True)
        qfa = ds_dir / "queries.fna"
        qfa.write_text("".join(f">{k}\n{v}\n" for k, v in queries.items()), encoding="utf-8", newline="\n")

        # Run blastn against each target species interval in WSL
        wsl_root = "/mnt/c/Users/user/bio-stack/runs/enamel_relics"
        lines = [
            "#!/usr/bin/env bash",
            "set -euo pipefail",
            f'Q="{wsl_root}/level4/{ds_name}/queries.fna"',
            f'O="{wsl_root}/level4/{ds_name}"',
            'FMT="6 qseqid sseqid pident length qstart qend sstart send evalue bitscore"',
        ]
        for sp in target_sps:
            lines.append(
                f'blastn -task {task} {extra_flags} -query "$Q" '
                f'-subject "{wsl_root}/intervals_enamel/interval_{sp}.fna" '
                f'-evalue 1e-4 -outfmt "$FMT" > "$O/blast_{sp}.tsv" || true'
            )
        run_wsl_cmd("\n".join(lines) + "\n", ds_dir / "run_blast.sh")

        hit_info: dict[str, dict[str, Any]] = {}
        msa_in_seqs: dict[str, str] = dict(croc_refs)

        for sp in target_sps:
            tsv = ds_dir / f"blast_{sp}.tsv"
            hits = []
            if tsv.exists():
                for c in csv.reader(tsv.open(encoding="utf-8"), delimiter="\t"):
                    if len(c) >= 10:
                        s1, s2 = int(c[6]), int(c[7])
                        hits.append({
                            "q": c[0], "pid": float(c[2]), "alen": int(c[3]),
                            "q1": int(c[4]), "q2": int(c[5]),
                            "s_lo": min(s1, s2), "s_hi": max(s1, s2),
                            "strand": "+" if s1 < s2 else "-",
                            "evalue": float(c[8]), "bits": float(c[9]),
                        })
            if not hits:
                continue
            # Pick best query + dominant strand cluster within 5 kb of top hit
            hits.sort(key=lambda h: -h["bits"])
            top = hits[0]
            cluster = [h for h in hits if h["strand"] == top["strand"] and abs(h["s_lo"] - top["s_lo"]) <= 5000]
            lo = min(h["s_lo"] for h in cluster)
            hi = max(h["s_hi"] for h in cluster)
            # Sanity cap against runaway spans (> 4x query length)
            max_qlen = max(len(v) for v in queries.values())
            if (hi - lo + 1) > max(4000, 3 * max_qlen):
                lo, hi = top["s_lo"], top["s_hi"]
            sp_seq, wstart, sacc = intervals[sp]
            dna = extract_subseq(sp_seq, lo - include_flank_bp, hi + include_flank_bp, top["strand"])
            msa_in_seqs[sp] = dna
            hit_info[sp] = {
                "query_used": top["q"],
                "seq_acc": sacc,
                "abs_start": wstart + lo - 1,
                "abs_end": wstart + hi - 1,
                "strand": top["strand"],
                "bp": hi - lo + 1,
                "top_pid": top["pid"],
                "top_bits": top["bits"],
                "top_evalue": top["evalue"],
                "n_hsps": len(cluster),
            }

        if len(msa_in_seqs) <= len(croc_refs) + 1:
            return {}, {}, hit_info

        in_fa = ds_dir / "msa_in.fna"
        in_fa.write_text("".join(f">{k}\n{v}\n" for k, v in msa_in_seqs.items()), encoding="utf-8", newline="\n")
        aln_sh = (
            "#!/usr/bin/env bash\nset -euo pipefail\n"
            f'O="{wsl_root}/level4/{ds_name}"\n'
            'mafft --localpair --maxiterate 1000 --thread 8 --quiet "$O/msa_in.fna" > "$O/mafft.fna"\n'
            'muscle -align "$O/msa_in.fna" -output "$O/muscle.fna" -threads 8 >/dev/null 2>&1\n'
        )
        run_wsl_cmd(aln_sh, ds_dir / "run_msa.sh")
        return read_fasta(ds_dir / "mafft.fna"), read_fasta(ds_dir / "muscle.fna"), hit_info

    all_shared_rows: list[dict[str, Any]] = []
    all_summaries: list[dict[str, Any]] = []
    all_splice_rows: list[dict[str, Any]] = []

    # ------------------------------------------------------------------
    # DATASET 1: Avian ENAM Exon 8 (full terminal coding exon, 2,634 bp)
    # ------------------------------------------------------------------
    print("\n=== [3/6] Dataset 1: Avian ENAM Exon 8 (full-length terminal exon, codons 201-1078) ===")
    shab_enam_exons = parse_miniprot_exons(OUT / "miniprot" / "mp_ENAM_Amis_Shab.paf")
    shab_ex8_span = (min(x["t_lo_0"] for x in shab_enam_exons), max(x["t_hi_0"] for x in shab_enam_exons), shab_enam_exons[0]["strand"])
    shab_ex8_dna = extract_subseq(intervals["Shab"][0], shab_ex8_span[0] + 1, shab_ex8_span[1], shab_ex8_span[2])
    arow_ex8 = max(arow_exons, key=lambda x: x["bp"])
    arow_ex8_dna = extract_subseq(intervals["Arow"][0], arow_ex8["t_lo_0"] + 1, arow_ex8["t_hi_0"], arow_ex8["strand"])

    m1, u1, h1 = build_and_align_dataset(
        "ENAM_exon8_aves",
        {"Shab_ENAM_ex8": shab_ex8_dna, "Arow_ENAM_ex8": arow_ex8_dna},
        AVES,
        croc_enam_ex8,
    )
    rows1, sum1 = score_dual_aligner_dataset(
        "ENAM_exon8_aves", "ENAM", "exon8", "Aves", m1, u1, AVIAN_CLADES,
        codon_offset_1based=enam_ref_exons[7]["codon_start_1"],
        reading_frame_shift=enam_ref_exons[7]["frame_shift"],
    )
    all_shared_rows.extend(rows1)
    all_summaries.append({**sum1, "blast_hits": h1})
    print(f"  Aligned birds ({sum1['n_target_species_aligned']}): {sum1['target_species']} | croc_ctrl={sum1['croc_control_lesions']} | "
          f"robust shared={sum1['shared_robust_both']} (congruent={sum1['congruent_count']}, incongruent={sum1['incongruent_count']})")
    for r in rows1:
        print(f"    codon {r['XP_019346377_or_386_codon']:4d} {r['mutation_type']:4s}({r['size_nt']:2d}nt) "
              f"carriers={r['carriers']:20s} without={r['covered_without']:15s} "
              f"anc={r['ancestral_croc_codon']}->{r['carrier_stop_codons']:8s} CpG={str(r['at_CpG_dinucleotide']):5s} {r['classification']}")

    # ------------------------------------------------------------------
    # DATASET 2: Avian ENAM Exon 6 + Intron 6 (GT->AT in Arow!) + Exon 7
    # ------------------------------------------------------------------
    print("\n=== [4/6] Dataset 2: Avian ENAM Exons 6-7 & Intron 6 GT->AT splice donor check ===")
    # Extract Arow Exon 6 + 60 bp flanking and Exon 7 + 60 bp flanking
    arow_ex6, arow_ex7 = arow_exons[0], arow_exons[1]
    arow_ex6_pad = extract_subseq(intervals["Arow"][0], arow_ex6["t_lo_0"] + 1 - 60, arow_ex6["t_hi_0"] + 60, arow_ex6["strand"])
    arow_ex7_pad = extract_subseq(intervals["Arow"][0], arow_ex7["t_lo_0"] + 1 - 60, arow_ex7["t_hi_0"] + 60, arow_ex7["strand"])
    # Also include Amis Exon 6 + 60 bp and Exon 7 + 60 bp
    ex6_s, ex6_e = enam_ref_exons[5]["abs_start"] - amis_win_start + 1, enam_ref_exons[5]["abs_end"] - amis_win_start + 1
    amis_ex6_pad = extract_subseq(amis_seq, ex6_s - 60, ex6_e + 60, "-")
    m2, u2, h2 = build_and_align_dataset(
        "ENAM_exon6_splice_aves",
        {"Arow_ENAM_ex6_pad": arow_ex6_pad},
        AVES,
        {"Amis": amis_ex6_pad},
        task="blastn",
        extra_flags="-word_size 7",
    )
    print(f"  Avian Exon 6 (+-60bp) hits: {list(h2.keys())}")
    if m2:
        # Locate exact canonical GT donor right after Exon 6 in Amis!
        # In amis_ex6_pad: 60 bp upstream acceptor + 48 bp Exon 6 + 60 bp downstream donor
        ref_al = m2["Amis"]
        # Find position of nucleotide index 60+48 = 108 and 109 (0-based) in ungapped amis_ex6_pad!
        col_of_p, p_idx = {}, -1
        for col_i, ch in enumerate(ref_al):
            if ch != "-":
                p_idx += 1
                col_of_p[p_idx] = col_i
        acc_cols = [col_of_p[58], col_of_p[59]]
        don_cols = [col_of_p[108], col_of_p[109]]
        for sp_k, seq_al in sorted(m2.items()):
            acc_sp = seq_al[acc_cols[0]] + seq_al[acc_cols[1]]
            don_sp = seq_al[don_cols[0]] + seq_al[don_cols[1]]
            print(f"    ENAM Exon 6 splice borders [{sp_k:4s}]: acceptor={acc_sp}  donor={don_sp}")
            all_splice_rows.append({
                "gene": "ENAM", "exon": "exon6", "clade_group": "Aves", "species_code": sp_k,
                "acceptor_2bp": acc_sp, "donor_2bp": don_sp,
                "acceptor_canonical": acc_sp == "AG", "donor_canonical": don_sp in ("GT", "GC"),
            })

    # ------------------------------------------------------------------
    # DATASET 3: Testudine ENAM Exon 8 (4 turtles + 3 crocodylians)
    # ------------------------------------------------------------------
    print("\n=== [5/6] Dataset 3: Testudine ENAM Exon 8 (4 turtles + 3 crocodylians) ===")
    cpic_ex8 = max(cpic_exons, key=lambda x: x["bp"])
    cpic_ex8_dna = extract_subseq(intervals["Cpic"][0], cpic_ex8["t_lo_0"] + 1, cpic_ex8["t_hi_0"], cpic_ex8["strand"])
    cmyd_enam_exons = parse_miniprot_exons(OUT / "miniprot" / "mp_ENAM_Amis_Cmyd.paf")
    cmyd_ex8 = max(cmyd_enam_exons, key=lambda x: x["bp"])
    cmyd_ex8_dna = extract_subseq(intervals["Cmyd"][0], cmyd_ex8["t_lo_0"] + 1, cmyd_ex8["t_hi_0"], cmyd_ex8["strand"])

    m3, u3, h3 = build_and_align_dataset(
        "ENAM_exon8_turtles",
        {"Cpic_ENAM_ex8": cpic_ex8_dna, "Cmyd_ENAM_ex8": cmyd_ex8_dna},
        TURTLES,
        croc_enam_ex8,
    )
    rows3, sum3 = score_dual_aligner_dataset(
        "ENAM_exon8_turtles", "ENAM", "exon8", "Testudines", m3, u3, TURTLE_CLADES,
        codon_offset_1based=enam_ref_exons[7]["codon_start_1"],
        reading_frame_shift=enam_ref_exons[7]["frame_shift"],
    )
    all_shared_rows.extend(rows3)
    all_summaries.append({**sum3, "blast_hits": h3})
    print(f"  Aligned turtles ({sum3['n_target_species_aligned']}): {sum3['target_species']} | croc_ctrl={sum3['croc_control_lesions']} | "
          f"robust shared={sum3['shared_robust_both']} (congruent={sum3['congruent_count']}, incongruent={sum3['incongruent_count']})")
    for r in rows3:
        print(f"    codon {r['XP_019346377_or_386_codon']:4d} {r['mutation_type']:4s}({r['size_nt']:2d}nt) "
              f"carriers={r['carriers']:20s} without={r['covered_without']:15s} "
              f"anc={r['ancestral_croc_codon']}->{r['carrier_stop_codons']:8s} CpG={str(r['at_CpG_dinucleotide']):5s} {r['classification']}")

    # ------------------------------------------------------------------
    # DATASET 4 & 5: AMBN Exon-by-Exon (Avian AND Testudine!)
    # ------------------------------------------------------------------
    print("\n=== [6/6] Datasets 4 & 5: Avian & Testudine AMBN coding exons (+ splice-site audit) ===")
    dnov_ambn_exons = parse_miniprot_exons(OUT / "miniprot" / "mp_AMBN_Amis_Dnov.paf")
    cmyd_ambn_exons = parse_miniprot_exons(OUT / "miniprot" / "mp_AMBN_Amis_Cmyd.paf")
    cpic_ambn_exons = parse_miniprot_exons(OUT / "miniprot" / "mp_AMBN_Amis_Cpic.paf")

    print(f"  Dnov AMBN miniprot exons: {[(e['q_start_aa'], e['q_end_aa'], e['bp']) for e in dnov_ambn_exons]}")
    print(f"  Cmyd AMBN miniprot exons: {[(e['q_start_aa'], e['q_end_aa'], e['bp']) for e in cmyd_ambn_exons]}")
    print(f"  Cpic AMBN miniprot exons: {[(e['q_start_aa'], e['q_end_aa'], e['bp']) for e in cpic_ambn_exons]}")

    # Map major AMBN exons: Exon 3 (117 bp, aa 33-71), Exon 4 (240 bp, aa 72-151), Exon 9 (525 bp, aa 216-390)
    # In Dnov: exon index 2 = aa 32-71 (117 bp), exon index 3 = aa 71-151 (~240 bp)
    # In Cmyd: last exon = terminal exon 9 (aa 216-390, ~525 bp)
    ambn_target_exons = [
        (3, "exon3", 32, 71),
        (4, "exon4", 71, 151),
        (9, "exon9", 215, 390),
    ]

    for ex_num, ex_lbl, qlo, qhi in ambn_target_exons:
        ref_ex = ambn_ref_exons[ex_num - 1]
        croc_sub = croc_ambn_exons[ex_num]
        # --- Avian AMBN exon ---
        # Find Dnov exon closest to [qlo, qhi]
        d_cand = [e for e in dnov_ambn_exons if min(e["q_end_aa"], qhi) - max(e["q_start_aa"], qlo) >= 20]
        if d_cand:
            d_ex = max(d_cand, key=lambda e: min(e["q_end_aa"], qhi) - max(e["q_start_aa"], qlo))
            d_dna = extract_subseq(intervals["Dnov"][0], d_ex["t_lo_0"] + 1, d_ex["t_hi_0"], d_ex["strand"])
            m_a, u_a, h_a = build_and_align_dataset(
                f"AMBN_{ex_lbl}_aves",
                {f"Dnov_AMBN_{ex_lbl}": d_dna},
                AVES,
                croc_sub,
                task="blastn" if len(d_dna) < 200 else "dc-megablast",
                extra_flags="-word_size 7" if len(d_dna) < 200 else "",
            )
            if m_a:
                r_a, s_a = score_dual_aligner_dataset(
                    f"AMBN_{ex_lbl}_aves", "AMBN", ex_lbl, "Aves", m_a, u_a, AVIAN_CLADES,
                    codon_offset_1based=ref_ex["codon_start_1"],
                    reading_frame_shift=ref_ex["frame_shift"],
                )
                all_shared_rows.extend(r_a)
                all_summaries.append({**s_a, "blast_hits": h_a})
                print(f"  [AMBN {ex_lbl} Aves] birds={s_a['target_species']} croc_ctrl={s_a['croc_control_lesions']} robust_shared={s_a['shared_robust_both']}")
                for r in r_a:
                    print(f"    codon {r['XP_019346377_or_386_codon']:4d} {r['mutation_type']:4s}({r['size_nt']:2d}nt) "
                          f"carriers={r['carriers']:20s} without={r['covered_without']:15s} "
                          f"anc={r['ancestral_croc_codon']}->{r['carrier_stop_codons']:8s} {r['classification']}")
                # Also inspect splice borders for every hit bird on this exon!
                for sp_h, hi in sorted(h_a.items()):
                    sp_seq, wstart, _ = intervals[sp_h]
                    loc_lo = hi["abs_start"] - wstart + 1
                    loc_hi = hi["abs_end"] - wstart + 1
                    if hi["strand"] == "-":
                        acc_2 = revcomp(sp_seq[loc_hi : loc_hi + 2])
                        don_2 = revcomp(sp_seq[loc_lo - 3 : loc_lo - 1])
                    else:
                        acc_2 = sp_seq[loc_lo - 3 : loc_lo - 1]
                        don_2 = sp_seq[loc_hi : loc_hi + 2]
                    all_splice_rows.append({
                        "gene": "AMBN", "exon": ex_lbl, "clade_group": "Aves", "species_code": sp_h,
                        "acceptor_2bp": acc_2, "donor_2bp": don_2,
                        "acceptor_canonical": acc_2 == "AG", "donor_canonical": don_2 in ("GT", "GC"),
                    })

        # --- Testudine AMBN exon ---
        t_queries: dict[str, str] = {}
        for tsp, texs in (("Cmyd", cmyd_ambn_exons), ("Cpic", cpic_ambn_exons)):
            tc = [e for e in texs if min(e["q_end_aa"], qhi) - max(e["q_start_aa"], qlo) >= 20]
            if tc:
                te = max(tc, key=lambda e: min(e["q_end_aa"], qhi) - max(e["q_start_aa"], qlo))
                t_queries[f"{tsp}_AMBN_{ex_lbl}"] = extract_subseq(intervals[tsp][0], te["t_lo_0"] + 1, te["t_hi_0"], te["strand"])
        if t_queries:
            m_t, u_t, h_t = build_and_align_dataset(
                f"AMBN_{ex_lbl}_turtles",
                t_queries,
                TURTLES,
                croc_sub,
                task="blastn" if min(len(v) for v in t_queries.values()) < 200 else "dc-megablast",
                extra_flags="-word_size 7" if min(len(v) for v in t_queries.values()) < 200 else "",
            )
            if m_t:
                r_t, s_t = score_dual_aligner_dataset(
                    f"AMBN_{ex_lbl}_turtles", "AMBN", ex_lbl, "Testudines", m_t, u_t, TURTLE_CLADES,
                    codon_offset_1based=ref_ex["codon_start_1"],
                    reading_frame_shift=ref_ex["frame_shift"],
                )
                all_shared_rows.extend(r_t)
                all_summaries.append({**s_t, "blast_hits": h_t})
                print(f"  [AMBN {ex_lbl} Turtles] turtles={s_t['target_species']} croc_ctrl={s_t['croc_control_lesions']} robust_shared={s_t['shared_robust_both']}")
                for r in r_t:
                    print(f"    codon {r['XP_019346377_or_386_codon']:4d} {r['mutation_type']:4s}({r['size_nt']:2d}nt) "
                          f"carriers={r['carriers']:20s} without={r['covered_without']:15s} "
                          f"anc={r['ancestral_croc_codon']}->{r['carrier_stop_codons']:8s} {r['classification']}")
                for sp_h, hi in sorted(h_t.items()):
                    sp_seq, wstart, _ = intervals[sp_h]
                    loc_lo = hi["abs_start"] - wstart + 1
                    loc_hi = hi["abs_end"] - wstart + 1
                    if hi["strand"] == "-":
                        acc_2 = revcomp(sp_seq[loc_hi : loc_hi + 2])
                        don_2 = revcomp(sp_seq[loc_lo - 3 : loc_lo - 1])
                    else:
                        acc_2 = sp_seq[loc_lo - 3 : loc_lo - 1]
                        don_2 = sp_seq[loc_hi : loc_hi + 2]
                    all_splice_rows.append({
                        "gene": "AMBN", "exon": ex_lbl, "clade_group": "Testudines", "species_code": sp_h,
                        "acceptor_2bp": acc_2, "donor_2bp": don_2,
                        "acceptor_canonical": acc_2 == "AG", "donor_canonical": don_2 in ("GT", "GC"),
                    })

    # Write clean tables & summary JSON
    if all_shared_rows:
        with (TABLES / "level4_shared_mutations.tsv").open("w", encoding="utf-8", newline="") as f:
            w = csv.DictWriter(f, fieldnames=list(all_shared_rows[0]), delimiter="\t", lineterminator="\n")
            w.writeheader()
            w.writerows(all_shared_rows)
    if all_splice_rows:
        with (TABLES / "level4_splice_borders.tsv").open("w", encoding="utf-8", newline="") as f:
            w = csv.DictWriter(f, fieldnames=list(all_splice_rows[0]), delimiter="\t", lineterminator="\n")
            w.writeheader()
            w.writerows(all_splice_rows)
    (OUT / "phase2_step4_level4_summary.json").write_text(json.dumps(all_summaries, indent=2), encoding="utf-8")

    # Append provenance
    prov_ev = {
        "ts": datetime.now(timezone.utc).isoformat(),
        "kind": "enamel-relics-phase2-step4-level4",
        "script": "scripts/enamel_relics/04_level4_shared_mutations.py",
        "n_datasets": len(all_summaries),
        "n_shared_mutations": len(all_shared_rows),
        "outputs": {
            "level4_shared_mutations.tsv": sha256_file(TABLES / "level4_shared_mutations.tsv") if (TABLES / "level4_shared_mutations.tsv").exists() else "",
            "level4_splice_borders.tsv": sha256_file(TABLES / "level4_splice_borders.tsv") if (TABLES / "level4_splice_borders.tsv").exists() else "",
        },
        "elapsed_s": round((datetime.now(timezone.utc) - t0).total_seconds(), 2),
    }
    for log_p in (OUT / "provenance_run.jsonl", ROOT / "runs" / "provenance.jsonl"):
        with log_p.open("a", encoding="utf-8") as f:
            f.write(json.dumps(prov_ev, ensure_ascii=False) + "\n")

    print(f"\nLevel 4 complete in {(datetime.now(timezone.utc) - t0).total_seconds():.1f}s -> {len(all_shared_rows)} robust shared mutations logged.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
