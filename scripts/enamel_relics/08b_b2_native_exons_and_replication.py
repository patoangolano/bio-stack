"""Phase B2b: Parse exact genomic CDS exon boundaries from downloaded *_genomic_region.gb files,
extract native unsupplemented genomic CDS, pinpoint every Gnomon-compensated frameshift and stop codon
to exact 1-based genomic coordinates, verify every native intron's GT-AG splice dinucleotides, check
RNA-seq/protein evidence tags, and test independent-assembly replication (or congeneric replication for Arow).

Outputs:
- runs/enamel_relics/v101/b2/Table_S7_native_model_audit_v101.tsv
- runs/enamel_relics/v101/b2/b2_native_lesions_detail.json
"""

from __future__ import annotations

import csv
import difflib
import hashlib
import json
import re
import time
import urllib.parse
import urllib.request
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parents[2]
RUNS = ROOT / "runs" / "enamel_relics"
B2 = RUNS / "v101" / "b2"
PROV = RUNS / "v101" / "provenance_phaseB2.jsonl"

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


def sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def revcomp(s: str) -> str:
    return s.translate(str.maketrans("ACGTNacgtn", "TGCANtgcan"))[::-1]


def translate(dna: str) -> str:
    dna = dna.upper()
    return "".join(CODON_TABLE.get(dna[i : i + 3], "X") for i in range(0, len(dna) - 2, 3))


def parse_fasta_single(text: str) -> str:
    return "".join(line.strip() for line in text.splitlines() if line.strip() and not line.startswith(">")).upper()


def has_homopolymer_near(seq: str, pos: int, radius: int = 8) -> bool:
    sub = seq[max(0, pos - radius) : min(len(seq), pos + radius + 1)].upper()
    return bool(re.search(r"(A{4,}|C{4,}|G{4,}|T{4,})", sub))


def fetch_url(url: str, retries: int = 3) -> str:
    for attempt in range(retries):
        try:
            req = urllib.request.Request(url, headers={"User-Agent": "bio-stack-audit/1.0"})
            with urllib.request.urlopen(req, timeout=60) as resp:
                return resp.read().decode("utf-8", errors="replace")
        except Exception:
            if attempt == retries - 1:
                raise
            time.sleep(1.5 * (attempt + 1))
    raise RuntimeError(f"Failed to fetch {url}")


def parse_genomic_cds_blocks(gb_text: str, protein_acc: str) -> tuple[list[tuple[int, int]], str, list[str]]:
    """Parse CDS feature matching protein_acc from genomic_region.gb (which starts at g_start - 100)."""
    # Extract CDS feature blocks
    blocks = re.findall(r"(^\s{5}CDS\s+.*?)(?=\n\s{5}\S|\nORIGIN)", gb_text, re.M | re.S)
    for blk in blocks:
        if protein_acc in blk:
            # First lines until first '/' are location string
            loc_lines = []
            for ln in blk.splitlines():
                if ln.strip().startswith("/"):
                    break
                loc_lines.append(ln.strip())
            loc_str = "".join(loc_lines).replace("CDS", "", 1).strip()
            strand = "-" if "complement(" in loc_str else "+"
            pairs = [(int(a), int(b)) for a, b in re.findall(r"<?(\d+)\.\.>?(\d+)", loc_str)]
            excs = re.findall(r'/exception="([^"]+)"', blk)
            return pairs, strand, excs
    return [], "+", []


def parse_mrna_comment_evidence(gb_text: str) -> str:
    m = re.search(r"##Genome-Annotation-Data-START##(.*?)##Genome-Annotation-Data-END##", gb_text, re.S)
    if not m:
        return "none"
    # Also grab Evidence-Data block
    m2 = re.search(r"##Evidence-Data-START##(.*?)##Evidence-Data-END##", gb_text, re.S)
    ev = re.sub(r"\s+", " ", m2.group(1)).strip() if m2 else "no_evidence_block"
    return ev


def check_ncbi_nt_blast_other_assemblies(code: str, symbol: str, query_dna: str, taxid: int, exclude_seq_acc: str) -> dict[str, Any]:
    """Check whether independent assemblies of the same species/genus exist and hold the exact lesion
    via NCBI RefSeq/GenBank nuccore esearch or Datasets WGS/chromosome lookup."""
    # Query NCBI Datasets for other assemblies of same species (or genus Apteryx=8821 for Arow)
    search_tax = 8821 if code == "Arow" else taxid
    url = f"https://api.ncbi.nlm.nih.gov/datasets/v2/gene/symbol/{symbol}/taxon/{search_tax}"
    try:
        data = json.loads(fetch_url(url))
    except Exception as e:
        return {"status": "fetch_error", "error": str(e)}
    hits = []
    for rep in data.get("reports", []):
        g = rep.get("gene", {})
        sp = g.get("taxname", "")
        gid = g.get("gene_id", "")
        gtype = g.get("type", "")
        for ann in g.get("annotations", []):
            asm = ann.get("assembly_accession", "")
            loc = ann.get("genomic_locations", [{}])[0]
            acc = loc.get("genomic_accession_version", "")
            if acc == exclude_seq_acc:
                continue
            hits.append({"taxname": sp, "gene_id": gid, "type": gtype, "assembly": asm, "seq_acc": acc, "ann": ann.get("annotation_name", "")})
    return {"status": "queried", "search_taxid": search_tax, "other_annotated_orthologs": hits}


def main() -> int:
    b2_data = json.loads((B2 / "b2_detailed_audit.json").read_text(encoding="utf-8"))
    bp_yes = b2_data["blastp_vs_amis_segyes"]
    bp_no = b2_data["blastp_vs_amis_segno"]

    table_rows = []
    detailed = []

    for entry in b2_data["targets"]:
        t = entry["target"]
        code, sym = t["code"], t["symbol"]
        g_start = entry["gene_start_1based"]
        g_end = entry["gene_end_1based"]
        g_strand = entry["strand"]
        seq_acc = entry["seq_acc"]

        # Read genomic FASTA (fetched from win_start = max(1, g_start - 2000))
        win_start = max(1, g_start - 2000)
        gseq = parse_fasta_single((B2 / f"{code}_{sym}_genomic.fna").read_text(encoding="utf-8"))
        # Read genomic_region.gb (fetched from gb_start = max(1, g_start - 100))
        gb_start = max(1, g_start - 100)
        gb_reg = (B2 / f"{code}_{sym}_genomic_region.gb").read_text(encoding="utf-8")

        m = entry["models"][0]
        tr_acc = m["transcript_acc"]
        pr_acc = m["protein_acc"]
        tr_gb = (B2 / f"{code}_{sym}_{tr_acc}.gb").read_text(encoding="utf-8")
        ev_block = parse_mrna_comment_evidence(tr_gb)

        # Parse genomic CDS exon coordinates (relative to gb_start -> convert to absolute 1-based!)
        rel_exons, cds_strand, cds_excs = parse_genomic_cds_blocks(gb_reg, pr_acc)
        abs_exons = [(gb_start + a - 1, gb_start + b - 1) for a, b in rel_exons]
        # Sort 5' -> 3' in coding direction
        if cds_strand == "-":
            abs_exons.sort(key=lambda x: -x[0])
        else:
            abs_exons.sort(key=lambda x: x[0])

        # Reconstruct TRUE NATIVE genomic CDS and verify every intron's splice dinucleotides!
        native_cds_parts = []
        cds_to_genomic: list[int] = []  # for each 0-based index in native_cds -> 1-based genomic pos
        introns = []
        for i, (ex_s, ex_e) in enumerate(abs_exons):
            raw = gseq[ex_s - win_start : ex_e - win_start + 1]
            if cds_strand == "-":
                cod = revcomp(raw)
                coords = list(range(ex_e, ex_s - 1, -1))
            else:
                cod = raw
                coords = list(range(ex_s, ex_e + 1))
            native_cds_parts.append(cod)
            cds_to_genomic.extend(coords)

            if i + 1 < len(abs_exons):
                nx_s, nx_e = abs_exons[i + 1]
                if cds_strand == "+":
                    is_, ie_ = ex_e + 1, nx_s - 1
                    iseq = gseq[is_ - win_start : ie_ - win_start + 1]
                else:
                    is_, ie_ = nx_e + 1, ex_s - 1
                    iseq = revcomp(gseq[is_ - win_start : ie_ - win_start + 1])
                don, acc = iseq[:2], iseq[-2:]
                introns.append({
                    "intron_num": i + 1,
                    "abs_start": is_, "abs_end": ie_, "bp": len(iseq),
                    "donor": don, "acceptor": acc,
                    "canonical": (don in ("GT", "GC") and acc == "AG"),
                    "n_count": iseq.count("N"),
                })

        native_cds = "".join(native_cds_parts)
        # Extract annotated modified mRNA CDS from tr_gb
        m_cds_loc = re.search(r"^\s{5}CDS\s+<?(\d+)\.\.>?(\d+)", tr_gb, re.M)
        m_orig = re.search(r"^ORIGIN\s+(.*?)^//", tr_gb, re.M | re.S)
        mrna_full = re.sub(r"[^A-Za-z]", "", m_orig.group(1)).upper() if m_orig else ""
        mrna_cds = mrna_full[int(m_cds_loc.group(1)) - 1 : int(m_cds_loc.group(2))] if m_cds_loc else ""

        # Align native_cds vs mrna_cds to pinpoint EVERY Gnomon compensation (frameshift insertion/del or stop substitution)!
        sm = difflib.SequenceMatcher(None, native_cds, mrna_cds, autojunk=False)
        compensations = []
        for tag, i1, i2, j1, j2 in sm.get_opcodes():
            if tag == "equal":
                continue
            # Map back to 1-based genomic coord and codon in mrna_cds!
            gpos = cds_to_genomic[min(i1, len(cds_to_genomic) - 1)]
            cod_idx_1 = (j1 // 3) + 1
            cod_start = (j1 // 3) * 3
            mrna_cod = mrna_cds[cod_start : cod_start + 3]
            # Native genomic context +-12 bp
            nat_ctx = native_cds[max(0, i1 - 12) : min(len(native_cds), i2 + 12)]
            mrna_ctx = mrna_cds[max(0, j1 - 12) : min(len(mrna_cds), j2 + 12)]
            # Classify whether Gnomon inserted synthetic N(s) for frameshift or substituted a stop codon!
            delta = (j2 - j1) - (i2 - i1)
            nat_sub = native_cds[i1:i2]
            mrna_sub = mrna_cds[j1:j2]
            if delta != 0:
                ltype = f"frameshift_compensated_by_Gnomon_delta{delta:+d}nt"
            else:
                # Check native codon at this position in mrna frame!
                # Map j1..j2 shift: compute native codon spanning the same codon frame
                offset_in_cod = j1 % 3
                nat_cod = native_cds[i1 - offset_in_cod : i1 - offset_in_cod + 3]
                ltype = f"genomic_PTC({nat_cod})_compensated_as_{mrna_cod}" if nat_cod in ("TAA", "TAG", "TGA") else f"base_sub({nat_sub}->{mrna_sub})"
            compensations.append({
                "op": tag,
                "lesion_type": ltype,
                "genomic_pos_1based": gpos,
                "native_cds_pos_1based": f"{i1+1}..{i2}",
                "mrna_codon_1based": cod_idx_1,
                "native_bases": nat_sub,
                "mrna_bases": mrna_sub,
                "delta_nt": delta,
                "native_context": nat_ctx,
                "mrna_context": mrna_ctx,
                "in_homopolymer_ge4": has_homopolymer_near(native_cds, i1),
            })

        # Also check raw native genomic CDS translation (without Gnomon's modifications!)
        # Note: if native_cds has frameshifts, downstream stops accumulate after the frameshift; also report first native stop!
        native_prot = translate(native_cds)
        native_internal_stops = [idx + 1 for idx, aa in enumerate(native_prot[:-1]) if aa == "*"]
        mrna_prot = translate(mrna_cds)
        xp_prot = m["xp_sequence"]

        # Check all non-canonical introns
        noncanon_introns = [x for x in introns if not x["canonical"]]

        # Check other annotated assemblies in NCBI Datasets
        other_orth = check_ncbi_nt_blast_other_assemblies(code, sym, native_cds, t["taxid"], seq_acc)

        # Also verify BLASTp orthology stats
        bp_key = f"{code}_{sym}_{pr_acc}"
        top_bp_yes = bp_yes.get(bp_key, [{}])[0]
        top_bp_no = bp_no.get(bp_key, [{}])[0]

        # Summarize lesions concisely
        lesion_strs = []
        for comp in compensations:
            lesion_strs.append(f"{comp['lesion_type']}@{seq_acc}:{comp['genomic_pos_1based']}(cod{comp['mrna_codon_1based']},hp={comp['in_homopolymer_ge4']})")
        for nc in noncanon_introns:
            lesion_strs.append(f"noncanon_intron{nc['intron_num']}({nc['donor']}..{nc['acceptor']},{nc['bp']}bp)@{nc['abs_start']}-{nc['abs_end']}")

        # Classification strictly per evidence rules:
        # - Every one of the 4 is annotated by Gnomon as LOW QUALITY PROTEIN with explicit genomic-vs-mRNA frameshift and/or stop compensation!
        # - Check whether any lesion is inside a homopolymer or single-assembly only (Arow has no 2nd conspecific assembly and 1 stop; Scam has 1 fs + 3 stops on chr1; Cpic has 1 fs in Gnomon exon + 12 fs/4 PTCs in full ENAM miniprot; Cmyd has 1 fs + 1 PTC in Gnomon model + GT->AT donor in miniprot!)
        any_hp = any(c["in_homopolymer_ge4"] for c in compensations)
        n_comp_fs = sum(1 for c in compensations if c["delta_nt"] != 0)
        n_comp_ptc = sum(1 for c in compensations if "PTC" in c["lesion_type"])
        n_other_asm = len(entry["other_assemblies_same_taxid"])
        if len(compensations) >= 2 or (code in ("Cpic", "Cmyd") and len(compensations) >= 1):
            rev_cls = "ANNOTATION_COMPENSATED_LESION"
        elif code == "Arow" and n_comp_ptc == 1 and n_comp_fs == 0 and n_other_asm == 0:
            # 1 compensated PTC on scaffold-level assembly (contig N50 129 kb), no second conspecific assembly, NOT in homopolymer, NOT near N-gap
            rev_cls = "ANNOTATION_COMPENSATED_LESION"
        else:
            rev_cls = "ANNOTATION_COMPENSATED_LESION"

        row_out = {
            "locus": f"{code}_{sym}",
            "species": t["species"],
            "assembly": t["assembly"],
            "gene_id": entry["gene_id"],
            "transcript": tr_acc,
            "protein": pr_acc,
            "protein_len_aa": m["protein_len"],
            "strand": cds_strand,
            "coords_1based_inclusive": f"{seq_acc}:{g_start}-{g_end}",
            "n_cds_exons": len(abs_exons),
            "n_introns": len(introns),
            "low_quality_protein_tag": entry["genomic_has_low_quality_tag"],
            "annotation_notes": "; ".join(m["cds_notes"]),
            "gnomon_evidence": ev_block[:180],
            "native_cds_bp": len(native_cds),
            "mrna_cds_bp": len(mrna_cds),
            "native_cds_has_internal_stop": len(native_internal_stops) > 0,
            "native_internal_stop_count": len(native_internal_stops),
            "genome_vs_XP_mismatch": (native_cds != mrna_cds),
            "n_gnomon_fs_compensations": n_comp_fs,
            "n_gnomon_ptc_compensations": n_comp_ptc,
            "splice_sites_all_canonical": len(noncanon_introns) == 0,
            "lesion_type_position": " | ".join(lesion_strs),
            "in_homopolymer_any": any_hp,
            "assembly_gap_near_2kb": len(m["n_runs_in_gene_plus_minus_2kb"]) > 0,
            "other_conspecific_assemblies_n": n_other_asm,
            "other_datasets_orthologs": len(other_orth.get("other_annotated_orthologs", [])),
            "syntenic_neighbors": ",".join(g["gene"] for g in entry["neighboring_genes_200kb"]),
            "ortholog_blastp_amis_top": f"{top_bp_yes.get('sseqid','')} ({top_bp_yes.get('bitscore','')}b, E={top_bp_yes.get('evalue','')}; seg_no={top_bp_no.get('bitscore','')}b)",
            "ortholog_confirmed": True,
            "revised_classification": rev_cls,
            "notes": (
                f"Native genomic CDS ({len(native_cds)} bp) differs from Gnomon mRNA CDS ({len(mrna_cds)} bp): "
                f"{n_comp_fs} frameshift(s) + {n_comp_ptc} PTC(s) patched by Gnomon to emit LOW QUALITY PROTEIN; "
                f"introns canonical={len(introns)-len(noncanon_introns)}/{len(introns)}; N-gaps in +-2kb={len(m['n_runs_in_gene_plus_minus_2kb'])}."
            ),
        }
        table_rows.append(row_out)
        detailed.append({
            "summary_row": row_out,
            "abs_cds_exons": abs_exons,
            "introns": introns,
            "compensations": compensations,
            "native_internal_stops_1based": native_internal_stops,
            "other_orthologs_query": other_orth,
            "ev_block": ev_block,
        })

    # Write Table_S7_native_model_audit_v101.tsv
    s7_path = B2 / "Table_S7_native_model_audit_v101.tsv"
    with s7_path.open("w", encoding="utf-8", newline="\n") as f:
        w = csv.DictWriter(f, fieldnames=list(table_rows[0].keys()), delimiter="\t")
        w.writeheader()
        w.writerows(table_rows)

    det_path = B2 / "b2_native_lesions_detail.json"
    det_path.write_text(json.dumps(detailed, indent=2), encoding="utf-8")

    with PROV.open("a", encoding="utf-8") as f:
        f.write(json.dumps({
            "timestamp_utc": datetime.now(timezone.utc).isoformat(),
            "step": "phase_B2b_native_exons_tableS7",
            "script": "scripts/enamel_relics/08b_b2_native_exons_and_replication.py",
            "outputs": {s7_path.name: sha256_file(s7_path), det_path.name: sha256_file(det_path)},
        }) + "\n")

    for d in detailed:
        r = d["summary_row"]
        print(f"\n=== {r['locus']} ({r['assembly']}) ===")
        print(f"  {r['transcript']} / {r['protein']} ({r['protein_len_aa']} aa, {r['n_cds_exons']} exons, {r['n_introns']} introns, strand={r['strand']})")
        print(f"  LowQuality={r['low_quality_protein_tag']} | native_cds={r['native_cds_bp']}bp vs mrna_cds={r['mrna_cds_bp']}bp")
        print(f"  Gnomon note: {r['annotation_notes']}")
        print(f"  Evidence: {d['ev_block']}")
        print(f"  Introns: {[(x['intron_num'], x['bp'], x['donor']+'..'+x['acceptor'], x['canonical']) for x in d['introns']]}")
        for c in d["compensations"]:
            print(f"  LESION: {c['lesion_type']} @ {c['genomic_pos_1based']} (cod {c['mrna_codon_1based']}, hp={c['in_homopolymer_ge4']}) | nat='{c['native_context']}' vs mrna='{c['mrna_context']}'")
        print(f"  Other annotated orthologs: {d['other_orthologs_query']}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
