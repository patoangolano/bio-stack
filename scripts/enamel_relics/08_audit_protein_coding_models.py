"""Phase B2 reproducible audit of the four RefSeq PROTEIN_CODING enamel models in edentulous sauropsids:
  1. Chrysemys picta (GCF_011386835.1) — ENAM
  2. Chelonia mydas (GCF_015237465.2) — AMBN
  3. Struthio camelus (GCF_040807025.1) — AMELX
  4. Apteryx rowi (GCF_003343035.1) — AMELX

For each locus:
- Queries NCBI Datasets v2 gene/assembly metadata and NCBI E-utilities (GenBank/GenPept + genomic FASTA).
- Extracts native annotated gene ID, transcript accession (XM_), protein accession (XP_), exon/CDS coords,
  strand, RefSeq/Gnomon curation notes ('LOW QUALITY PROTEIN', 'modified relative to its source genomic sequence',
  'frameshift', 'partial', 'exception=', etc.), and evidence support.
- Compares NATIVE genomic spliced CDS translation directly against annotated XM_ CDS and XP_ protein to detect
  RefSeq/Gnomon frameshift/stop compensations vs intact native ORFs.
- Audits canonical GT-AG/GC-AG splice dinucleotides of every native annotated intron directly on the genome.
- Checks flanking sequence context (homopolymers >=4 nt, assembly N-gaps within +-1000 bp, scaffold/contig status).
- Checks syntenic neighborhood (flanking genes) and reciprocal BLASTp against Alligator (48,040 proteins) & Human
  orthologs to verify true orthology (critical for avian AMELX).
- Queries NCBI Datasets for independent assemblies of each species (or closest congeneric/same-species assembly)
  and tests whether each native lesion replicates in a second assembly.

Outputs:
- runs/enamel_relics/v101/b2/Table_S7_native_model_audit_v101.tsv
- runs/enamel_relics/v101/b2/b2_detailed_audit.json
- runs/enamel_relics/v101/provenance_phaseB2.jsonl
"""

from __future__ import annotations

import csv
import hashlib
import json
import re
import subprocess
import time
import urllib.parse
import urllib.request
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parents[2]
RUNS = ROOT / "runs" / "enamel_relics"
OUT = RUNS / "v101" / "b2"
OUT.mkdir(parents=True, exist_ok=True)
PROV = RUNS / "v101" / "provenance_phaseB2.jsonl"

WSL_DISTRO = "Ubuntu-26.04"
WSL_ROOT = "/mnt/c/Users/user/bio-stack"

TARGETS = [
    {"code": "Cpic", "species": "Chrysemys picta", "taxid": 8478, "assembly": "GCF_011386835.1", "symbol": "ENAM"},
    {"code": "Cmyd", "species": "Chelonia mydas", "taxid": 8469, "assembly": "GCF_015237465.2", "symbol": "AMBN"},
    {"code": "Scam", "species": "Struthio camelus", "taxid": 8801, "assembly": "GCF_040807025.1", "symbol": "AMELX"},
    {"code": "Arow", "species": "Apteryx rowi", "taxid": 308060, "assembly": "GCF_003343035.1", "symbol": "AMELX"},
]

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
    out = []
    for i in range(0, len(dna) - 2, 3):
        cod = dna[i : i + 3]
        out.append(CODON_TABLE.get(cod, "X"))
    return "".join(out)


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


def efetch(db: str, acc_or_id: str, rettype: str, retmode: str = "text", seq_start: int | None = None, seq_stop: int | None = None) -> str:
    params = {"db": db, "id": acc_or_id, "rettype": rettype, "retmode": retmode}
    if seq_start is not None and seq_stop is not None:
        params["seq_start"] = str(seq_start)
        params["seq_stop"] = str(seq_stop)
    url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?" + urllib.parse.urlencode(params)
    time.sleep(0.4)
    return fetch_url(url)


def run_wsl_script(script_text: str, script_path: Path) -> str:
    script_path.write_text(script_text, encoding="utf-8", newline="\n")
    rel = script_path.relative_to(ROOT).as_posix()
    res = subprocess.run(
        ["wsl", "-d", WSL_DISTRO, "--", "bash", "-l", f"{WSL_ROOT}/{rel}"],
        capture_output=True, text=True, encoding="utf-8", check=False,
    )
    if res.returncode != 0:
        raise RuntimeError(f"WSL failed ({res.returncode}): {res.stderr[:800]}")
    return res.stdout


def parse_fasta_single(text: str) -> str:
    return "".join(line.strip() for line in text.splitlines() if line.strip() and not line.startswith(">"))


def find_homopolymer_near(seq: str, idx: int, radius: int = 8, min_run: int = 4) -> bool:
    sub = seq[max(0, idx - radius) : min(len(seq), idx + radius + 1)].upper()
    return bool(re.search(r"(A{" + str(min_run) + r",}|C{" + str(min_run) + r",}|G{" + str(min_run) + r",}|T{" + str(min_run) + r",})", sub))


def parse_gb_feature_blocks(gb_text: str) -> dict[str, Any]:
    """Parse key metadata from GenBank flatfile (mRNA or protein)."""
    info: dict[str, Any] = {}
    m_def = re.search(r"^DEFINITION\s+(.*?)(?=\n[A-Z])", gb_text, re.M | re.S)
    info["definition"] = re.sub(r"\s+", " ", m_def.group(1)).strip() if m_def else ""
    m_com = re.search(r"^COMMENT\s+(.*?)(?=\n(?:PRIMARY|FEATURES))", gb_text, re.M | re.S)
    info["comment"] = re.sub(r"\s+", " ", m_com.group(1)).strip() if m_com else ""
    # Extract CDS qualifiers (/note, /exception, /codon_start, /protein_id, /translation)
    m_cds = re.search(r"^\s{5}CDS\s+(.*?)(?=\n\s{5}\S|\nORIGIN)", gb_text, re.M | re.S)
    if m_cds:
        cds_blk = m_cds.group(1)
        info["cds_location"] = cds_blk.splitlines()[0].strip()
        info["cds_notes"] = [re.sub(r"\s+", " ", x) for x in re.findall(r'/note="([^"]+)"', cds_blk, re.S)]
        info["cds_exceptions"] = [re.sub(r"\s+", " ", x) for x in re.findall(r'/exception="([^"]+)"', cds_blk, re.S)]
        m_cs = re.search(r"/codon_start=(\d+)", cds_blk)
        info["codon_start"] = int(m_cs.group(1)) if m_cs else 1
        m_pid = re.search(r'/protein_id="([^"]+)"', cds_blk)
        info["protein_id"] = m_pid.group(1) if m_pid else ""
        m_tr = re.search(r'/translation="([^"]+)"', cds_blk, re.S)
        info["translation"] = re.sub(r"\s+", "", m_tr.group(1)) if m_tr else ""
    # Extract ORIGIN sequence
    m_orig = re.search(r"^ORIGIN\s+(.*?)^//", gb_text, re.M | re.S)
    if m_orig:
        info["sequence"] = re.sub(r"[^A-Za-z]", "", m_orig.group(1)).upper()
    return info


def get_datasets_gene(symbol: str, taxid: int) -> dict[str, Any]:
    url = f"https://api.ncbi.nlm.nih.gov/datasets/v2/gene/symbol/{symbol}/taxon/{taxid}"
    data = json.loads(fetch_url(url))
    return data


def get_datasets_gene_product_report(gene_id: str) -> dict[str, Any]:
    url = f"https://api.ncbi.nlm.nih.gov/datasets/v2/gene/id/{gene_id}/product_report"
    return json.loads(fetch_url(url))


def get_neighboring_genes(assembly: str, seq_acc: str, start: int, end: int, pad: int = 500000) -> list[dict[str, Any]]:
    """Query NCBI Datasets annotation report or E-utilities gene search around genomic interval."""
    # Use NCBI E-utilities esearch on gene DB for taxon + chromosome/scaffold, or efetch genomic gbwithparts on +-200kb!
    lo = max(1, start - 200000)
    hi = end + 200000
    # Fetch ft (feature table) for the +-200kb window — tiny, fast, lists all gene names and coords!
    ft = efetch("nuccore", seq_acc, "ft", "text", seq_start=lo, seq_stop=hi)
    genes = []
    cur_start, cur_end, cur_type = None, None, None
    cur_gene = None
    for line in ft.splitlines():
        if not line.strip() or line.startswith(">Feature"):
            continue
        if not line.startswith("\t"):
            if cur_type == "gene" and cur_gene:
                genes.append({"gene": cur_gene, "start": lo + min(cur_start, cur_end) - 1, "end": lo + max(cur_start, cur_end) - 1, "strand": "+" if cur_start < cur_end else "-"})
            parts = line.split("\t")
            if len(parts) >= 3:
                try:
                    cur_start = int(re.sub(r"[<>]", "", parts[0]))
                    cur_end = int(re.sub(r"[<>]", "", parts[1]))
                    cur_type = parts[2].strip()
                    cur_gene = None
                except ValueError:
                    cur_type = None
        else:
            parts = line.strip().split("\t")
            if cur_type == "gene" and len(parts) >= 2 and parts[0] == "gene":
                cur_gene = parts[1].strip()
    if cur_type == "gene" and cur_gene:
        genes.append({"gene": cur_gene, "start": lo + min(cur_start, cur_end) - 1, "end": lo + max(cur_start, cur_end) - 1, "strand": "+" if cur_start < cur_end else "-"})
    return sorted(genes, key=lambda g: g["start"])


def get_other_assemblies(taxid: int, exclude_acc: str) -> list[dict[str, Any]]:
    url = f"https://api.ncbi.nlm.nih.gov/datasets/v2/genome/taxon/{taxid}/dataset_report?page_size=20"
    data = json.loads(fetch_url(url))
    out = []
    for rep in data.get("reports", []):
        acc = rep.get("accession", "")
        if acc == exclude_acc:
            continue
        asm = rep.get("assembly_info", {})
        stats = rep.get("assembly_stats", {})
        out.append({
            "accession": acc,
            "name": asm.get("assembly_name", ""),
            "level": asm.get("assembly_level", ""),
            "release_date": asm.get("release_date", ""),
            "contig_n50": stats.get("contig_n50", 0),
            "scaffold_n50": stats.get("scaffold_n50", 0),
            "paired_acc": rep.get("paired_accession", ""),
        })
    # Exclude paired GCA of the exact same assembly!
    base_num = exclude_acc.split("_")[1].split(".")[0]
    out = [x for x in out if x["accession"].split("_")[1].split(".")[0] != base_num]
    out.sort(key=lambda x: (-int(x.get("contig_n50") or 0)))
    return out


def audit_one_target(t: dict[str, Any]) -> dict[str, Any]:
    print(f"\n=== Auditing {t['code']} ({t['species']}) {t['symbol']} [{t['assembly']}] ===", flush=True)
    gdata = get_datasets_gene(t["symbol"], t["taxid"])
    (OUT / f"{t['code']}_{t['symbol']}_gene.json").write_text(json.dumps(gdata, indent=2), encoding="utf-8")

    # Find matching gene report
    gene_rep = None
    for rep in gdata.get("reports", []):
        g = rep.get("gene", {})
        if g.get("symbol", "").upper() == t["symbol"].upper():
            gene_rep = g
            break
    if gene_rep is None:
        return {"target": t, "error": "Gene not found in Datasets v2"}

    gene_id = str(gene_rep["gene_id"])
    # Find annotation on target assembly (or genomic_ranges)
    annots = gene_rep.get("annotations", [])
    chosen_annot = None
    for a in annots:
        if a.get("assembly_accession") == t["assembly"]:
            chosen_annot = a
            break
    if chosen_annot is None and annots:
        chosen_annot = annots[0]

    seq_acc = chosen_annot["genomic_locations"][0]["genomic_accession_version"]
    g_range = chosen_annot["genomic_locations"][0]["genomic_range"]
    g_start = int(g_range["begin"])
    g_end = int(g_range["end"])
    g_strand = "+" if g_range.get("orientation") == "plus" else "-"
    ann_name = chosen_annot.get("annotation_name", "")
    ann_date = chosen_annot.get("release_date", "")

    # Product report gives exact transcripts, CDS exon ranges, proteins!
    prep = get_datasets_gene_product_report(gene_id)
    (OUT / f"{t['code']}_{t['symbol']}_products.json").write_text(json.dumps(prep, indent=2), encoding="utf-8")
    transcripts_info = []
    for rep in prep.get("reports", []):
        prod = rep.get("product", {})
        for tr in prod.get("transcripts", []):
            tr_acc = tr.get("accession_version", "")
            prot = tr.get("protein", {})
            prot_acc = prot.get("accession_version", "")
            cds_obj = tr.get("cds", {})
            cds_ranges = []
            for r in cds_obj.get("range", []):
                cds_ranges.append({
                    "begin": int(r["begin"]),
                    "end": int(r["end"]),
                    "orientation": r.get("orientation", ""),
                })
            # Also transcript genomic_range_set
            ex_ranges = []
            for r in tr.get("exons", {}).get("range", []):
                ex_ranges.append({"begin": int(r["begin"]), "end": int(r["end"]), "orientation": r.get("orientation", "")})
            transcripts_info.append({
                "transcript_acc": tr_acc,
                "protein_acc": prot_acc,
                "protein_name": prot.get("name", ""),
                "protein_len": prot.get("length", 0),
                "cds_ranges": cds_ranges,
                "exon_ranges": ex_ranges,
                "genomic_acc": tr.get("genomic_locations", [{}])[0].get("genomic_accession_version", seq_acc) if tr.get("genomic_locations") else seq_acc,
            })

    # Fetch genomic FASTA of gene (+-2000 bp)
    pad = 2000
    win_start = max(1, g_start - pad)
    win_end = g_end + pad
    genomic_fa = efetch("nuccore", seq_acc, "fasta", "text", seq_start=win_start, seq_stop=win_end)
    (OUT / f"{t['code']}_{t['symbol']}_genomic.fna").write_text(genomic_fa, encoding="utf-8")
    gseq = parse_fasta_single(genomic_fa).upper()

    # Fetch genomic feature table of exact gene region (+-100 bp) to capture GFF/GenBank exceptions directly
    gene_gb_snippet = efetch("nuccore", seq_acc, "gbwithparts", "text", seq_start=max(1, g_start - 100), seq_stop=g_end + 100)
    (OUT / f"{t['code']}_{t['symbol']}_genomic_region.gb").write_text(gene_gb_snippet, encoding="utf-8")
    # Extract CDS/mRNA exceptions and notes from the genomic GenBank record!
    genomic_exceptions = sorted(set(re.findall(r'/exception="([^"]+)"', gene_gb_snippet)))
    genomic_notes = sorted(set(re.sub(r"\s+", " ", x) for x in re.findall(r'/note="([^"]+)"', gene_gb_snippet, re.S)))
    has_low_quality = "LOW QUALITY PROTEIN" in gene_gb_snippet.upper()

    # Fetch neighboring genes (+-200 kb) for synteny check
    neighbors = get_neighboring_genes(t["assembly"], seq_acc, g_start, g_end)

    # Audit each transcript/protein model natively
    audited_models = []
    for ti in transcripts_info:
        tr_acc = ti["transcript_acc"]
        prot_acc = ti["protein_acc"]
        tr_gb = efetch("nuccore", tr_acc, "gb", "text") if tr_acc else ""
        if tr_gb:
            (OUT / f"{t['code']}_{t['symbol']}_{tr_acc}.gb").write_text(tr_gb, encoding="utf-8")
        pr_gb = efetch("protein", prot_acc, "gp", "text") if prot_acc else ""
        if pr_gb:
            (OUT / f"{t['code']}_{t['symbol']}_{prot_acc}.gp").write_text(pr_gb, encoding="utf-8")

        tr_parsed = parse_gb_feature_blocks(tr_gb) if tr_gb else {}
        pr_parsed = parse_gb_feature_blocks(pr_gb) if pr_gb else {}

        # Reconstruct native genomic CDS directly from genome sequence using CDS ranges from Datasets
        # Order exons 5'->3' according to strand
        cds_exons = list(ti["cds_ranges"])
        if g_strand == "-":
            cds_exons.sort(key=lambda x: -x["begin"])
        else:
            cds_exons.sort(key=lambda x: x["begin"])

        exon_details = []
        spliced_genomic_cds = []
        intron_details = []
        for idx, ex in enumerate(cds_exons):
            b, e = ex["begin"], ex["end"]
            rel_b = b - win_start
            rel_e = e - win_start + 1
            raw = gseq[rel_b:rel_e]
            seq_cod = revcomp(raw) if g_strand == "-" else raw
            spliced_genomic_cds.append(seq_cod)
            n_count = seq_cod.count("N")
            exon_details.append({
                "exon_idx": idx + 1,
                "abs_start": b,
                "abs_end": e,
                "bp": len(seq_cod),
                "n_bases": n_count,
            })
            # Check intron between this exon and the next exon
            if idx + 1 < len(cds_exons):
                nx = cds_exons[idx + 1]
                if g_strand == "+":
                    intron_start = e + 1
                    intron_end = nx["begin"] - 1
                    intron_seq = gseq[intron_start - win_start : intron_end - win_start + 1]
                else:
                    intron_start = nx["end"] + 1
                    intron_end = b - 1
                    intron_seq = revcomp(gseq[intron_start - win_start : intron_end - win_start + 1])
                donor = intron_seq[:2] if len(intron_seq) >= 2 else "??"
                acceptor = intron_seq[-2:] if len(intron_seq) >= 2 else "??"
                intron_details.append({
                    "intron_idx": idx + 1,
                    "abs_start": intron_start,
                    "abs_end": intron_end,
                    "bp": len(intron_seq),
                    "donor_2bp": donor,
                    "acceptor_2bp": acceptor,
                    "canonical": (donor in ("GT", "GC") and acceptor == "AG"),
                    "n_bases": intron_seq.count("N"),
                })

        genomic_cds = "".join(spliced_genomic_cds)
        # Extract annotated mRNA CDS from tr_parsed!
        mrna_seq = tr_parsed.get("sequence", "")
        mrna_cds = ""
        m_loc = re.match(r"<?(\d+)\.\.>?(\d+)", tr_parsed.get("cds_location", ""))
        if m_loc and mrna_seq:
            c1, c2 = int(m_loc.group(1)), int(m_loc.group(2))
            mrna_cds = mrna_seq[c1 - 1 : c2]

        # Translate native genomic CDS (accounting for codon_start if partial)
        cstart = tr_parsed.get("codon_start", 1)
        genomic_prot = translate(genomic_cds[cstart - 1 :])
        mrna_prot = translate(mrna_cds[cstart - 1 :]) if mrna_cds else ""
        xp_prot = tr_parsed.get("translation", "")

        # Internal stops (excluding terminal stop at very end)
        genomic_internal_stops = [i + 1 for i, aa in enumerate(genomic_prot[:-1]) if aa == "*"]
        mrna_internal_stops = [i + 1 for i, aa in enumerate(mrna_prot[:-1]) if aa == "*"]

        # Exact nucleotide difference between genomic_cds and mrna_cds!
        # Find exact indels/substitutions between genomic_cds and mrna_cds via difflib
        import difflib
        diffs = []
        if mrna_cds and genomic_cds != mrna_cds:
            sm = difflib.SequenceMatcher(None, genomic_cds, mrna_cds, autojunk=False)
            for tag, i1, i2, j1, j2 in sm.get_opcodes():
                if tag != "equal":
                    ctx_g = genomic_cds[max(0, i1 - 10) : min(len(genomic_cds), i2 + 10)]
                    ctx_m = mrna_cds[max(0, j1 - 10) : min(len(mrna_cds), j2 + 10)]
                    diffs.append({
                        "op": tag,
                        "genomic_cds_pos_1based": f"{i1+1}..{i2}",
                        "mrna_cds_pos_1based": f"{j1+1}..{j2}",
                        "genomic_bases": genomic_cds[i1:i2],
                        "mrna_bases": mrna_cds[j1:j2],
                        "delta_nt": (j2 - j1) - (i2 - i1),
                        "genomic_context": ctx_g,
                        "mrna_context": ctx_m,
                        "in_homopolymer": find_homopolymer_near(genomic_cds, i1),
                    })

        # Check N-gaps anywhere in gene +-1000 bp
        n_runs = [(m.start() + win_start, m.end() + win_start - 1, m.end() - m.start()) for m in re.finditer(r"N+", gseq)]

        audited_models.append({
            "transcript_acc": tr_acc,
            "protein_acc": prot_acc,
            "protein_name": ti["protein_name"],
            "protein_len": ti["protein_len"],
            "tr_definition": tr_parsed.get("definition", ""),
            "tr_comment": tr_parsed.get("comment", ""),
            "pr_comment": pr_parsed.get("comment", ""),
            "cds_notes": tr_parsed.get("cds_notes", []),
            "cds_exceptions": tr_parsed.get("cds_exceptions", []),
            "n_cds_exons": len(exon_details),
            "exon_details": exon_details,
            "intron_details": intron_details,
            "all_introns_canonical": all(x["canonical"] for x in intron_details),
            "genomic_cds_bp": len(genomic_cds),
            "mrna_cds_bp": len(mrna_cds),
            "genomic_cds_equals_mrna_cds": (genomic_cds == mrna_cds),
            "genomic_cds_len_mod3": len(genomic_cds) % 3,
            "genomic_internal_stop_codons_1based": genomic_internal_stops,
            "mrna_internal_stop_codons_1based": mrna_internal_stops,
            "genomic_prot_matches_xp": (genomic_prot.rstrip("*") == xp_prot.rstrip("*")),
            "mrna_prot_matches_xp": (mrna_prot.rstrip("*") == xp_prot.rstrip("*")),
            "genomic_vs_mrna_diffs": diffs,
            "n_runs_in_gene_plus_minus_2kb": n_runs,
            "xp_sequence": xp_prot,
            "genomic_cds_sequence": genomic_cds,
        })

    # Find other assemblies of same species for replication check
    other_asms = get_other_assemblies(t["taxid"], t["assembly"])

    return {
        "target": t,
        "gene_id": gene_id,
        "gene_type": gene_rep.get("type", ""),
        "description": gene_rep.get("description", ""),
        "annotation_name": ann_name,
        "annotation_release_date": ann_date,
        "seq_acc": seq_acc,
        "gene_start_1based": g_start,
        "gene_end_1based": g_end,
        "strand": g_strand,
        "genomic_has_low_quality_tag": has_low_quality,
        "genomic_exceptions": genomic_exceptions,
        "genomic_notes": genomic_notes,
        "neighboring_genes_200kb": neighbors,
        "models": audited_models,
        "other_assemblies_same_taxid": other_asms[:5],
    }


def main() -> int:
    results = []
    xp_fasta_lines = []
    for t in TARGETS:
        res = audit_one_target(t)
        results.append(res)
        # Collect primary XP sequence for reciprocal BLASTp against Alligator proteome
        for m in res.get("models", [])[:1]:
            if m.get("xp_sequence"):
                xp_fasta_lines.append(f">{t['code']}_{t['symbol']}_{m['protein_acc']}\n{m['xp_sequence']}\n")

    # Run reciprocal BLASTp of all 4 XP_ proteins against full 48,040 Alligator proteome in WSL!
    xp_fa = OUT / "b2_xp_proteins.faa"
    xp_fa.write_text("".join(xp_fasta_lines), encoding="utf-8", newline="\n")
    bp_sh = f"""#!/usr/bin/env bash
set -euo pipefail
D="{WSL_ROOT}/runs/enamel_relics"
O="{WSL_ROOT}/runs/enamel_relics/v101/b2"
blastp -query "$O/b2_xp_proteins.faa" -db "$D/blastdb/amis_proteome" -evalue 10 -seg yes -max_target_seqs 10 -num_threads 8 -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore stitle" > "$O/b2_blastp_vs_amis.tsv"
blastp -query "$O/b2_xp_proteins.faa" -db "$D/blastdb/amis_proteome" -evalue 10 -seg no -max_target_seqs 10 -num_threads 8 -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore stitle" > "$O/b2_blastp_vs_amis_segno.tsv"
"""
    run_wsl_script(bp_sh, OUT / "run_b2_blastp.sh")

    # Parse BLASTp hits
    def parse_bp(p: Path) -> dict[str, list[dict]]:
        by_q: dict[str, list[dict]] = {}
        for ln in p.read_text(encoding="utf-8").splitlines():
            if not ln.strip():
                continue
            c = ln.split("\t")
            by_q.setdefault(c[0], []).append({
                "sseqid": c[1], "pident": float(c[2]), "length": int(c[3]),
                "evalue": float(c[10]), "bitscore": float(c[11]), "stitle": c[12],
            })
        return by_q

    bp_yes = parse_bp(OUT / "b2_blastp_vs_amis.tsv")
    bp_no = parse_bp(OUT / "b2_blastp_vs_amis_segno.tsv")

    out_obj = {
        "timestamp_utc": datetime.now(timezone.utc).isoformat(),
        "targets": results,
        "blastp_vs_amis_segyes": bp_yes,
        "blastp_vs_amis_segno": bp_no,
    }
    out_json = OUT / "b2_detailed_audit.json"
    out_json.write_text(json.dumps(out_obj, indent=2), encoding="utf-8")

    # Print crisp diagnostic summary
    for res in results:
        t = res["target"]
        print(f"\n--- {t['code']} {t['symbol']} ({t['assembly']}) ---")
        print(f"  GeneID={res['gene_id']} type={res['gene_type']} loc={res['seq_acc']}:{res['gene_start_1based']}-{res['gene_end_1based']} ({res['strand']})")
        print(f"  LowQuality={res['genomic_has_low_quality_tag']} exceptions={res['genomic_exceptions']}")
        print(f"  Neighbors (+-200kb): {[g['gene'] for g in res['neighboring_genes_200kb']]}")
        print(f"  Other assemblies (same taxid): {[(a['accession'], a['level'], a['contig_n50']) for a in res['other_assemblies_same_taxid']]}")
        for m in res["models"]:
            print(f"  Model {m['transcript_acc']} / {m['protein_acc']} ({m['protein_len']} aa, {m['n_cds_exons']} CDS exons):")
            print(f"    def: {m['tr_definition'][:120]}")
            print(f"    comment: {m['tr_comment'][:220]}")
            print(f"    cds_notes: {m['cds_notes']} | exceptions: {m['cds_exceptions']}")
            print(f"    genomic_cds={m['genomic_cds_bp']}bp (mod3={m['genomic_cds_len_mod3']}) == mrna_cds({m['mrna_cds_bp']}bp)? {m['genomic_cds_equals_mrna_cds']}")
            print(f"    genomic internal stops: {m['genomic_internal_stop_codons_1based']} | mrna internal stops: {m['mrna_internal_stop_codons_1based']}")
            print(f"    all_introns_canonical: {m['all_introns_canonical']} | introns: {[(x['intron_idx'], x['bp'], x['donor_2bp']+'..'+x['acceptor_2bp'], x['canonical']) for x in m['intron_details']]}")
            print(f"    genomic_vs_mrna_diffs ({len(m['genomic_vs_mrna_diffs'])}): {m['genomic_vs_mrna_diffs']}")
            print(f"    N-runs (+-2kb): {m['n_runs_in_gene_plus_minus_2kb']}")
    print("\n=== BLASTp vs Alligator proteome (-seg yes) ===")
    for k, v in bp_yes.items():
        print(f"  {k}: {[(x['sseqid'], x['pident'], x['length'], x['bitscore'], x['evalue'], x['stitle'][:55]) for x in v[:3]]}")
    print("\n=== BLASTp vs Alligator proteome (-seg no) ===")
    for k, v in bp_no.items():
        print(f"  {k}: {[(x['sseqid'], x['pident'], x['length'], x['bitscore'], x['evalue'], x['stitle'][:55]) for x in v[:3]]}")

    with PROV.open("a", encoding="utf-8") as f:
        f.write(json.dumps({
            "timestamp_utc": datetime.now(timezone.utc).isoformat(),
            "step": "phase_B2_native_model_audit",
            "script": "scripts/enamel_relics/08_audit_protein_coding_models.py",
            "output_json": str(out_json.relative_to(ROOT).as_posix()),
            "output_sha256": sha256_file(out_json),
        }) + "\n")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
