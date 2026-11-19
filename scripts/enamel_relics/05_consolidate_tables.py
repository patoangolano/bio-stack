"""Phase 2.5 — Consolidate all verified Phase 2 outputs into publication-ready
supplementary tables (runs/enamel_relics/supplementary/) and SHA-256 manifest.

Usage:
    uv run python scripts/enamel_relics/05_consolidate_tables.py
"""

from __future__ import annotations

import csv
import hashlib
import json
import shutil
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / "runs" / "enamel_relics"
TABLES = OUT / "tables"
SUPP = OUT / "supplementary"
LEGACY = ROOT / "runs" / "atavismo"

TOOTH_GENES = ["ENAM", "AMBN", "AMTN", "AMELX", "MMP20", "DSPP", "ODAM", "ODAPH", "ACP4"]
CTRL_GENES = ["NRG1", "LPL", "PSD3", "IBSP", "SPP1", "SPARCL1", "CDKL2", "THAP6", "RCHY1"]


def sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def write_tsv(path: Path, rows: list[dict[str, Any]]) -> None:
    if not rows:
        return
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0]), delimiter="\t", lineterminator="\n")
        w.writeheader()
        w.writerows(rows)


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8") as f:
        return list(csv.DictReader(f, delimiter="\t"))


def main() -> int:
    SUPP.mkdir(parents=True, exist_ok=True)
    meta = json.loads((OUT / "audit_metadata.json").read_text(encoding="utf-8"))

    # --- Table S1: Assemblies & RefSeq annotations across 18 study species + GRCh38 ---
    s1_rows: list[dict[str, Any]] = []
    for code, asm in meta["assemblies"].items():
        ann = meta["annotations"][code]
        tooth_status = {}
        n_tooth = 0
        for g in TOOTH_GENES:
            rec = ann.get(g)
            if rec:
                n_tooth += 1
                tooth_status[f"RefSeq_{g}"] = f"{rec['type']}({rec['seq']})"
            else:
                tooth_status[f"RefSeq_{g}"] = "absent"
        ctrl_status = {f"RefSeq_{g}": (ann[g]["seq"] if ann.get(g) else "absent") for g in CTRL_GENES}
        note = ""
        if code in ("Cpic", "Cmyd", "Scam", "Arow"):
            note = "Automated RefSeq protein-coding annotation on pseudogenized enamel locus (false positive at Level 1)"
        elif code == "Ggal":
            note = "Ensembl Release 116 GTF also audited: 0/9 tooth genes annotated; 0 overlaps with psi-ENAM or psi-AMBN"
        elif code == "Amis":
            note = "ODAPH present as uncharacterized LOC106740521 (XP_019349100.2, NC_081825.1:164672404-164676187)"
        s1_rows.append({
            "species_code": code,
            "scientific_name": asm["scientific_name"],
            "common_name": asm["common_name"],
            "taxid": asm["taxid"],
            "clade": asm["clade"],
            "subclade": asm["subclade"],
            "dentition": asm["dentition"],
            "assembly_accession": asm["assembly_accession"],
            "assembly_name": asm["assembly_name"],
            "assembly_level": asm["assembly_level"],
            "contig_n50_bp": asm["contig_n50"],
            "scaffold_n50_bp": asm["scaffold_n50"],
            "annotation_name": asm["annotation_name"],
            "annotation_release_date": asm["annotation_release_date"],
            "n_tooth_genes_annotated_of_9": n_tooth,
            **tooth_status,
            **ctrl_status,
            "audit_note": note,
        })
    write_tsv(SUPP / "Table_S1_assemblies_and_annotations.tsv", s1_rows)

    # --- Table S2: Synteny intervals (18 species) + Whole-proteome reciprocal BLASTx ---
    shutil.copy2(OUT / "enamel_intervals.tsv", SUPP / "Table_S2a_synteny_intervals_18sp.tsv")
    shutil.copy2(TABLES / "whole_proteome_reciprocal.tsv", SUPP / "Table_S2b_whole_proteome_reciprocal_blastx.tsv")

    # --- Table S3: Compositional / SEG / Null-corrected indices ---
    shutil.copy2(TABLES / "query_composition_seg.tsv", SUPP / "Table_S3a_query_composition_seg.tsv")
    shutil.copy2(TABLES / "retention_indices_18sp.tsv", SUPP / "Table_S3b_retention_indices_18sp.tsv")
    shutil.copy2(TABLES / "tblastn_per_gene_seg_comparison.tsv", SUPP / "Table_S3c_per_gene_seg_null_18sp.tsv")

    # --- Table S4: miniprot Level 3 coding lesions & GT-AG splice-site audit ---
    s4_rows = read_tsv(TABLES / "miniprot_level3_splice_audit.tsv")
    for r in s4_rows:
        if r["query_gene"] == "ENAM_Amis" and r["species_code"] == "Cpic":
            r["qc_note"] = "Spurious 58.4-kb upstream exon jump (NC_088796.1:129180576-129180957, reciprocal BLASTx = protein SON, e=0.010) filtered by 25-kb max-intron rule"
        else:
            r["qc_note"] = ""
    write_tsv(SUPP / "Table_S4_miniprot_level3_and_splice_audit.tsv", s4_rows)

    # --- Table S5: Level 4 shared mutations + splice borders ---
    shutil.copy2(TABLES / "level4_shared_mutations.tsv", SUPP / "Table_S5a_level4_shared_mutations.tsv")
    shutil.copy2(TABLES / "level4_splice_borders.tsv", SUPP / "Table_S5b_level4_splice_borders.tsv")

    # --- Table S6: Auxiliary verified panels (ODAPH/ODAM, multiz100way, recombination, Psi-AMEL) ---
    # All strings normalized to clean academic English (zero legacy Portuguese).
    aux_rows: list[dict[str, Any]] = []

    def _en_odaph_odam(r: dict[str, str]) -> tuple[str, str, str, str]:
        st = r["estado"]
        nv = r["nivel"]
        coord = r["coordenada"].replace("(vao)", "(gap)").replace("sem vao", "0-bp gap (CDKL2-SCLT1 adjacent)")
        ctrl = r["controle"]
        # Translate status
        st = (
            st.replace("nao achado no intervalo NRG1..LPL (melhor tblastn ", "Not detected in NRG1-LPL interval (best tBLASTn ")
            .replace("; miniprot cobre ", "; miniprot covers ")
            .replace("% da consulta)", "% of query)")
            .replace("ORF sem lesão no trecho alinhado", "Intact ORF across aligned span")
            .replace("nao achado no vao (nenhum HSP); vao de ", "Not detected in syntenic gap (0 HSPs); gap length ")
            .replace("nao achado no vao (melhor ", "Not detected in syntenic gap (best ")
            .replace("nao achado no vao", "Not detected in syntenic gap")
            .replace("); vao de ", "); gap length ")
            .replace(" do gene no jacare", " of Alligator ODAPH gene span")
            .replace("nao testado", "Not tested (unannotated flanking anchors in GavGan_comp1)")
            .replace("presente: tblastn ", "Present: tBLASTn ")
            .replace("; integridade do ORF nao medida (miniprot < 50% da consulta)", "; full-ORF integrity not assessed (miniprot coverage < 50%)")
        )
        # Translate level
        nv_map = {
            "1 (negativo sem controle positivo em Sauropsida)": "Level 1 (uninformative negative; no sauropsid positive control)",
            "presente": "Present (toothed positive control)",
            "2": "Level 2 (synteny-delimited absence / gap collapse)",
            "1 (sem anotacao)": "Level 1 (flanking anchors unannotated)",
        }
        nv = nv_map.get(nv, nv)
        # Translate control note
        ctrl = (
            ctrl.replace("ancoras E<=1e-3: ", "Positive-control anchors (E<=1e-3): ")
            .replace("nenhuma", "none")
            .replace("ODAM nao localizado em nenhuma especie COM dente; sensibilidade no ODAM do Xenopus:", "ODAM not recovered in any toothed sauropsid; Xenopus positive-control sensitivity:")
            .replace("ancoras: ", "Flanking anchors: ")
            .replace("mesmas consultas acham ODAPH em ", "same queries recover ODAPH in ")
        )
        if r["gene"] == "ODAPH" and r["especie"] == "Psin":
            ctrl += "; gap coordinates from prior Pelodiscus RefSeq assembly GCF_049634645.1 (NC_134715.1)"
        return st, nv, coord, ctrl

    # 6a. ODAPH/ODAM
    for r in read_tsv(LEGACY / "odam_odaph" / "s32_odam_odaph.tsv"):
        st_en, nv_en, coord_en, ctrl_en = _en_odaph_odam(r)
        aux_rows.append({
            "module": "ODAPH_ODAM",
            "gene": r["gene"],
            "species_code": r["especie"],
            "scientific_name": r["nome"],
            "clade": r["grupo"],
            "dentition": "toothed" if r["dente"] == "COM" else "edentulous",
            "metric_summary": st_en,
            "evidence_level": nv_en,
            "coordinates": coord_en,
            "control_detail": ctrl_en,
        })
    # 6b. multiz100way
    multiz_sp_map = {
        "allMis1": ("Alligator mississippiensis", "Crocodylia", "toothed"),
        "anoCar2": ("Anolis carolinensis", "Squamata", "toothed"),
        "cheMyd1": ("Chelonia mydas", "Testudines", "edentulous"),
        "chrPic2": ("Chrysemys picta", "Testudines", "edentulous"),
        "amaVit1": ("Amazona vittata", "Aves", "edentulous"),
        "falChe1": ("Falco cherrug", "Aves", "edentulous"),
        "falPer1": ("Falco peregrinus", "Aves", "edentulous"),
        "ficAlb2": ("Ficedula albicollis", "Aves", "edentulous"),
        "zonAlb1": ("Zonotrichia albicollis", "Aves", "edentulous"),
    }
    for r in read_tsv(LEGACY / "multiz" / "s36_multiz_resumo.tsv"):
        sci, clade_en, dent_en = multiz_sp_map.get(r["montagem"], (r["montagem"], r["grupo"], "edentulous"))
        aux_rows.append({
            "module": "multiz100way_hg38",
            "gene": r["gene"],
            "species_code": r["montagem"],
            "scientific_name": sci,
            "clade": clade_en,
            "dentition": dent_en,
            "metric_summary": f"CDS_aligned_fraction={r['cobertura']}; raw_multiz_frameshifts={r['frameshifts']} (discarded by positive control)",
            "evidence_level": "Independent WGA coverage only",
            "coordinates": "hg38 chr4 (29 coding exons)",
            "control_detail": "Functional allMis1 ENAM yields 54 raw frameshifts against hg38; restricted to coverage-only support",
        })
    # 6c. Chicken chrZ recombination (Elferink 2010, galGal3)
    for r in read_tsv(LEGACY / "pendencias" / "s37_reliquia_galgal3.tsv"):
        aux_rows.append({
            "module": "chrZ_recombination_Elferink2010",
            "gene": f"psi-{r['reliquia']}",
            "species_code": "Ggal",
            "scientific_name": "Gallus gallus",
            "clade": "Aves",
            "dentition": "edentulous",
            "metric_summary": f"bin={r['bin']}; male_rate={r['macho_cmmb']} cM/Mb (pop1={r['pop1_macho_cmmb']}, pop2={r['pop2_macho_cmmb']}); dist_to_bin1783_edge={r['dist_borda_kb']} kb",
            "evidence_level": "Linkage map (32-mer anchored liftOver)",
            "coordinates": f"GRCg7b={r['grcg7b']} -> galGal3 chrZ:{r['centro_galgal3']} ({r['concordantes']}/{r['ancoras_32mer']} 32-mers)",
            "control_detail": "NRG1-LPL window weighted male rate=2.27 cM/Mb (54th percentile of 122 chrZ bins; chrZ mean=3.00, autosomes=2.97)",
        })
    # 6d. Psi-AMEL blind-spot calibration (GRCg7b chr1)
    aux_rows.append({
        "module": "PsiAMEL_blindspot_calibration",
        "gene": "psi-AMEL",
        "species_code": "Ggal",
        "scientific_name": "Gallus gallus",
        "clade": "Aves",
        "dentition": "edentulous",
        "metric_summary": "BLASTn vs EU340348: 100.0% identity across 4 exon blocks (513 bp); Alligator AMELX tBLASTn: 0 HSPs (-seg yes), max 25.8 bits (-seg no, noise floor)",
        "evidence_level": "Empirical blind-spot calibration",
        "coordinates": "NC_052532.1:124410142-124413840(-)",
        "control_detail": "Flanking anchors MSL3 (121 bits) and MID1 (455 bits) detected in same interval",
    })
    write_tsv(SUPP / "Table_S6_auxiliary_odaph_multiz_recomb_psiamel.tsv", aux_rows)

    # --- SHA-256 Manifest ---
    manifest_rows = []
    for p in sorted(SUPP.glob("*.tsv")):
        if p.name == "MANIFEST_SHA256.tsv":
            continue
        manifest_rows.append({
            "file": p.name,
            "bytes": p.stat().st_size,
            "sha256": sha256_file(p),
            "generated_utc": datetime.now(timezone.utc).isoformat(),
        })
    write_tsv(SUPP / "MANIFEST_SHA256.tsv", manifest_rows)
    print(f"Consolidated {len(manifest_rows)} supplementary TSVs into {SUPP}:")
    for r in manifest_rows:
        print(f"  {r['file']:52s} {r['bytes']:>8,} bytes  {r['sha256'][:16]}...")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
