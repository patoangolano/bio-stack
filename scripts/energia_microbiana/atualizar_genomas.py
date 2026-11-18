"""Inventário direcionado de montagens NCBI; não infere eficiência de genomas.

Executar no WSL com shell de login. A lista é um painel de candidatos, não um
censo de toda a biodiversidade. Retém até dez resultados por consulta e informa
o total anunciado pelo NCBI. Não baixa sequências nem altera dados brutos.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import subprocess
from collections import Counter
from datetime import UTC, datetime
from pathlib import Path


def now() -> str:
    """Return a UTC ISO timestamp."""
    return datetime.now(UTC).isoformat()


def flatten(candidate: dict, report: dict) -> dict:
    """Extract metadata without equating an annotation with a phenotype."""
    organism = report.get("organism", {})
    info = report.get("assembly_info", {})
    stats = report.get("assembly_stats", {})
    annotation = report.get("annotation_info", {})
    biosample = info.get("biosample", {})
    strain = organism.get("infraspecific_names", {}).get("strain") or biosample.get("strain", "")
    return {
        "id": candidate["id"],
        "grupo": candidate["grupo"],
        "consulta_taxon": candidate["taxon"],
        "filtro_search": candidate["search"],
        "accession": report.get("accession", ""),
        "organism": organism.get("organism_name", ""),
        "tax_id": organism.get("tax_id", ""),
        "strain": strain,
        "level": info.get("assembly_level", ""),
        "assembly_status": info.get("assembly_status", ""),
        "reference_category": info.get("refseq_category", ""),
        "size_bp": stats.get("total_sequence_length", ""),
        "gc_percent": stats.get("gc_percent", ""),
        "protein_coding": annotation.get("stats", {})
        .get("gene_counts", {})
        .get("protein_coding", ""),
        "annotation_date": annotation.get("release_date", ""),
        "sequencing_tech": info.get("sequencing_tech", ""),
        "assembly_type": info.get("assembly_type", ""),
        "biosample_models": "; ".join(biosample.get("models", [])),
        "url": "https://www.ncbi.nlm.nih.gov/datasets/genome/" + report.get("accession", "") + "/",
    }


def rank(row: dict) -> tuple:
    """Rank only the returned subset; do not claim a global best assembly."""
    levels = {"Complete Genome": 0, "Chromosome": 1, "Scaffold": 2, "Contig": 3}
    return (
        row["assembly_status"] != "current",
        levels.get(row["level"], 9),
        not row["accession"].startswith("GCF_"),
        row["reference_category"] != "reference genome",
        row["accession"],
    )


def write_tsv(path: Path, rows: list[dict], fields: list[str]) -> None:
    """Write a UTF-8 table with a header, including an empty result set."""
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def main() -> None:
    """Query the panel with a bounded CLI contract and preserve raw responses."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--candidates", type=Path, required=True)
    parser.add_argument("--out", type=Path, required=True)
    parser.add_argument("--only", nargs="*")
    args = parser.parse_args()
    args.out.mkdir(parents=True, exist_ok=True)
    raw_dir = args.out / "ncbi_raw"
    raw_dir.mkdir(exist_ok=True)
    with args.candidates.open(encoding="utf-8") as handle:
        candidates = list(csv.DictReader(handle, delimiter="\t"))
    all_rows: list[dict] = []
    selected: list[dict] = []
    failures: list[dict] = []
    ledger: list[dict] = []
    cli_version = subprocess.run(
        ["datasets", "version"], capture_output=True, text=True, check=False
    ).stdout.strip()
    for candidate in candidates:
        cache = raw_dir / f"{candidate['id']}.json"
        if args.only and candidate["id"] not in args.only and not cache.exists():
            continue
        cmd = [
            "datasets",
            "summary",
            "genome",
            "taxon",
            candidate["taxon"],
            "--limit",
            "10",
        ]
        if candidate["search"]:
            cmd.extend(["--search", candidate["search"]])
        if cache.exists():
            envelope = json.loads(cache.read_text(encoding="utf-8"))
        else:
            try:
                proc = subprocess.run(cmd, capture_output=True, text=True, timeout=100, check=False)
                data = json.loads(proc.stdout) if proc.stdout.strip() else {}
                envelope = {
                    "retrieved_at": now(),
                    "command": cmd,
                    "returncode": proc.returncode,
                    "stderr": proc.stderr,
                    "data": data,
                }
            except (subprocess.TimeoutExpired, json.JSONDecodeError) as exc:
                envelope = {
                    "retrieved_at": now(),
                    "command": cmd,
                    "returncode": -1,
                    "stderr": type(exc).__name__,
                    "data": {},
                }
            cache.write_text(json.dumps(envelope, ensure_ascii=False, indent=2), encoding="utf-8")
        data = envelope.get("data", {})
        reports = data.get("reports", [])
        rows = [flatten(candidate, report) for report in reports]
        for row in rows:
            row["retrieved_at"] = envelope["retrieved_at"]
        all_rows.extend(rows)
        current_rows = [row for row in rows if row["assembly_status"] == "current"]
        if current_rows:
            chosen = dict(sorted(current_rows, key=rank)[0])
            chosen["interpretacao"] = (
                "Montagem representativa entre os resultados devolvidos; "
                "conferir linhagem no artigo antes de atribuir fenótipos."
            )
            selected.append(chosen)
        else:
            failures.append(
                {
                    "id": candidate["id"],
                    "taxon": candidate["taxon"],
                    "search": candidate["search"],
                    "status": (
                        "consulta_falhou"
                        if envelope["returncode"] != 0
                        else "sem_montagem_atual_na_consulta"
                    ),
                    "stderr": envelope.get("stderr", "")[:1500],
                }
            )
        ledger.append(
            {
                "id": candidate["id"],
                "retrieved_at": envelope["retrieved_at"],
                "command": envelope["command"],
                "returncode": envelope["returncode"],
                "total_count": data.get("total_count"),
                "retrieved": len(reports),
                "current_in_subset": len(current_rows),
                "cache_sha256": hashlib.sha256(cache.read_bytes()).hexdigest(),
                "exhaustive": False,
            }
        )
        print(
            f"{candidate['id']}: {len(reports)} registros; {len(current_rows)} atuais",
            flush=True,
        )
    if all_rows:
        fields = list(all_rows[0])
        write_tsv(args.out / "montagens_retornadas.tsv", all_rows, fields)
        write_tsv(
            args.out / "genomas_representativos.tsv",
            selected,
            fields + ["interpretacao"],
        )
    summary = {
        "generated_at": now(),
        "datasets_version": cli_version,
        "panel_entries": len(candidates),
        "queries_in_run": len(ledger),
        "returned_assemblies": len(all_rows),
        "representatives_current": len(selected),
        "representative_levels": dict(Counter(row["level"] for row in selected)),
        "failed_or_unresolved": len(failures),
        "failures": failures,
        "ledger": ledger,
        "limitations": [
            "Painel direcionado, não censo de todas as espécies ou montagens.",
            "Até dez registros por consulta, sem paginação exaustiva.",
            "Nomes e campos de linhagem podem mudar; resultados retêm ambos.",
            "Montagem completa não implica cultura pura, expressão ou eficiência.",
            "Contagem CDS não é capacidade metabólica experimental.",
        ],
    }
    (args.out / "inventario_ncbi.json").write_text(
        json.dumps(summary, indent=2, ensure_ascii=False), encoding="utf-8"
    )
    print(json.dumps({k: v for k, v in summary.items() if k != "ledger"}, indent=2))


if __name__ == "__main__":
    main()
