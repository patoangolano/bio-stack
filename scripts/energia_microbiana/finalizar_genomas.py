"""Consolida o painel e verifica acessos de três linhagens históricas.

Todas as respostas originais são preservadas. Substituição por outra linhagem
é registrada como correspondência de espécie, não de fenótipo experimental.
"""

from __future__ import annotations

import csv
import json
import subprocess
from collections import Counter
from datetime import UTC, datetime
from pathlib import Path

from atualizar_genomas import flatten, rank, write_tsv

ROOT = Path("runs/energia_microbiana_2026-11-17/genomas")
PINS = {
    "rcaps": "GCF_000021865.1",
    "hcrun": "GCF_000012605.1",
    "mbark": "GCF_000195895.1",
}


def read_table(path: Path) -> list[dict]:
    """Read a generated TSV."""
    with path.open(encoding="utf-8") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def main() -> None:
    """Retain actual organism names even when taxonomy has changed."""
    candidates = read_table(Path("scripts/energia_microbiana/candidatos.tsv"))
    by_id = {row["id"]: row for row in candidates}
    rows = read_table(ROOT / "montagens_retornadas.tsv")
    fallback = read_table(ROOT / "resolucao_nomes/montagens_retornadas.tsv")
    available: dict[str, list[dict]] = {}
    for row in rows + fallback:
        available.setdefault(row["id"], []).append(row)
    cache = ROOT / "linhagens_historicas.json"
    if cache.exists():
        envelope = json.loads(cache.read_text(encoding="utf-8"))
    else:
        command = [
            "datasets",
            "summary",
            "genome",
            "accession",
            *PINS.values(),
        ]
        proc = subprocess.run(command, capture_output=True, text=True, timeout=120, check=True)
        envelope = {
            "command": command,
            "retrieved_at": datetime.now(UTC).isoformat(),
            "data": json.loads(proc.stdout),
        }
        cache.write_text(json.dumps(envelope, ensure_ascii=False, indent=2), encoding="utf-8")
    pinned = {}
    for report in envelope["data"].get("reports", []):
        key = next(
            (key for key, accession in PINS.items() if accession == report["accession"]),
            None,
        )
        if key:
            row = flatten(by_id[key], report)
            row["retrieved_at"] = envelope["retrieved_at"]
            if row["assembly_status"] == "current":
                pinned[key] = row
    chosen = []
    for candidate in candidates:
        key = candidate["id"]
        current = [row for row in available.get(key, []) if row["assembly_status"] == "current"]
        if key in pinned:
            row = dict(pinned[key])
            row["resolucao"] = "Acesso histórico de linhagem consultado diretamente."
        elif current:
            row = dict(sorted(current, key=rank)[0])
            row["resolucao"] = (
                "Consulta com filtro original."
                if row["filtro_search"] or not candidate["search"]
                else "Filtro de linhagem removido: representante da espécie, não da linhagem pedida."
            )
        else:
            raise RuntimeError(f"Sem representante atual para {key}.")
        row["taxon_pedido"] = candidate["taxon"]
        row["linhagem_filtro_pedido"] = candidate["search"]
        if key == "hvolc":
            row["grupo"] = "Haloarqueia heterotrófica — controle de extrapolação"
        chosen.append(row)
    fields = list(chosen[0])
    write_tsv(ROOT / "painel_final.tsv", chosen, fields)
    summary = {
        "generated_at": datetime.now(UTC).isoformat(),
        "entries": len(chosen),
        "distinct_accessions": len({row["accession"] for row in chosen}),
        "levels": dict(Counter(row["level"] for row in chosen)),
        "pinned_results": {key: row["organism"] for key, row in pinned.items()},
        "strain_filter_relaxed": [
            row["id"] for row in chosen if row["resolucao"].startswith("Filtro de linhagem")
        ],
        "interpretation": (
            "Existência de sequenciamento no painel; funções da literatura "
            "só podem ser atribuídas às linhagens efetivamente estudadas."
        ),
    }
    (ROOT / "painel_final_resumo.json").write_text(
        json.dumps(summary, ensure_ascii=False, indent=2), encoding="utf-8"
    )
    lines = [
        "# Inventário de genomas",
        "",
        "Montagens consultadas no NCBI; não são uma medição de eficiência.",
        "CG = Complete Genome; Chr = Chromosome; Sc = Scaffold; Ct = Contig.",
        "Resultados de espécie não transferem fenótipos entre linhagens.",
        "",
        "| # | Nome devolvido pelo NCBI | Linhagem | Acesso | Nível | Mb |",
        "|---|---|---|---|---|---|",
    ]
    abbreviations = {
        "Complete Genome": "CG",
        "Chromosome": "Chr",
        "Scaffold": "Sc",
        "Contig": "Ct",
    }
    for index, row in enumerate(chosen, 1):
        size = float(row["size_bp"]) / 1e6
        lines.append(
            f"| {index} | {row['organism']} | {row['strain'] or 'não informada'} "
            f"| [{row['accession']}]({row['url']}) "
            f"| {abbreviations[row['level']]} | {size:.3f} |"
        )
    (ROOT / "INVENTARIO.md").write_text("\n".join(lines) + "\n", encoding="utf-8")
    print(json.dumps(summary, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
