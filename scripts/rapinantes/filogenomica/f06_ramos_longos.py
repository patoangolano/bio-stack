"""
f06 — triagem de ramos terminais anomalos nas arvores de genes (suspeita de
paralogo, contaminacao ou erro de extracao; Regra 5 do agente de filogenia).

Para cada taxon, a referencia e a mediana do proprio ramo terminal entre os
locos; um ramo e marcado quando passa de 5x essa mediana e de 0,3 substituicao
por sitio. Conta tambem por fonte (RefSeq vs miniprot), para ver se a extracao
por miniprot produz mais ramos anomalos.

    uv run python scripts/rapinantes/filogenomica/f06_ramos_longos.py <run_id> <saida.tsv>
"""

from __future__ import annotations

import csv
import io
import re
import statistics
import sys
from collections import defaultdict
from pathlib import Path

from Bio import Phylo

RAIZ = Path(__file__).resolve().parents[3]
FATOR, MIN_ABS = 5.0, 0.3


def main() -> None:
    run = RAIZ / "runs" / "filogenia" / sys.argv[1]
    saida = Path(sys.argv[2])
    fonte = {t["taxon"]: t["fonte"] for t in csv.DictReader(
        (RAIZ / "data/proc/rapinantes/filogenomica/taxons.tsv").open(encoding="utf-8"), delimiter="\t")}
    term: dict[str, dict[str, float]] = defaultdict(dict)
    for f in sorted((run / "04_arvores_genes").glob("*.treefile")):
        txt = re.sub(r"\)[^:,;()]+", ")", f.read_text(encoding="utf-8").strip())
        arv = Phylo.read(io.StringIO(txt), "newick")  # pyright: ignore[reportPrivateImportUsage]
        for t in arv.get_terminals():
            term[t.name][f.stem] = t.branch_length or 0.0
    med = {t: statistics.median(v.values()) for t, v in term.items()}
    marcados = []
    for t, v in term.items():
        for loco, bl in v.items():
            if bl > FATOR * med[t] and bl > MIN_ABS:
                marcados.append({"loco": loco, "taxon": t, "fonte": fonte.get(t, "?"), "ramo": f"{bl:.4f}",
                                 "mediana_taxon": f"{med[t]:.4f}", "razao": f"{bl / med[t]:.1f}"})
    marcados.sort(key=lambda r: -float(r["razao"]))
    saida.parent.mkdir(parents=True, exist_ok=True)
    with saida.open("w", encoding="utf-8", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=["loco", "taxon", "fonte", "ramo", "mediana_taxon", "razao"],
                           delimiter="\t", lineterminator="\n")
        w.writeheader()
        w.writerows(marcados)
    n_locos = len({loco for v in term.values() for loco in v})
    por_fonte: dict[str, list[int]] = defaultdict(lambda: [0, 0])
    for t, v in term.items():
        por_fonte[fonte.get(t, "?")][1] += len(v)
    for r in marcados:
        por_fonte[r["fonte"]][0] += 1
    print(f"{n_locos} arvores; {len(marcados)} ramos terminais marcados (> {FATOR}x a mediana do taxon e > {MIN_ABS}); "
          + "; ".join(f"{k}: {a}/{b} ({100 * a / b:.2f} %)" for k, (a, b) in por_fonte.items()))
    for r in marcados[:10]:
        print("\t".join(r.values()))


if __name__ == "__main__":
    main()
