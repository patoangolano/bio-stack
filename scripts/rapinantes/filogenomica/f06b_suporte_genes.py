"""
f06b — resumo do suporte nas arvores de genes: por loco, fracao de ramos internos
com SH-aLRT >= 80 e UFBoot >= 95 (criterio do IQ-TREE para arvore de gene unico),
e a distribuicao entre locos. Tambem resume os ramos terminais anomalos do f06.

    uv run python scripts/rapinantes/filogenomica/f06b_suporte_genes.py <run_id> <ramos_anomalos.tsv> <saida.tsv>
"""

from __future__ import annotations

import csv
import re
import statistics
import sys
from collections import Counter
from pathlib import Path

RAIZ = Path(__file__).resolve().parents[3]


def main() -> None:
    run = RAIZ / "runs" / "filogenia" / sys.argv[1]
    linhas = []
    for f in sorted((run / "04_arvores_genes").glob("*.treefile")):
        rot = re.findall(r"\)([0-9.]+)/([0-9.]+)/([0-9.]+)", f.read_text(encoding="utf-8"))
        n = len(rot)
        fortes = sum(1 for sh, _ab, ufb in rot if float(sh) >= 80 and float(ufb) >= 95)
        linhas.append({"loco": f.stem, "ramos_internos": n, "fortes_sh80_ufb95": fortes,
                       "fracao": f"{fortes / n:.3f}" if n else ""})
    with Path(sys.argv[3]).open("w", encoding="utf-8", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(linhas[0]), delimiter="\t", lineterminator="\n")
        w.writeheader()
        w.writerows(linhas)
    fr = sorted(float(r["fracao"]) for r in linhas if r["fracao"])
    q = statistics.quantiles(fr, n=4)
    print(f"{len(linhas)} arvores de genes; fracao de ramos internos com SH-aLRT>=80 e UFBoot>=95: "
          f"mediana {statistics.median(fr):.3f} (IQR {q[0]:.3f}-{q[2]:.3f}; min {fr[0]:.3f}, max {fr[-1]:.3f})")
    an = list(csv.DictReader(Path(sys.argv[2]).open(encoding="utf-8"), delimiter="\t"))
    graves = [r for r in an if float(r["ramo"]) > 1.0]
    print(f"ramos terminais anomalos: {len(an)} em {len({r['loco'] for r in an})} locos; "
          f"> 1 subst./sitio: {len(graves)} em {len({r['loco'] for r in graves})} locos")
    print("por taxon (top 8):", Counter(r["taxon"] for r in an).most_common(8))


if __name__ == "__main__":
    main()
