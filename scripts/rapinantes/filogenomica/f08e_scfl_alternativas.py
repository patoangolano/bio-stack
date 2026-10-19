"""
f08e — sCF/sDF do ramo que cada arvore restrita forca (saida do f08d), lido como biparticao, como no f08.

    uv run python scripts/rapinantes/filogenomica/f08e_scfl_alternativas.py <dir das hipoteses> <saida.tsv>
"""

from __future__ import annotations

import csv
import importlib.util
import sys
from pathlib import Path

_spec = importlib.util.spec_from_file_location("f08", Path(__file__).with_name("f08_nos_chave.py"))
assert _spec and _spec.loader
f08 = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(f08)

RAMOS: list[tuple[str, str, str]] = [
    ("H2", "Strigiformes + Accipitrimorphae (Hieraves)", "STRIG+ACCMORPH"),
    ("H1b", "Accipitrimorphae + Coraciimorphae", "ACCMORPH+CORAC"),
    ("H7b", "Cariamiformes + Falconiformes", "CARI+FALC"),
]


def main() -> None:
    hip, saida = Path(sys.argv[1]), Path(sys.argv[2])
    linhas = []
    for h, nome, expr in RAMOS:
        taxa, bs = f08.biparticoes(f08.ler(hip / f"{h}.scfl.cf.tree"))
        stat = f08.casar(f08.ler_stat(hip / f"{h}.scfl.cf.stat"), bs)
        r = stat.get(f08.canon(f08.conj(expr, taxa), taxa))
        if r is None:
            sys.exit(f"{h}: ramo {expr} nao encontrado na arvore restrita")
        linhas.append(
            {
                "hipotese": h,
                "ramo_forcado": nome,
                "sCF": r["sCF"],
                "sDF1": r["sDF1"],
                "sDF2": r["sDF2"],
                "sN": r["sN"],
                "p_sDF1_vs_sDF2": f08.pval(r["sDF1_N"], r["sDF2_N"]),
            }
        )
    with saida.open("w", encoding="utf-8", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(linhas[0]), delimiter="\t", lineterminator="\n")
        w.writeheader()
        w.writerows(linhas)
    for r in linhas:
        print("\t".join(str(v) for v in r.values()))


if __name__ == "__main__":
    main()
