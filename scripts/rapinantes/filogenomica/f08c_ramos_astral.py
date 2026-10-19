"""
f08c — todos os ramos internos do wASTRAL (nao so os nos-chave pre-definidos),
com localPP, q1/q2/q3, gCF/gDF1/gDF2/gDFP e qui-quadrados. O clado de cada ramo
e descrito pelos grupos de f07 (lado sem a raiz Gallus/Anas).

    uv run python scripts/rapinantes/filogenomica/f08c_ramos_astral.py <run_id> <dir_f08a> <saida.tsv>
"""

from __future__ import annotations

import csv
import importlib.util
import sys
from pathlib import Path

from scipy.stats import chisquare

RAIZ = Path(__file__).resolve().parents[3]
_spec = importlib.util.spec_from_file_location("f08", Path(__file__).with_name("f08_nos_chave.py"))
assert _spec and _spec.loader
f08 = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(f08)

ROTULOS = ["CATH", "SAG", "PAN", "ACCIPITRIDAE", "STRIG", "CORAC", "FALC", "PSIT", "PASS", "CARI",
           "STRISORES"]


def descreve(lado: frozenset[str]) -> str:
    resto = set(lado)
    partes = []
    for g in sorted(ROTULOS, key=lambda g: -len(f08.G[g])):
        s = f08.G[g]
        if s and s <= resto and len(s) > 1:
            partes.append(g)
            resto -= s
    partes += sorted(resto)
    return " + ".join(partes)


def main() -> None:
    run = RAIZ / "runs" / "filogenia" / sys.argv[1]
    dq = Path(sys.argv[2])
    taxa, b_a = f08.biparticoes(f08.ler(run / "07_arvore_especies" / "especies.tre"))
    raiz = f08.G["RAIZ"] & taxa
    pp = {b[0]: b[1] for b in b_a}
    _, b_q = f08.biparticoes(f08.ler(dq / "quartetos_astral.tre"))
    q = {b[0]: dict(kv.split("=", 1) for kv in b[1].strip("[]").split(";") if "=" in kv) for b in b_q}
    _, b_id = f08.biparticoes(f08.ler(dq / "gcf_astral.cf.branch"))
    id2 = {b[1]: b[0] for b in b_id}
    gcf = {id2[r["ID"]]: r for r in f08.ler_stat(dq / "gcf_astral.cf.stat") if r["ID"] in id2}
    linhas = []
    for k in pp:
        lado = frozenset(taxa - k) if k & raiz else k
        qq, g = q.get(k, {}), gcf.get(k, {})
        f2, f3 = float(qq.get("f2", 0)), float(qq.get("f3", 0))
        d1, d2 = float(g.get("gDF1_N", 0) or 0), float(g.get("gDF2_N", 0) or 0)
        linhas.append({
            "n": len(lado), "clado": descreve(lado), "localPP": pp[k],
            "q1": qq.get("q1", ""), "q2": qq.get("q2", ""), "q3": qq.get("q3", ""),
            "p_q2_vs_q3": f"{chisquare([f2, f3]).pvalue:.3g}" if f2 + f3 else "",
            "gCF": g.get("gCF", ""), "gDF1": g.get("gDF1", ""), "gDF2": g.get("gDF2", ""), "gDFP": g.get("gDFP", ""),
            "p_gDF1_vs_gDF2": f"{chisquare([d1, d2]).pvalue:.3g}" if d1 + d2 else "",
        })
    linhas.sort(key=lambda r: -int(r["n"]))
    with Path(sys.argv[3]).open("w", encoding="utf-8", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(linhas[0]), delimiter="\t", lineterminator="\n")
        w.writeheader()
        w.writerows(linhas)
    print(f"{len(linhas)} ramos -> {sys.argv[3]}")


if __name__ == "__main__":
    main()
