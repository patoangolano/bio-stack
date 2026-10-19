"""
f08b — tabela dos nos-chave so com arvores de genes e wASTRAL (antes da supermatriz).

Por no-chave (biparticao) e conjunto:
  * wASTRAL do pipeline (modo hibrido, ponderado pelo aBayes): presente?, localPP;
  * quartetos nao ponderados na mesma topologia (f08a: wastral --mode 4 -C -u 2):
    q1/q2/q3 e qui-quadrado f2 x f3 (H0: as duas alternativas na mesma frequencia);
  * gCF das arvores de genes sobre a topologia do wASTRAL (f08a, IQ-TREE --gcf):
    gCF, gDF1, gDF2, gDFP, gN e qui-quadrado gDF1_N x gDF2_N;
  * arvores de genes com o no (restrito aos taxons do gene) a SH-aLRT >= 80 e
    UFBoot >= 95, e com um ramo nesse criterio que o contradiz.
No ausente do wASTRAL: as colunas de quarteto e gCF ficam vazias.

    uv run python scripts/rapinantes/filogenomica/f08b_tabela_genes_astral.py <saida.tsv> rotulo=run_id=dir_f08a [...]
"""

from __future__ import annotations

import csv
import importlib.util
import json
import sys
from pathlib import Path

from scipy.stats import chisquare

RAIZ = Path(__file__).resolve().parents[3]
_spec = importlib.util.spec_from_file_location("f08", Path(__file__).with_name("f08_nos_chave.py"))
assert _spec and _spec.loader
f08 = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(f08)


def attrs(rotulo: str) -> dict[str, str]:
    return dict(kv.split("=", 1) for kv in rotulo.strip("[]").split(";") if "=" in kv)


def p2(a: float, b: float) -> str:
    return f"{chisquare([a, b]).pvalue:.3g}" if a + b > 0 else ""


def main() -> None:
    saida = Path(sys.argv[1])
    linhas = []
    for arg in sys.argv[2:]:
        rotulo, run_id, d = arg.split("=", 2)
        run = RAIZ / "runs" / "filogenia" / run_id
        dq = Path(d)
        rot_sup = json.loads((run / "resumo.json").read_text(encoding="utf-8")).get("rotulos_suporte") if (run / "resumo.json").is_file() else None
        rot_sup = rot_sup or ["SH-aLRT", "aBayes", "UFBoot"]
        # o resumo.json grava alrt/abayes/ufboot; aqui os nomes de exibicao
        rot_sup = [{"alrt": "SH-aLRT", "abayes": "aBayes", "ufboot": "UFBoot"}.get(x, x) for x in rot_sup]
        a_taxa, b_a = f08.biparticoes(f08.ler(run / "07_arvore_especies" / "especies.tre"))
        localpp = {b[0]: b[1] for b in b_a}
        _, b_q = f08.biparticoes(f08.ler(dq / "quartetos_astral.tre"))
        quart = {b[0]: attrs(b[1]) for b in b_q}
        _, b_id = f08.biparticoes(f08.ler(dq / "gcf_astral.cf.branch"))
        id2bip = {b[1]: b[0] for b in b_id}
        gcf = {id2bip[r["ID"]]: r for r in f08.ler_stat(dq / "gcf_astral.cf.stat") if r["ID"] in id2bip}
        genes = []
        for gf in sorted((run / "04_arvores_genes").glob("*.treefile")):
            genes.append(f08.biparticoes(f08.ler(gf)))
        for nome, expr in f08.NOS:
            c = f08.conj(expr, a_taxa)
            if len(c) < 2 or len(a_taxa) - len(c) < 2:
                continue
            k = f08.canon(c, a_taxa)
            r: dict[str, object] = {"conjunto": rotulo, "run_id": run_id, "no": nome, "n_taxons_no": len(c),
                                    "no_wASTRAL": "sim" if k in localpp else "nao", "localPP": localpp.get(k, "")}
            q = quart.get(k)
            if q:
                for col in ("q1", "q2", "q3"):
                    r[col] = f"{float(q[col]):.3f}"
                r["p_q2_vs_q3"] = p2(float(q["f2"]), float(q["f3"]))
            g = gcf.get(k)
            if g:
                for col in ("gCF", "gDF1", "gDF2", "gDFP", "gN"):
                    r[col] = g[col]
                r["gCF_N/gDF1_N/gDF2_N"] = f"{g['gCF_N']}/{g['gDF1_N']}/{g['gDF2_N']}"
                r["p_gDF1_vs_gDF2"] = p2(float(g["gDF1_N"]), float(g["gDF2_N"]))
            inf = fav = con = 0
            for gt, gb in genes:
                cg = frozenset(c & gt)
                if len(cg) < 2 or len(gt) - len(cg) < 2:
                    continue
                inf += 1
                kg = f08.canon(cg, gt)
                fortes = []
                for lado, lab, _ in gb:
                    v = lab.split("/")
                    try:
                        if float(v[rot_sup.index("SH-aLRT")]) >= 80 and float(v[rot_sup.index("UFBoot")]) >= 95:
                            fortes.append(lado)
                    except (ValueError, IndexError):
                        continue
                if kg in fortes:
                    fav += 1
                elif any(not f08.compativel(kg, f, gt) for f in fortes):
                    con += 1
            r.update(genes_informativos=inf, genes_fortes_a_favor=fav, genes_fortes_contra=con)
            linhas.append(r)
    campos: list[str] = []
    for r in linhas:
        campos += [k for k in r if k not in campos]
    saida.parent.mkdir(parents=True, exist_ok=True)
    with saida.open("w", encoding="utf-8", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=campos, delimiter="\t", lineterminator="\n", restval="")
        w.writeheader()
        w.writerows(linhas)
    print(f"{len(linhas)} linhas -> {saida}")


if __name__ == "__main__":
    main()
