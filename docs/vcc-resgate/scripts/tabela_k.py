"""Os seis membros por k, mais as referencias, todos do MESMO bundle."""

import csv
import os

W = os.environ["HOME"] + "/vcc-avaliador/tres"


def ler(caminho: str) -> dict[str, float | None]:
    with open(caminho, encoding="utf-8") as fh:
        return {r["metric"]: (float(r["from_replicate"])
                              if r["from_replicate"] not in ("", "nan") else None)
                for r in csv.DictReader(fh)}


def bruto(caminho: str) -> dict[str, float]:
    with open(caminho, encoding="utf-8") as fh:
        for r in csv.DictReader(fh):
            if r["statistic"] == "mean":
                return {k: float(v) for k, v in r.items()
                        if k != "statistic" and v not in ("", "nan")}
    return {}


ARMS = [
    ("k=1 (v4)", f"{W}/scored_nossa_v4.csv", f"{W}/run_nossa_v4/agg_results.csv"),
    ("k=1.5", f"{W}/scored_kd1.5.csv", f"{W}/run_kd1.5/agg_results.csv"),
    ("k=2", f"{W}/scored_kd2.csv", f"{W}/run_kd2/agg_results.csv"),
    ("k=3", f"{W}/scored_kd3.csv", f"{W}/run_kd3/agg_results.csv"),
    ("k=4", f"{W}/scored_kd4.csv", f"{W}/run_kd4/agg_results.csv"),
    ("k=5", f"{W}/scored_kd5.csv", f"{W}/run_kd5/agg_results.csv"),
    ("k=6", f"{W}/scored_kd6.csv", f"{W}/run_kd6/agg_results.csv"),
    ("k=8", f"{W}/scored_kd8.csv", f"{W}/run_kd8/agg_results.csv"),
]
REFS = [
    ("controle", f"{W}/scored_controle_unico.csv"),
    ("oraculo", f"{W}/scored_oraculo_unico.csv"),
]
MEMBROS = [
    ("de_wilcoxon_direction_fidelity_yield_raw", "fidelity"),
    ("de_wilcoxon_sig_jaccard", "jaccard"),
    ("de_wilcoxon_direction_reach_raw", "reach"),
    ("de_wilcoxon_lfc_nmae", "lfc_nmae"),
    ("expr_mse_unbiased_capped_norm", "expr_mse"),
    ("pds_cosine", "pds_cosine"),
    ("avg_score", "AVG_SCORE"),
]

presentes = [(r, s, a) for r, s, a in ARMS if os.path.exists(s)]
dados = {r: ler(s) for r, s, _ in presentes}
refs = {r: ler(s) for r, s in REFS if os.path.exists(s)}

larg = 9
print("ESCALADO (from_replicate), mesmo bundle proxy-D-orion-k40, 40 alvos do contexto D")
print(f"{'membro':12s}" + "".join(f"{r:>{larg}s}" for r, _, _ in presentes)
      + "".join(f"{r:>{larg}s}" for r in refs))
for chave, rot in MEMBROS:
    linha = f"{rot:12s}"
    for r, _, _ in presentes:
        v = dados[r].get(chave)
        linha += f"{v:{larg}.4f}" if v is not None else f"{'—':>{larg}s}"
    for r in refs:
        v = refs[r].get(chave)
        linha += f"{v:{larg}.4f}" if v is not None else f"{'—':>{larg}s}"
    print(linha)

print()
print("BRUTO (statistic=mean), antes do reescalonamento")
brutos = {r: bruto(a) for r, _, a in presentes}
print(f"{'membro':12s}" + "".join(f"{r:>{larg}s}" for r, _, _ in presentes))
for chave, rot in MEMBROS[:-1]:
    linha = f"{rot:12s}"
    for r, _, _ in presentes:
        v = brutos[r].get(chave)
        linha += f"{v:{larg}.4f}" if v is not None else f"{'—':>{larg}s}"
    print(linha)

print()
if "avg_score" in next(iter(dados.values()), {}):
    piso = refs.get("controle", {}).get("avg_score")
    teto = refs.get("oraculo", {}).get("avg_score")
    base = dados.get("k=1 (v4)", {}).get("avg_score")
    for r, _, _ in presentes:
        v = dados[r]["avg_score"]
        if v is None:
            continue
        extra = ""
        if piso is not None and teto is not None:
            extra += f"  {(v - piso) / (teto - piso):6.1%} do intervalo"
        if base is not None:
            extra += f"  delta vs k=1: {v - base:+.4f}"
        print(f"  {r:10s} avg_score={v:+.4f}{extra}")
