"""Controle de composicao (critica externa, item b): o tBLASTn e cego para
proteina desordenada nas aves?

ENAM, AMBN e AMTN sao fosfoproteinas secretadas, desordenadas, ricas em
Pro/Gln — a classe que o SEG mascara e que diverge rapido. Se o metodo nao
enxerga ESSA CLASSE nas aves, o "zero" das aves no s33 seria cegueira, nao
perda. O controle e uma proteina da mesma classe que as aves MANTEM funcional:
IBSP (sialoproteina ossea, cluster SIBLING, anotada em todas as especies-base).
Ela e ancora dos intervalos do DSPP do s27 (runs/atavismo/painel/
intervalo_DSPP_<rot>.fna), entao a janela ja existe.

Consulta = IBSP do Anolis, a mesma distancia de crocodilianos, aves e
tartaruga (Archelosauria x Lepidosauria). Mesmo nulo do s33 (20 embaralhadas,
excesso sobre o percentil 95).

    uv run python scripts/atavismo/s34_composicao.py preparar
    wsl -d Ubuntu-24.04 -- bash -l /mnt/c/Users/Manec/bio-stack/scripts/atavismo/s34_alinhar.sh
    uv run python scripts/atavismo/s34_composicao.py analisar
"""

from __future__ import annotations

import csv
import importlib.util
import sys
from pathlib import Path

import numpy as np

AQUI = Path(__file__).resolve().parent
_spec = importlib.util.spec_from_file_location("s33", AQUI / "s33_indice_nulo.py")
if _spec is None or _spec.loader is None:
    raise SystemExit("nao carreguei s33_indice_nulo.py")
s33 = importlib.util.module_from_spec(_spec)
sys.modules["s33"] = s33
_spec.loader.exec_module(s33)

PAINEL = s33.RAIZ / "runs" / "atavismo" / "painel"
OUT = s33.RAIZ / "runs" / "atavismo" / "composicao"
GENE = "IBSP"
ESPECIES = {"Amis": "COM", "Ggan": "COM", "Acar": "COM", "Cpic": "sem", "Ggal": "sem", "Shab": "sem"}


def preparar() -> int:
    OUT.mkdir(parents=True, exist_ok=True)
    p = s33.s27.gene_pos(GENE, 28377)
    if not p:
        raise SystemExit("IBSP nao anotado no Anolis")
    genes = s33.s27.tabela_features(p["acc"], p["ini"], p["fim"])
    g = next((x for x in genes if x.sym.upper() == GENE and x.proteina), None)
    if not (g and g.proteina):
        raise SystemExit("IBSP do Anolis sem proteina na tabela de features")
    seq = "".join(s33.s27.proteinas({GENE: g.proteina}).split("\n", 1)[1].split())
    rng = np.random.default_rng(s33.SEMENTE)
    linhas = [f">{GENE}_Acar {g.proteina}\n{seq}\n"]
    for i in range(s33.N_EMBARALHADAS):
        linhas.append(f">{GENE}_Acar~emb{i:02d}\n{''.join(rng.permutation(list(seq)))}\n")
    (OUT / "consultas.faa").write_text("".join(linhas), encoding="utf-8", newline="\n")
    print(f"IBSP do Anolis: {g.proteina}, {len(seq)} aa, + {s33.N_EMBARALHADAS} embaralhadas")
    return 0


def analisar() -> int:
    q = f"{GENE}_Acar"
    linhas = []
    print(f"{'esp':5s} {'dente':5s} {'IBSP obs':>9s} {'nulo95':>7s} {'excesso':>8s}")
    for rot, dente in ESPECIES.items():
        tsv = OUT / f"s34_{rot}.tsv"
        if not tsv.exists():
            continue
        melhor: dict[str, float] = {}
        with tsv.open(encoding="utf-8") as f:
            for c in csv.reader(f, delimiter="\t"):
                melhor[c[0]] = max(melhor.get(c[0], 0.0), float(c[9]))
        obs = melhor.get(q, 0.0)
        n95 = float(np.percentile([melhor.get(f"{q}~emb{i:02d}", 0.0)
                                   for i in range(s33.N_EMBARALHADAS)], 95))
        print(f"{rot:5s} {dente:5s} {obs:9.1f} {n95:7.1f} {max(0.0, obs - n95):8.1f}")
        linhas.append({"especie": rot, "dente": dente, "ibsp_obs": obs, "nulo95": n95,
                       "excesso": round(max(0.0, obs - n95), 1)})
    with (OUT / "s34_ibsp.tsv").open("w", encoding="utf-8", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(linhas[0]), delimiter="\t", lineterminator="\n")
        w.writeheader()
        w.writerows(linhas)
    return 0


if __name__ == "__main__":
    modo = sys.argv[1] if len(sys.argv) > 1 else ""
    if modo not in ("preparar", "analisar"):
        raise SystemExit("uso: s34_composicao.py preparar|analisar")
    raise SystemExit(preparar() if modo == "preparar" else analisar())
