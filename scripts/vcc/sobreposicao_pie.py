"""
Quanto dos 300 alvos do VCC 2026 (contextos D/E/F) ja tem efeito medido em
dado publico de Perturb-seq, e em quais linhagens.

O desafio de 2026 nao tem conjunto de treino: o efeito de cada perturbacao em
D/E/F tem de ser transferido de outro contexto. Esta conta decide a
estrategia — alvo com efeito medido em alguma linhagem pode ser transferido;
alvo sem medida nenhuma depende de conhecimento previo (embeddings,
rede, DepMap) ou cai na resposta media.

Entradas: data/raw/pie/<repo>/meta.json, ctx_ids.npy, pert_ids.npy, baixados de
huggingface.co/datasets/arcinstitute/<repo>/preprocessed. Os .npy sao lidos
com allow_pickle=False: .npy vindo da rede com pickle e execucao de codigo.
"""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd

ALVOS = Path("data/raw/vcc/controls-final/pert_counts.csv")
GENES_VCC = Path("data/raw/vcc/controls-final/gene_names.csv")
PIE = Path("data/raw/pie")
SAIDA = Path("runs/vcc/caracterizacao/sobreposicao_pie.tsv")
REPOS = ("PIE_arc_vcc_25", "PIE_replogle_nadig_essential", "PIE_x_atlas_orion", "PIE_jiang")


def main() -> None:
    alvos = set(pd.read_csv(ALVOS)["target_gene"].astype(str))
    genes_vcc = set(pd.read_csv(GENES_VCC)["gene_name"].astype(str))
    linhas: list[dict[str, str | int]] = []

    for repo in REPOS:
        base = PIE / "arrays" / repo / "preprocessed"
        meta = json.loads((base / "meta.json").read_text(encoding="utf-8"))
        ctx_ids = np.load(base / "ctx_ids.npy", allow_pickle=False)
        pert_ids = np.load(base / "pert_ids.npy", allow_pickle=False)
        id_para_ctx = {v: k for k, v in meta["context_to_id"].items()}
        id_para_pert = {v: k for k, v in meta["pert_to_id"].items()}

        genes_pie = set(meta["genes"])
        print(f"\n== {repo} ==")
        print(f"  eixo de genes: {len(genes_pie)} | em comum com o painel VCC: {len(genes_pie & genes_vcc)}")
        todos = {id_para_pert[int(p)] for p in pert_ids}
        print(f"  alvos do VCC com efeito medido em alguma linhagem: {len(alvos & todos)} de {len(alvos)}")

        por_ctx: dict[str, set[str]] = {}
        for c, p in zip(ctx_ids.tolist(), pert_ids.tolist(), strict=True):
            por_ctx.setdefault(id_para_ctx[int(c)], set()).add(id_para_pert[int(p)])
        for ctx, perts in sorted(por_ctx.items()):
            comuns = sorted(alvos & perts)
            print(f"    {ctx:12s} {len(perts):5d} perturbacoes | {len(comuns):3d} em comum")
            for alvo in comuns:
                linhas.append({"repo": repo, "contexto": ctx, "alvo": alvo})

    cobertura = pd.DataFrame(linhas)
    SAIDA.parent.mkdir(parents=True, exist_ok=True)
    cobertura.to_csv(SAIDA, sep="\t", index=False)

    cobertos = set(cobertura["alvo"]) if len(cobertura) else set()
    sem_medida = sorted(alvos - cobertos)
    print(f"\n== Total: {len(cobertos)} de {len(alvos)} alvos com efeito medido em pelo menos uma fonte ==")
    print(f"   sem nenhuma medida ({len(sem_medida)}): {', '.join(sem_medida[:40])}{' ...' if len(sem_medida) > 40 else ''}")


if __name__ == "__main__":
    main()
