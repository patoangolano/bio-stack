# scripts/vcc/avaliador_recorte.py
"""Recorta uma submissao real para o formato que o `avaliador_proxy` pontua.

Etapa 5 de docs/VCC-2026-ROTA.md. O `submissao_v*.h5ad` do `emissor_400` tem
360 000 celulas, 300 alvos x 3 contextos, e nenhuma linha `non-targeting` —
exatamente o contrato do desafio. O avaliador precisa de outra coisa:

1. **O mesmo conjunto de perturbacoes nos dois lados.** `cell_eval2.io.
   validate_pair` aborta quando os conjuntos diferem, e aborta inclusive por
   causa do rotulo `non-targeting`, que o contrato de submissao proibe. O passo
   de preparacao da Arc reinjeta os controles reais na predicao; este script faz
   o mesmo, copiando o pool de controle do proprio lado real.
2. **Um recorte que caiba em 15 GiB.** Medido: 18 400 controles + 12 800 celulas
   perturbadas (32 alvos) levam o passo `baseline` a 10,9 GiB. Os 120 000 de um
   contexto inteiro nao cabem.

O lado real contra o qual o recorte e pontuado vem de
`scripts/vcc/avaliador_proxy.py --so-alvos-vcc`, que rotula as perturbacoes com
nomes do painel VCC e aplica a eles efeitos MEDIDOS pela Arc em `ARC_H1`.

⚠️ O QUE ISSO MEDE, E O QUE NAO MEDE. A predicao e nossa de verdade — saiu do
`transfer_fit` e do `emissor_400`. A "verdade" nao: e o efeito de `ARC_H1`
(celula H1 embrionaria) aplicado a celulas de D/E/F. Entao a pontuacao responde
"quanto a nossa predicao concorda com uma resposta de formato H1", nao "quanto
ela acerta D". Serve para dois propositos, ambos reais: provar que o caminho de
pontuacao aceita a submissao de verdade ponta a ponta, e dar um numero de
referencia contra o qual comparar a PROXIMA versao do modelo no mesmo proxy.
Nao serve para estimar o placar.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import scipy.sparse as sp

COLUNA_PERT = "target_gene"
COLUNA_CTX = "context"
ROTULO_CONTROLE = "non-targeting"


def _rotulos(adata: ad.AnnData, coluna: str) -> np.ndarray:
    return adata.obs[coluna].astype(str).to_numpy()


def main() -> int:
    p = argparse.ArgumentParser(description=__doc__ or "",
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--submissao", required=True, help="submissao_v*.h5ad do emissor_400")
    p.add_argument("--proxy-real", required=True,
                   help="real.h5ad de avaliador_proxy.py; define os alvos e o pool de controle")
    p.add_argument("--contexto", default="D", choices=["D", "E", "F"])
    p.add_argument("--saida", required=True, help="pred.h5ad a gravar")
    p.add_argument("--resumo", default=None, help="JSON com o que foi recortado")
    args = p.parse_args()

    # O lado real primeiro: ele decide quais alvos entram e quais celulas de
    # controle a predicao carrega. Ler o real BACKED nao serve — as linhas de
    # controle sao copiadas para a predicao, entao a matriz e necessaria.
    real = ad.read_h5ad(args.proxy_real)
    rot_real = _rotulos(real, COLUNA_PERT)
    e_ctrl = rot_real == ROTULO_CONTROLE
    alvos: list[str] = sorted({str(a) for a in rot_real[~e_ctrl]})
    if not alvos:
        raise SystemExit(f"{args.proxy_real} nao tem perturbacao alguma fora do controle")

    sub = ad.read_h5ad(args.submissao, backed="r")
    if list(sub.var_names) != list(real.var_names):
        raise SystemExit("o eixo de genes da submissao nao bate com o do lado real")
    if COLUNA_CTX not in sub.obs.columns:
        raise SystemExit(f"a submissao nao tem a coluna {COLUNA_CTX!r}; "
                         f"presentes: {list(sub.obs.columns)}")

    sub_ctx = _rotulos(sub, COLUNA_CTX)
    sub_pert = _rotulos(sub, COLUNA_PERT)
    conjunto: set[str] = set(alvos)
    mascara = (sub_ctx == args.contexto) & np.isin(sub_pert, alvos)
    faltam = sorted(conjunto - set(sub_pert[sub_ctx == args.contexto]))
    if faltam:
        raise SystemExit(
            f"a submissao nao cobre {len(faltam)} alvo(s) do lado real no contexto "
            f"{args.contexto}: {faltam[:10]}. Rode o proxy com --so-alvos-vcc para que os "
            "rotulos sejam alvos do painel."
        )

    # Indices ORDENADOS: o h5py recusa indexacao fancy fora de ordem, e a CSR do
    # emissor e escrita bloco por bloco, logo as linhas de um (contexto, alvo) ja
    # sao contiguas — ordenar nao embaralha nada.
    linhas = np.flatnonzero(mascara)
    recorte = sub[linhas].to_memory()
    X_pred = sp.csr_matrix(recorte.X)
    pert_pred = _rotulos(recorte, COLUNA_PERT)
    del recorte, sub

    X_ctrl = sp.csr_matrix(real[e_ctrl].X)
    n_ctrl = int(e_ctrl.sum())
    del real

    X = sp.vstack([X_ctrl, X_pred], format="csr").astype(np.float32)
    pert = np.concatenate([np.full(n_ctrl, ROTULO_CONTROLE, dtype=object), pert_pred])
    obs = pd.DataFrame(
        {COLUNA_PERT: pd.Categorical(pert)},
        index=pd.Index([f"c{i:07d}" for i in range(pert.size)], name=None),
    )
    eixo = pd.read_csv(
        Path(__file__).resolve().parents[2] / "data" / "raw" / "vcc" /
        "controls-final" / "gene_names.csv")["gene_name"].astype(str).to_numpy()
    pred = ad.AnnData(X=X, obs=obs, var=pd.DataFrame(index=pd.Index(eixo, name=None)))

    saida = Path(args.saida)
    saida.parent.mkdir(parents=True, exist_ok=True)
    pred.write_h5ad(saida, compression=None)

    _, contagens = np.unique(pert_pred, return_counts=True)
    resumo = {
        "submissao": str(args.submissao),
        "proxy_real": str(args.proxy_real),
        "contexto": args.contexto,
        "alvos": alvos,
        "n_alvos": len(alvos),
        "celulas_por_alvo": sorted({int(c) for c in contagens}),
        "celulas_perturbadas": int(pert_pred.size),
        "celulas_de_controle_copiadas_do_real": n_ctrl,
        "celulas_total": int(pert.size),
        "nnz": int(X.nnz),
        "nnz_por_celula": round(X.nnz / pert.size, 1),
        "saida": str(saida),
    }
    if args.resumo:
        Path(args.resumo).write_text(json.dumps(resumo, indent=2, ensure_ascii=False) + "\n",
                                     encoding="utf-8")
    print(json.dumps(resumo, indent=2, ensure_ascii=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
