"""A razao `orcamento/pedido` do `expr_mse`: nossa correcao de amostragem e paga ou descontada?

## O que esta em jogo

`expr_mse_unbiased_capped_norm` subtrai da distancia predita uma **correcao de
amostragem** `C_p` — o jackknife delete-1 do pseudobulk do grupo, que estima
quanto da distancia e so ruido de amostrar 400 celulas. Sem ela a predicao paga
por ruido que nao e erro.

Mas a correcao e **limitada**. O `cell-eval2` 0.19.0 compara o que a submissao
PEDE com um orcamento derivado da propria submissao (issue #348):

    row_w  = 1 / n_row              n_row = G-1 na linha cujo alvo resolve, G se nao
    pedido = sum_p C_p * row_w[p]   com C_p ja limitado por min(C_p, K * C_p_real)
    orcamento = sum_g sum_{p in R_g} w_p (mu_pred[p,g] - wmean mu_pred[.,g])**2
    razao  = orcamento / pedido
    r      = min(1, razao)          e `C_p` efetivo e `r * C_p`

O orcamento e a dispersao ENTRE perturbacoes da predicao. A logica: se a
predicao tem pouca dispersao entre perturbacoes mas pede uma correcao grande,
o espalhamento por celula dela CANCELA no pseudobulk — que nao e o que uma
amostra i.i.d. de uma populacao predita de celulas parece. Abaixo de
`PRED_CORRECTION_BUDGET_FLAG_RATIO = 0,7` o avaliador registra um aviso.

Registrado das corridas de 22/10: v4 (k=1) razao **0,179**, abaixo do flag; v5
(k=4) razao **2,854**, `r = 1,0`, limpa. A diferenca e `k_desvio`: subir a
dispersao entre perturbacoes sobe o orcamento e nao sobe `C_p`, que depende do
espalhamento DENTRO do grupo. Entao a mesma escolha que resolveu `fidelity`
resolveu isto de lado — e vale conferir em cada submissao nova, porque e
possivel sair do regime seguro sem nada falhar.

## Por que a razao medida aqui e um LIMITE INFERIOR, e por que isso basta

O `pedido` do avaliador usa `C_p` **ja limitado**:

    pred_term = np.minimum(pred_term, PRED_TRACE_CAP_K * real_tn[rows])

`real_tn` e a correcao do LADO REAL, que em D/E/F e retido e nao existe aqui.
Mas limitar so pode REDUZIR o pedido, e reduzir o pedido AUMENTA a razao. Logo:

    razao_medida_aqui  <=  razao_do_avaliador

Entao `razao_sem_limite >= 1` certifica `r = 1` sem o lado real: a correcao esta
sendo paga por inteiro. E um certificado de um lado so, e e o lado que importa.
Abaixo de 0,7 o veredito e "em risco" e nao "sinalizado", porque o limite ainda
pode resgatar — e o relatorio diz qual dos dois.

## Nada e portado

`jackknife_correction`, `_across_pert_budget` e `_row_weights` sao CHAMADOS do
`cell_eval2` instalado em `envs/cell-eval2`. Uma reimplementacao concordaria com
o original ate o dia em que o original mudasse, e o erro apareceria como uma
razao plausivel e errada. Por isso este script roda no python daquele ambiente:

    ./envs/cell-eval2/.venv/Scripts/python.exe -I scripts/vcc/orcamento.py \
        --submissao runs/vcc/submissao/submissao_v6.h5ad

## Memoria: um grupo por vez, e isso e exato

`jackknife_correction` e `O(n*G)` DENSO em blocos de linhas, e a submissao tem
682 milhoes de nao nulos por contexto — 5,5 GB em CSR, que nao cabe junto com os
temporarios. Mas o jackknife de um grupo so olha aquele grupo, entao varrer
grupo a grupo da o MESMO vetor: verificado identico bit a bit contra a chamada
com todos os grupos juntos (`test_jackknife_por_grupo_e_identico_ao_conjunto`).
Cada grupo e ~2,3 milhoes de nao nulos, ~18 MB.

O `chunk` fica no padrao 512 de proposito. O docstring do `jackknife_correction`
mede 32 como mais rapido E menor (89 MiB / 116 ms contra 349 MiB / 375 ms), e
diz que o bloco e invariante na resposta a ~1 ULP. Mas aqui o objetivo e
reproduzir o numero que o avaliador vai calcular, e o avaliador usa o padrao.
Fidelidade antes de velocidade quando o produto e uma auditoria.
"""

from __future__ import annotations

import argparse
import json
import sys
from dataclasses import asdict, dataclass
from pathlib import Path

import anndata as ad
import h5py
import numpy as np
import pandas as pd
import scipy.sparse as sp
from cell_eval2 import moments
from cell_eval2.metrics import delta

sys.path.insert(0, str(Path(__file__).resolve().parent))
from pseudobulk import ALVO_SOMA, bulk_lognorm  # noqa: E402

CONTROLE = "non-targeting"
QUEM = "orcamento"


@dataclass(frozen=True)
class Laudo:
    """O veredito de um contexto."""

    contexto: str
    grupos: int
    genes: int
    alvos_pontuados: int
    alvos_com_alvo_resolvido: int
    orcamento: float
    pedido_sem_limite: float
    razao_sem_limite: float
    r_sem_limite: float
    flag: float
    veredito: str
    jk_mediano: float
    jk_min: float
    jk_max: float


def grupos_do_contexto(caminho: Path, contexto: str | None) -> tuple[list[str], np.ndarray, list[str]]:
    """Rotulos de grupo em ordem de aparicao, codigo por celula e eixo de genes.

    O codigo e por celula e nao por bloco porque a submissao vir ordenada por
    grupo e uma propriedade do nosso emissor, nao do formato: um arquivo de
    terceiro pode intercalar, e `np.flatnonzero` sobre o codigo resolve os dois.
    """
    with h5py.File(caminho, "r") as f:
        g_obs, g_var = f["obs"], f["var"]
        assert isinstance(g_obs, h5py.Group) and isinstance(g_var, h5py.Group)
        obs, var = ad.io.read_elem(g_obs), ad.io.read_elem(g_var)
    assert isinstance(obs, pd.DataFrame) and isinstance(var, pd.DataFrame)
    genes = var.index.astype(str).tolist()
    if contexto is not None:
        if "context" not in obs.columns:
            raise ValueError(f"{caminho}: pediu o contexto {contexto} e obs nao tem a coluna 'context'")
        obs = obs.loc[obs["context"].astype(str) == contexto]
        if obs.empty:
            raise ValueError(f"{caminho}: nenhuma celula no contexto {contexto}")
    rotulos_por_celula = obs["target_gene"].astype(str).to_numpy()
    codigos, unicos = pd.factorize(rotulos_por_celula)
    return [str(u) for u in unicos], np.asarray(codigos), genes


def varrer(
    caminho: Path, contexto: str | None, chunk: int
) -> tuple[list[str], np.ndarray, np.ndarray, list[str]]:
    """Pseudobulk `bulk_lognorm` e jackknife `C_p` por grupo, um grupo por vez.

    Devolve (rotulos, pred_means [P, G], jk [P], genes).
    """
    rotulos, codigos, genes = grupos_do_contexto(caminho, contexto)
    # As linhas do contexto pedido, no indice do arquivo inteiro. O recorte de
    # obs acima perdeu a numeracao original, entao ela e refeita aqui.
    with h5py.File(caminho, "r") as f:
        g_obs = f["obs"]
        assert isinstance(g_obs, h5py.Group)
        obs_todo = ad.io.read_elem(g_obs)
        assert isinstance(obs_todo, pd.DataFrame)
        if contexto is None:
            linhas_do_contexto = np.arange(len(obs_todo), dtype=np.int64)
        else:
            linhas_do_contexto = np.flatnonzero(
                obs_todo["context"].astype(str).to_numpy() == contexto
            ).astype(np.int64)

        x = f["X"]
        if not isinstance(x, h5py.Group):
            raise ValueError(f"{caminho}: X nao e esparso; rode validar_submissao.py primeiro")
        dados, indices = x["data"], x["indices"]
        assert isinstance(dados, h5py.Dataset) and isinstance(indices, h5py.Dataset)
        indptr = np.asarray(x["indptr"]).astype(np.int64)
        n_genes = len(genes)

        somas = np.zeros((len(rotulos), n_genes), dtype=np.float64)
        jk = np.zeros(len(rotulos), dtype=np.float64)
        for p in range(len(rotulos)):
            linhas = linhas_do_contexto[codigos == p]
            if linhas.size == 0:
                continue
            # Um grupo do nosso emissor e um bloco contiguo; de um arquivo
            # qualquer pode nao ser, e aí a leitura vai por fatia minima
            # cobrindo as linhas e descarta o resto.
            r0, r1 = int(linhas.min()), int(linhas.max()) + 1
            a, b = int(indptr[r0]), int(indptr[r1])
            bloco = sp.csr_matrix(
                (
                    np.asarray(dados[a:b], dtype=np.float64),
                    np.asarray(indices[a:b]),
                    indptr[r0 : r1 + 1] - a,
                ),
                shape=(r1 - r0, n_genes),
            )
            if linhas.size != r1 - r0:
                # `sp.csr_matrix(...)` em volta do recorte: o verificador nao
                # estreita o `__getitem__` de matriz esparsa e le `None` como
                # possivel. Mesmo idioma de `pseudobulk.somas_por_grupo`.
                bloco = sp.csr_matrix(bloco[linhas - r0])
            somas[p] = np.asarray(bloco.sum(axis=0)).ravel()
            # `linhas.size` e nao `bloco.shape[0]`: depois do recorte acima o
            # verificador le a forma da esparsa como possivelmente `None`, e os
            # dois valores sao o mesmo numero por construcao nos dois ramos.
            correcao = np.asarray(
                moments.jackknife_correction(
                    bloco,
                    np.zeros(int(linhas.size), dtype=np.int64),
                    1,
                    ALVO_SOMA,
                    chunk=chunk,
                ),
                dtype=np.float64,
            )
            jk[p] = float(correcao[0])
            if (p + 1) % 50 == 0:
                print(f"    {p + 1}/{len(rotulos)} grupos", flush=True)

    return rotulos, bulk_lognorm(somas), jk, genes


def auditar(caminho: Path, contexto: str | None, chunk: int) -> Laudo:
    """Monta orcamento e pedido chamando o `cell_eval2` instalado."""
    rotulos, pred_means, jk, genes = varrer(caminho, contexto, chunk)
    n_genes = len(genes)

    # `keep`: o que o avaliador pontua — tudo menos o controle. Nossa submissao
    # nao emite controle, mas um arquivo de terceiro pode.
    keep = np.array([i for i, r in enumerate(rotulos) if r != CONTROLE], dtype=np.intp)
    if keep.size == 0:
        raise ValueError(f"{caminho}: nenhuma perturbacao fora do controle")
    rotulos_keep = [rotulos[i] for i in keep]

    # Do `cell_eval2`, nao reimplementado. `target_gene_map=None` deixa o
    # resolvedor casar o rotulo com o eixo de genes, que e o caso do painel VCC
    # (o alvo E um gene do eixo). `gate_labels` e o conjunto julgado pelo raise
    # de "nenhum alvo resolve", e aqui ele e o proprio conjunto pontuado.
    excl = delta._exclusion_cols(  # noqa: SLF001
        rotulos_keep,
        genes,
        target_gene_map=None,
        gate_labels=rotulos_keep,
        who=QUEM,
        n_genes=n_genes,
    )
    row_w = delta._row_weights(keep.size, n_genes, excl)  # noqa: SLF001
    orcamento = float(delta._across_pert_budget(pred_means[keep], excl, row_w))  # noqa: SLF001
    # SEM o limite: ver "Por que a razao medida aqui e um LIMITE INFERIOR".
    pedido = float((jk[keep] * row_w).sum())
    razao = orcamento / pedido if pedido > 0 else float("inf")
    flag = float(delta.PRED_CORRECTION_BUDGET_FLAG_RATIO)
    if razao >= 1.0:
        veredito = "paga por inteiro (r = 1,0 garantido, pois o limite so aumenta a razao)"
    elif razao >= flag:
        veredito = f"descontada: r <= {min(1.0, razao):.4f}, mas acima do flag de {flag:.2f}"
    else:
        veredito = (
            f"EM RISCO: razao sem limite {razao:.4f} abaixo do flag de {flag:.2f}. O limite "
            "pelo lado real pode resgatar, e so o avaliador com o lado real decide"
        )

    return Laudo(
        contexto=contexto or "(todos)",
        grupos=len(rotulos),
        genes=n_genes,
        alvos_pontuados=int(keep.size),
        alvos_com_alvo_resolvido=len(excl),
        orcamento=orcamento,
        pedido_sem_limite=pedido,
        razao_sem_limite=razao,
        r_sem_limite=min(1.0, razao),
        flag=flag,
        veredito=veredito,
        jk_mediano=float(np.median(jk[keep])),
        jk_min=float(jk[keep].min()),
        jk_max=float(jk[keep].max()),
    )


def main(argv: list[str] | None = None) -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--submissao", type=Path, required=True)
    ap.add_argument(
        "--contextos",
        nargs="*",
        default=None,
        help="contextos a auditar; omitir trata o arquivo inteiro como um conjunto",
    )
    ap.add_argument(
        "--chunk",
        type=int,
        default=512,
        help="bloco de linhas do jackknife. 512 e o padrao do cell-eval2, e e o que ele usara",
    )
    ap.add_argument("--relatorio", type=Path, default=None)
    args = ap.parse_args(argv)

    contextos: list[str | None] = list(args.contextos) if args.contextos else [None]
    laudos: list[Laudo] = []
    for contexto in contextos:
        print(f"== {args.submissao.name} contexto {contexto or '(todos)'}", flush=True)
        laudo = auditar(args.submissao, contexto, args.chunk)
        laudos.append(laudo)
        print(
            f"  {laudo.alvos_pontuados} alvos pontuados, {laudo.alvos_com_alvo_resolvido} com "
            f"o proprio alvo resolvido no eixo"
        )
        print(f"  C_p (jackknife): mediana {laudo.jk_mediano:.5f} "
              f"[{laudo.jk_min:.5f}, {laudo.jk_max:.5f}]")
        print(f"  orcamento        {laudo.orcamento:.6g}")
        print(f"  pedido sem limite {laudo.pedido_sem_limite:.6g}")
        print(f"  razao            {laudo.razao_sem_limite:.4f}  ->  r <= {laudo.r_sem_limite:.4f}")
        print(f"  {laudo.veredito}")

    if args.relatorio:
        args.relatorio.parent.mkdir(parents=True, exist_ok=True)
        args.relatorio.write_text(
            json.dumps(
                {
                    "submissao": str(args.submissao),
                    "chunk": args.chunk,
                    "bulk_target_sum": ALVO_SOMA,
                    "pred_trace_cap_k": float(delta.PRED_TRACE_CAP_K),
                    "flag": float(delta.PRED_CORRECTION_BUDGET_FLAG_RATIO),
                    "chamadas_do_cell_eval2": [
                        "moments.jackknife_correction",
                        "metrics.delta._exclusion_cols",
                        "metrics.delta._row_weights",
                        "metrics.delta._across_pert_budget",
                    ],
                    "o_limite_e_de_um_lado_so": (
                        "O pedido do avaliador usa C_p ja limitado por min(C_p, K*C_p_real), e "
                        "real_tn nao existe para D/E/F. Limitar reduz o pedido e portanto AUMENTA "
                        "a razao, entao a razao medida aqui e um limite inferior da do avaliador: "
                        "razao >= 1 certifica r = 1; razao < 0,7 diz 'em risco', nao 'sinalizado'."
                    ),
                    "laudos": [asdict(x) for x in laudos],
                },
                indent=2,
                ensure_ascii=False,
            ),
            encoding="utf-8",
        )
        print(f"\nrelatorio em {args.relatorio}")


if __name__ == "__main__":
    main()
