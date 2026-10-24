"""Quantos genes significativos a submissao realmente produz, pelo motor do avaliador.

## O risco aberto que isto fecha

A `varredura_k.py` conta genes com `|delta| > 0,01` e chama o numero de
aproximacao declarada: em k=4 deu 5.532, "2,02x os 2.739 reais". Mas `|delta|`
nao e o teste de DE. O avaliador roda Wilcoxon por gene contra os controles,
aplica o gate de CPM e corrige por BH dentro de cada alvo. Uma predicao pode ter
delta grande em muitos genes e `nsig` baixo (se o espalhamento por celula
afogar o sinal) ou delta modesto e `nsig` altissimo (se as celulas forem quase
identicas entre si). Os dois casos pontuam mal em `sig_jaccard` e em
`direction_reach`, e nenhum dos dois aparece numa contagem por limiar.

Este modulo mede `nsig` com o motor que vai pontuar.

## E e exato, nao um proxy — tres propriedades do preset fazem isso

1. **`control_source: real`.** O lado de controle do teste de DE vem do dado
   REAL, e os controles de D/E/F estao publicados em `controls-final`. Entao o
   teste que o avaliador roda na nossa submissao e inteiramente reproduzivel
   aqui: nossas 400 celulas preditas do alvo X contra os 18.400 controles reais
   do contexto. Nao ha lado retido nesta conta.

2. **`fdr_scope: per_pert`.** `de._bh_per_target` refaz o BH dentro de cada alvo
   sobre os p-valores DELE. Entao o `nsig` de um alvo nao depende de quais
   outros alvos estao na corrida, e amostrar alvos e EXATO para os amostrados —
   nao uma estimativa. E o que torna a auditoria viavel: 300 alvos x 3 contextos
   de Wilcoxon contra 18.400 celulas sao horas, e uma amostra responde a mesma
   pergunta com os mesmos numeros nos alvos que entraram.

3. **`filter_gene_min_cpm_cell: 5,0` com UM conjunto de genes para todos os
   alvos** (#351). O gate e `media aritmetica do CPM por celula da REFERENCIA >
   5,0`, e a referencia e o controle real. Logo o universo e fixo por contexto e
   independe da nossa predicao — o que casa com os 9.633 / 9.743 / 9.537
   medidos em `ctrl_bulk`. Um gate que dependesse da predicao poderia ser
   movido por ela, e este nao pode.

## Nada e portado, de novo

`de_compute.compute_de`, `de.prep_de_side` e `de.rank_de_side` sao CHAMADOS do
`cell_eval2` instalado. O `nsig` sai de duas contagens independentes sobre a
mesma tabela — `p_adj < 0,05` agrupado por alvo, e a largura da coluna do alvo
na matriz de postos que os membros de DE consomem. As duas tem de concordar; se
discordarem, o relatorio diz, porque ai a definicao de "significativo" que a
auditoria usa nao e a que o placar usa.

    ./envs/cell-eval2/.venv/Scripts/python.exe -I scripts/vcc/auditoria.py \
        --submissao runs/vcc/submissao/submissao_v6.h5ad --contexto D --n-alvos 40

## O numero de referencia

Mediana de ~2.739 genes significativos por perturbacao no ensaio real do
arc_vcc_25. Nao e o alvo a atingir — e a ordem de grandeza contra a qual uma
emissao de 400 celulas sinteticas pode estar absurda em qualquer das direcoes.
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
import polars as pl
import scipy.sparse as sp
from cell_eval2 import de, de_compute

sys.path.insert(0, str(Path(__file__).resolve().parent))

CONTROLE = "non-targeting"
COLUNA_ALVO = "target_gene"
# Do `configs/vcc2026.yaml`, verbatim. Copiados e nao lidos porque o preset e um
# YAML interno do pacote e a chave pode mudar de lugar entre versoes; o teste
# `test_os_parametros_batem_o_preset_instalado` compara os dois e falha se
# divergirem, que e a defesa certa para uma constante duplicada.
P_ADJ = 0.05
MIN_ABS_LOG2FC = 0.0
SORT_BY = "abs_log2_fold_change"
METODO = "wilcoxon"
NAN_LFC = "mask"
MEAN_CALC = "arithmetic"
EPSILON = 1.0e-09
FDR_SCOPE = "per_pert"
GATE_CPM = 5.0
TARGET_SUM = 1_000_000.0
TIPO_DE_ENTRADA = "counts"
# Mediana de genes significativos por perturbacao no ensaio real (arc_vcc_25).
DE_REAIS = 2739


@dataclass(frozen=True)
class Laudo:
    contexto: str
    alvos_auditados: int
    alvos_no_arquivo: int
    amostra: bool
    semente: int
    genes_no_eixo: int
    genes_apos_o_gate: int
    nsig_mediano: float
    nsig_p05: float
    nsig_p25: float
    nsig_p75: float
    nsig_p95: float
    nsig_min: int
    nsig_max: int
    alvos_sem_nenhum_gene: int
    razao_para_de_real: float
    duas_contagens_concordam: bool
    maior_divergencia: int


def ler_celulas(caminho: Path, contexto: str | None, alvos: set[str] | None) -> ad.AnnData:
    """Le do `.h5ad` apenas as celulas pedidas, como AnnData denso-esparso em memoria.

    Le por fatia de linhas e nao com `ad.read_h5ad`: a submissao inteira tem 2,1
    bilhoes de nao nulos e 16 GB em disco, e um alvo sao 2,3 milhoes.
    """
    with h5py.File(caminho, "r") as f:
        g_obs, g_var = f["obs"], f["var"]
        assert isinstance(g_obs, h5py.Group) and isinstance(g_var, h5py.Group)
        obs, var = ad.io.read_elem(g_obs), ad.io.read_elem(g_var)
        assert isinstance(obs, pd.DataFrame) and isinstance(var, pd.DataFrame)

        quer = np.ones(len(obs), dtype=bool)
        if contexto is not None:
            if "context" not in obs.columns:
                raise ValueError(f"{caminho}: pediu contexto {contexto} e obs nao tem 'context'")
            quer &= obs["context"].astype(str).to_numpy() == contexto
        if alvos is not None:
            quer &= np.isin(obs[COLUNA_ALVO].astype(str).to_numpy(), list(alvos))
        linhas = np.flatnonzero(quer).astype(np.int64)
        if linhas.size == 0:
            raise ValueError(f"{caminho}: nenhuma celula no recorte pedido")

        x = f["X"]
        if not isinstance(x, h5py.Group):
            raise ValueError(f"{caminho}: X nao e esparso; rode validar_submissao.py primeiro")
        dados, indices = x["data"], x["indices"]
        assert isinstance(dados, h5py.Dataset) and isinstance(indices, h5py.Dataset)
        indptr = np.asarray(x["indptr"]).astype(np.int64)
        n_genes = var.shape[0]

        # Fatias contiguas maximais: as celulas de um grupo saem juntas do nosso
        # emissor, entao 300 alvos pedidos viram ~300 leituras e nao 120.000.
        quebras = np.flatnonzero(np.diff(linhas) != 1) + 1
        pedacos: list[sp.csr_matrix] = []
        for parte in np.split(linhas, quebras):
            r0, r1 = int(parte[0]), int(parte[-1]) + 1
            a, b = int(indptr[r0]), int(indptr[r1])
            pedacos.append(
                sp.csr_matrix(
                    (
                        np.asarray(dados[a:b], dtype=np.float32),
                        np.asarray(indices[a:b]),
                        indptr[r0 : r1 + 1] - a,
                    ),
                    shape=(r1 - r0, n_genes),
                )
            )
        matriz = sp.vstack(pedacos, format="csr") if len(pedacos) > 1 else pedacos[0]

    obs_recorte = obs.iloc[linhas].copy()
    # `compute_de` agrupa por string; categoria com nivel nao presente viraria
    # grupo vazio e o backend reclamaria de um grupo sem celula.
    obs_recorte[COLUNA_ALVO] = obs_recorte[COLUNA_ALVO].astype(str)
    return ad.AnnData(X=matriz, obs=obs_recorte, var=var.copy())


def contar_por_alvo(tabela: pl.DataFrame) -> tuple[dict[str, int], dict[str, int]]:
    """Duas contagens independentes de `nsig`, pela cadeia do pacote.

    A primeira filtra `p_adj < P_ADJ` e agrupa. A segunda vai pela
    `de.rank_de_side`, que e a matriz de postos que os membros de DE realmente
    consomem, e conta a altura nao nula da coluna do alvo. Devolver as duas e o
    que permite afirmar que a definicao usada aqui e a do placar.
    """
    preparada, _ = de.prep_de_side(
        tabela, name="pred", sort_by=SORT_BY, nan_lfc_policy=NAN_LFC, min_abs_log2fc=MIN_ABS_LOG2FC
    )
    direto = {
        str(r["target"]): int(r["len"])
        for r in preparada.filter(pl.col("p_adj") < P_ADJ)
        .group_by("target")
        .len()
        .iter_rows(named=True)
    }
    postos = de.rank_de_side(preparada, sort_by=SORT_BY, p_adj_threshold=P_ADJ)
    pela_matriz: dict[str, int] = {}
    if postos.height:
        for coluna in postos.columns:
            pela_matriz[str(coluna)] = int(postos[coluna].len() - postos[coluna].null_count())
    return direto, pela_matriz


def auditar(
    submissao: Path,
    controles: Path,
    contexto: str,
    n_alvos: int,
    semente: int,
    threads: int,
) -> tuple[Laudo, dict[str, int]]:
    with h5py.File(submissao, "r") as f:
        g_obs = f["obs"]
        assert isinstance(g_obs, h5py.Group)
        obs = ad.io.read_elem(g_obs)
    assert isinstance(obs, pd.DataFrame)
    do_contexto = obs.loc[obs["context"].astype(str) == contexto] if "context" in obs.columns else obs
    todos = sorted({str(a) for a in do_contexto[COLUNA_ALVO].astype(str) if str(a) != CONTROLE})
    if not todos:
        raise ValueError(f"{submissao}: nenhum alvo no contexto {contexto}")

    amostra = 0 < n_alvos < len(todos)
    if amostra:
        rng = np.random.default_rng(semente)
        escolhidos = sorted(rng.choice(np.array(todos, dtype=object), size=n_alvos, replace=False))
    else:
        escolhidos = todos
    print(f"  {len(escolhidos)} de {len(todos)} alvos{' (amostra)' if amostra else ''}", flush=True)

    predito = ler_celulas(submissao, contexto, set(escolhidos))
    print(f"  predito: {predito.n_obs} celulas x {predito.n_vars} genes", flush=True)
    referencia = ler_celulas(controles, None, None)
    print(f"  controle real: {referencia.n_obs} celulas", flush=True)
    if list(predito.var_names) != list(referencia.var_names):
        raise ValueError("eixo de genes diferente entre a submissao e o controle real")

    # CONCATENAR e nao passar `reference=` como AnnData: o `compute_de` aceita um
    # pool de controle externo em memoria SO no backend gpudge, e levanta
    # `"pdex/scanpy must concatenate"` nos outros dois. Concatenar e exatamente o
    # que a corrida real faz sob `control_source: real`, entao o teste e o mesmo.
    referencia.obs[COLUNA_ALVO] = CONTROLE
    if "context" in predito.obs.columns:
        referencia.obs["context"] = predito.obs["context"].iloc[0]
    junto = ad.concat([predito, referencia], join="outer", index_unique="-")
    junto.obs[COLUNA_ALVO] = junto.obs[COLUNA_ALVO].astype(str)

    tabela = de_compute.compute_de(
        junto,
        backend="scanpy",
        groupby=COLUNA_ALVO,
        reference=CONTROLE,
        control_group=CONTROLE,
        mean_calc=MEAN_CALC,
        epsilon=EPSILON,
        input_type=TIPO_DE_ENTRADA,
        target_sum=TARGET_SUM,
        filter_gene_min_cpm_cell=GATE_CPM,
        fdr_scope=FDR_SCOPE,
        threads=threads,
    )
    genes_apos_gate = int(tabela["feature"].n_unique())
    direto, pela_matriz = contar_por_alvo(tabela)

    nsig = np.array([direto.get(a, 0) for a in escolhidos], dtype=np.int64)
    divergencias = [abs(direto.get(a, 0) - pela_matriz.get(a, 0)) for a in escolhidos]
    laudo = Laudo(
        contexto=contexto,
        alvos_auditados=len(escolhidos),
        alvos_no_arquivo=len(todos),
        amostra=amostra,
        semente=semente,
        genes_no_eixo=int(predito.n_vars),
        genes_apos_o_gate=genes_apos_gate,
        nsig_mediano=float(np.median(nsig)),
        nsig_p05=float(np.percentile(nsig, 5)),
        nsig_p25=float(np.percentile(nsig, 25)),
        nsig_p75=float(np.percentile(nsig, 75)),
        nsig_p95=float(np.percentile(nsig, 95)),
        nsig_min=int(nsig.min()),
        nsig_max=int(nsig.max()),
        alvos_sem_nenhum_gene=int((nsig == 0).sum()),
        razao_para_de_real=float(np.median(nsig)) / DE_REAIS,
        duas_contagens_concordam=max(divergencias) == 0,
        maior_divergencia=int(max(divergencias)),
    )
    return laudo, {a: int(v) for a, v in zip(escolhidos, nsig, strict=True)}


def main(argv: list[str] | None = None) -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--submissao", type=Path, required=True)
    ap.add_argument(
        "--controles",
        type=Path,
        default=None,
        help="o context_<C>.h5ad do controls-final; o padrao segue o contexto pedido",
    )
    ap.add_argument("--contexto", default="D")
    ap.add_argument(
        "--n-alvos",
        type=int,
        default=40,
        help="alvos a auditar; 0 roda todos. Amostrar e EXATO sob fdr_scope=per_pert",
    )
    ap.add_argument("--semente", type=int, default=0)
    ap.add_argument("--threads", type=int, default=-1)
    ap.add_argument("--relatorio", type=Path, default=None)
    args = ap.parse_args(argv)

    controles = args.controles or (
        Path("data/raw/vcc/controls-final") / f"context_{args.contexto}.h5ad"
    )
    if not controles.exists():
        raise FileNotFoundError(f"controle real ausente em {controles}")

    print(f"== {args.submissao.name} contexto {args.contexto}", flush=True)
    laudo, por_alvo = auditar(
        args.submissao, controles, args.contexto, args.n_alvos, args.semente, args.threads
    )
    print(
        f"  genes: {laudo.genes_no_eixo} no eixo, {laudo.genes_apos_o_gate} apos o gate de "
        f"CPM > {GATE_CPM}"
    )
    print(
        f"  nsig por alvo: mediana {laudo.nsig_mediano:.0f}  "
        f"[p05 {laudo.nsig_p05:.0f}, p25 {laudo.nsig_p25:.0f}, p75 {laudo.nsig_p75:.0f}, "
        f"p95 {laudo.nsig_p95:.0f}]  min {laudo.nsig_min}  max {laudo.nsig_max}"
    )
    print(
        f"  contra os {DE_REAIS} reais: {laudo.razao_para_de_real:.2f}x  |  "
        f"{laudo.alvos_sem_nenhum_gene} alvos sem nenhum gene significativo"
    )
    print(
        f"  duas contagens concordam: {laudo.duas_contagens_concordam} "
        f"(maior divergencia {laudo.maior_divergencia})"
    )

    if args.relatorio:
        args.relatorio.parent.mkdir(parents=True, exist_ok=True)
        args.relatorio.write_text(
            json.dumps(
                {
                    "submissao": str(args.submissao),
                    "controles": str(controles),
                    "parametros": {
                        "p_adj_threshold": P_ADJ,
                        "min_abs_log2fc": MIN_ABS_LOG2FC,
                        "method": METODO,
                        "nan_lfc_policy": NAN_LFC,
                        "mean_calc": MEAN_CALC,
                        "fdr_scope": FDR_SCOPE,
                        "filter_gene_min_cpm_cell": GATE_CPM,
                        "target_sum": TARGET_SUM,
                        "backend": "scanpy",
                    },
                    "chamadas_do_cell_eval2": [
                        "de_compute.compute_de",
                        "de.prep_de_side",
                        "de.rank_de_side",
                    ],
                    "por_que_amostrar_e_exato": (
                        "fdr_scope=per_pert: de._bh_per_target refaz o BH dentro de cada alvo "
                        "sobre os p-valores dele, entao o nsig de um alvo nao depende de quais "
                        "outros alvos entraram na corrida."
                    ),
                    "por_que_nao_ha_lado_retido": (
                        "control_source=real: o lado de controle do teste de DE vem do dado "
                        "real, e os controles de D/E/F estao publicados em controls-final."
                    ),
                    "de_reais_de_referencia": DE_REAIS,
                    "laudo": asdict(laudo),
                    "nsig_por_alvo": por_alvo,
                },
                indent=2,
                ensure_ascii=False,
            ),
            encoding="utf-8",
        )
        print(f"\nrelatorio em {args.relatorio}")


if __name__ == "__main__":
    main()
