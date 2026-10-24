"""
Etapa 2 de docs/VCC-2026-ROTA.md: `ctrl_bulk` — a linha de base de que toda
predicao e um delta, mais o universo de genes que a pontuacao chega a testar.

Quatro produtos, cada um consumido por uma etapa adiante:

`b_ctrl`        log1p(5e4 * P / sum(P)), P = pseudobulk de soma dos controles.
                E a origem do eixo: o `delta` que o modelo preve se soma aqui.
                CP50k e log1p porque e o que o cell-eval2 aplica antes das
                metricas de magnitude.

`se_jackknife`  erro padrao de `b_ctrl` por gene, por jackknife delete-1 sobre
                as 18.400 celulas. Gene cujo proprio controle e ruidoso nao
                sustenta chamada de DE, e e onde a predicao deve recuar para
                zero em vez de arriscar.

`universo_de`   genes com CPM medio no controle acima de 5. O Wilcoxon do
                cell-eval2 so testa esses, entao predizer efeito fora dali nao
                rende e predizer dentro e obrigatorio. Duas definicoes de "CPM
                medio" circulam e dao universos diferentes; as duas sao
                calculadas e o resumo diz qual reproduz o numero medido.

`biblioteca`    distribuicao de contagem total e de genes detectados por
                celula. O `emissor_400` tem de devolver contagens cruas com
                essa forma, nao uma matriz normalizada.

O jackknife delete-1 tem forma fechada sobre uma soma: P sem a celula i e
P - x_i, e o total sem i e T - t_i. Entao nao ha 18.400 recalculos — ha um
produto externo por blocos de genes, corrigido nas entradas nao nulas.

Saidas em runs/vcc/ctrl_bulk/. Rodar da raiz do projeto:
    ./.venv/Scripts/python.exe -I scripts/vcc/ctrl_bulk.py
"""

from __future__ import annotations

import json
from dataclasses import dataclass
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import scipy.sparse as sp

VCC = Path("data/raw/vcc/controls-final")
SAIDA = Path("runs/vcc/ctrl_bulk")
CONTEXTOS = ("D", "E", "F")

# O cell-eval2 normaliza o pseudobulk para 50.000 contagens antes do log1p.
ESCALA_PSEUDOBULK = 5e4
# Corte do universo de DE, em CPM medio no controle.
CPM_MINIMO = 5.0
# Blocos de genes no produto externo do jackknife: 512 x 18.400 float64 = 75 MB.
BLOCO_GENES = 512


@dataclass
class Controle:
    contexto: str
    n_celulas: int
    total_umi: float
    pseudobulk: np.ndarray
    b_ctrl: np.ndarray
    se_jackknife: np.ndarray
    cpm_pseudobulk: np.ndarray
    cpm_medio_celula: np.ndarray
    umi_por_celula: np.ndarray
    genes_por_celula: np.ndarray


def _vetor(m: object) -> np.ndarray:
    """Soma de matriz esparsa vem como np.matrix 2D; achata para 1D float64."""
    return np.asarray(m, dtype=np.float64).ravel()


def jackknife_delete1(
    indptr: np.ndarray,
    indices: np.ndarray,
    dados: np.ndarray,
    pseudobulk: np.ndarray,
    umi_por_celula: np.ndarray,
    total_umi: float,
) -> np.ndarray:
    """Erro padrao de `b_ctrl` por gene, deixando uma celula de fora por vez.

    Recebe os internos de um CSC (coluna = gene) em vez da matriz: e tudo que a
    conta usa, e evita que o tipo exato de esparsa do scipy entre na assinatura.

    Forma fechada: b sem a celula i vale log1p(ESCALA * (P_g - x_ig) / (T - t_i)).
    Para as entradas nulas (a esmagadora maioria) o numerador nao muda, so o
    denominador — isso permite varrer todas as celulas com um produto externo
    por bloco de genes e depois corrigir apenas as entradas nao nulas.
    """
    n_celulas = int(umi_por_celula.size)
    n_genes = int(indptr.size) - 1
    # u_i = ESCALA / (T - t_i): a parte que nao depende do gene.
    u = ESCALA_PSEUDOBULK / (total_umi - umi_por_celula)
    variancia = np.empty(n_genes, dtype=np.float64)

    for inicio in range(0, n_genes, BLOCO_GENES):
        fim = min(inicio + BLOCO_GENES, n_genes)
        p_bloco = pseudobulk[inicio:fim, None]
        # Valor que cada celula daria se sua contagem deste gene fosse zero.
        valores = np.log1p(p_bloco * u[None, :])
        soma = valores.sum(axis=1)
        soma_quadrados = np.square(valores).sum(axis=1)

        for deslocamento, gene in enumerate(range(inicio, fim)):
            ini, fin = indptr[gene], indptr[gene + 1]
            if ini == fin:
                continue
            celulas = indices[ini:fin]
            contagens = dados[ini:fin]
            # Troca o valor assumido (zero) pelo real nas celulas nao nulas.
            assumido = np.log1p(pseudobulk[gene] * u[celulas])
            real = np.log1p((pseudobulk[gene] - contagens) * u[celulas])
            soma[deslocamento] += real.sum() - assumido.sum()
            soma_quadrados[deslocamento] += np.square(real).sum() - np.square(assumido).sum()

        media = soma / n_celulas
        # Variancia do jackknife: (n-1)/n * soma dos desvios ao quadrado.
        desvio = soma_quadrados / n_celulas - np.square(media)
        variancia[inicio:fim] = (n_celulas - 1) * np.maximum(desvio, 0.0)

    return np.sqrt(variancia)


def medir(contexto: str, n_genes_esperado: int) -> Controle:
    adata = ad.read_h5ad(VCC / f"context_{contexto}.h5ad")
    if int(adata.n_vars) != n_genes_esperado:
        raise ValueError(
            f"contexto {contexto}: {adata.n_vars} genes, esperado {n_genes_esperado}"
        )
    x = sp.csr_matrix(adata.X)
    x.eliminate_zeros()
    if x.data.size and not np.allclose(x.data, np.round(x.data)):
        raise ValueError(f"contexto {contexto}: X nao e contagem inteira")

    n_celulas = int(adata.n_obs)
    umi_por_celula = _vetor(x.sum(axis=1))
    genes_por_celula = _vetor(x.getnnz(axis=1))
    if float(umi_por_celula.min()) <= 0.0:
        raise ValueError(f"contexto {contexto}: celula com zero UMI quebra a normalizacao")

    pseudobulk = _vetor(x.sum(axis=0))
    total_umi = float(pseudobulk.sum())
    b_ctrl = np.log1p(ESCALA_PSEUDOBULK * pseudobulk / total_umi)
    cpm_pseudobulk = 1e6 * pseudobulk / total_umi

    # Definicao alternativa: media, entre celulas, do CPM da propria celula.
    escalada = sp.csr_matrix(sp.diags(1e6 / umi_por_celula) @ x)
    cpm_medio_celula = _vetor(escalada.mean(axis=0))

    por_gene_csc = sp.csc_matrix(x)
    se = jackknife_delete1(
        np.asarray(por_gene_csc.indptr),
        np.asarray(por_gene_csc.indices),
        np.asarray(por_gene_csc.data, dtype=np.float64),
        pseudobulk,
        umi_por_celula,
        total_umi,
    )

    print(
        f"{contexto}: {n_celulas} celulas, {total_umi:.3e} UMI, "
        f"mediana {np.median(umi_por_celula):.0f} UMI/celula, "
        f"{np.median(genes_por_celula):.0f} genes/celula, "
        f"SE do b_ctrl mediana {np.median(se):.2e}"
    )
    return Controle(
        contexto=contexto,
        n_celulas=n_celulas,
        total_umi=total_umi,
        pseudobulk=pseudobulk,
        b_ctrl=b_ctrl,
        se_jackknife=se,
        cpm_pseudobulk=cpm_pseudobulk,
        cpm_medio_celula=cpm_medio_celula,
        umi_por_celula=umi_por_celula,
        genes_por_celula=genes_por_celula,
    )


def main() -> None:
    SAIDA.mkdir(parents=True, exist_ok=True)
    genes = pd.read_csv(VCC / "gene_names.csv")["gene_name"].astype(str).tolist()
    medidas = {ctx: medir(ctx, len(genes)) for ctx in CONTEXTOS}

    # Numeros publicados no memorando, para dizer qual definicao de CPM os gerou.
    alvo_memorando = {"D": 9633, "E": 9743, "F": 9537}
    universo: list[dict[str, object]] = []
    conferencia: dict[str, dict[str, int | bool]] = {}
    for ctx, m in medidas.items():
        por_pseudobulk = int((m.cpm_pseudobulk > CPM_MINIMO).sum())
        por_media_celula = int((m.cpm_medio_celula > CPM_MINIMO).sum())
        conferencia[ctx] = {
            "cpm_do_pseudobulk": por_pseudobulk,
            "cpm_medio_por_celula": por_media_celula,
            "memorando": alvo_memorando[ctx],
            "pseudobulk_reproduz": por_pseudobulk == alvo_memorando[ctx],
            "media_celula_reproduz": por_media_celula == alvo_memorando[ctx],
        }
        print(
            f"{ctx}: universo DE por pseudobulk {por_pseudobulk}, "
            f"por media de celula {por_media_celula}, memorando {alvo_memorando[ctx]}"
        )
        universo.append(
            {
                "contexto": ctx,
                "n_por_pseudobulk": por_pseudobulk,
                "n_por_media_celula": por_media_celula,
            }
        )

    # Uma linha por gene e contexto: o que a etapa 3 e a 6 leem.
    por_gene = pd.DataFrame({"gene": genes})
    for ctx, m in medidas.items():
        por_gene[f"b_ctrl_{ctx}"] = m.b_ctrl
        por_gene[f"se_{ctx}"] = m.se_jackknife
        por_gene[f"cpm_pb_{ctx}"] = m.cpm_pseudobulk
        por_gene[f"cpm_cel_{ctx}"] = m.cpm_medio_celula
        # Media do CPM por celula, NAO CPM do pseudobulk: medido, e a unica das
        # duas que reproduz 9633/9743/9537. A diferenca e de 68 genes em D e 38
        # em F, e cada um deles seria predito-sem-teste ou testado-sem-predicao,
        # o que bate direto em jaccard e reach.
        por_gene[f"no_universo_{ctx}"] = m.cpm_medio_celula > CPM_MINIMO
    por_gene["no_universo_em_todos"] = por_gene[
        [f"no_universo_{c}" for c in CONTEXTOS]
    ].all(axis=1)
    por_gene["no_universo_em_algum"] = por_gene[
        [f"no_universo_{c}" for c in CONTEXTOS]
    ].any(axis=1)
    por_gene.to_csv(SAIDA / "por_gene.tsv", sep="\t", index=False)

    biblioteca = pd.DataFrame(
        [
            {
                "contexto": ctx,
                "celulas": m.n_celulas,
                "umi_total": m.total_umi,
                "umi_p01": float(np.percentile(m.umi_por_celula, 1)),
                "umi_p25": float(np.percentile(m.umi_por_celula, 25)),
                "umi_mediana": float(np.median(m.umi_por_celula)),
                "umi_p75": float(np.percentile(m.umi_por_celula, 75)),
                "umi_p99": float(np.percentile(m.umi_por_celula, 99)),
                "umi_max": float(m.umi_por_celula.max()),
                "genes_mediana": float(np.median(m.genes_por_celula)),
                "nnz_por_celula_medio": float(m.genes_por_celula.mean()),
                "pseudobulk_max_gene": float(m.pseudobulk.max()),
            }
            for ctx, m in medidas.items()
        ]
    )
    biblioteca.to_csv(SAIDA / "biblioteca.tsv", sep="\t", index=False)
    print("\n== biblioteca por contexto ==")
    print(biblioteca.to_string(index=False))

    # Arrays binarios: o emissor_400 precisa das distribuicoes empiricas.
    np.savez_compressed(
        SAIDA / "ctrl_bulk.npz",
        genes=np.array(genes, dtype=object),
        **{f"b_ctrl_{c}": medidas[c].b_ctrl for c in CONTEXTOS},
        **{f"se_{c}": medidas[c].se_jackknife for c in CONTEXTOS},
        **{f"pseudobulk_{c}": medidas[c].pseudobulk for c in CONTEXTOS},
        **{f"umi_por_celula_{c}": medidas[c].umi_por_celula for c in CONTEXTOS},
        **{f"genes_por_celula_{c}": medidas[c].genes_por_celula for c in CONTEXTOS},
        allow_pickle=True,
    )

    resumo = {
        "escala_pseudobulk": ESCALA_PSEUDOBULK,
        "cpm_minimo": CPM_MINIMO,
        "contextos": {
            ctx: {
                "celulas": m.n_celulas,
                "umi_total": m.total_umi,
                "umi_mediana_por_celula": float(np.median(m.umi_por_celula)),
                "se_b_ctrl_mediana": float(np.median(m.se_jackknife)),
                "se_b_ctrl_p99": float(np.percentile(m.se_jackknife, 99)),
            }
            for ctx, m in medidas.items()
        },
        "universo_de": universo,
        "conferencia_contra_memorando": conferencia,
        "universo_usado": "cpm_medio_por_celula",
        "por_que": (
            "Das duas definicoes de 'CPM medio no controle', so a media do CPM por "
            "celula reproduz os tres numeros medidos (9633/9743/9537). A do pseudobulk "
            "da 9565/9744/9499. O cell-eval2 usa a primeira."
        ),
        "genes_no_universo_em_todos": int(por_gene["no_universo_em_todos"].sum()),
        "genes_no_universo_em_algum": int(por_gene["no_universo_em_algum"].sum()),
        "ressalva": (
            "Isto e exato se os 18.400 controles entregues forem os de pontuacao. "
            "manifest.json declara 138.400 celulas de verdade (300x400 + 18.400) e "
            "panel_id vcc2026-test-1, o que e indicio forte, mas a especificacao fala "
            "em 'model inputs only'. Se diferirem, b_ctrl desloca e tudo a jusante com ele."
        ),
    }
    (SAIDA / "resumo.json").write_text(
        json.dumps(resumo, indent=2, ensure_ascii=False), encoding="utf-8"
    )
    print(f"\nescrito em {SAIDA}")


if __name__ == "__main__":
    main()
