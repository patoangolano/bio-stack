"""
Caracterizacao dos contextos D/E/F do VCC 2026 a partir das celulas controle.

Tres perguntas, todas respondiveis sem uma unica celula perturbada:

1. Que tipo celular e cada contexto? Paineis de marcadores, exigindo dois ou
   mais concordantes antes de aceitar um rotulo. O CellTypist e treinado em
   atlas imune e nao serve para tipo nao imune; por isso paineis explicitos.
2. Como sao as bibliotecas? Contagens por celula, genes detectados, %MT e
   %ribo. A submissao e de contagens CRUAS, entao o gerador de predicoes tem
   de reproduzir essas distribuicoes. Nada e filtrado: estas celulas sao a
   distribuicao oficial do desafio, nao dado a limpar.
3. Cada um dos 300 alvos e expresso no contexto? CRISPRi reprime transcricao;
   alvo que o contexto nao transcreve nao tem o que derrubar, e a perturbacao
   tende ao controle.

Os marcadores de contexto sao um RANQUEAMENTO por especificidade, nao um
teste: celula nao e replica (ver skill rigor-estatistico, item 6), e um valor
de p calculado com 18 400 celulas por grupo seria inflado por ordens de
grandeza. Nenhum p-valor e produzido aqui.

Saidas em runs/vcc/caracterizacao/.
"""

from __future__ import annotations

import json
from dataclasses import dataclass
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import scipy.sparse as sp

BASE = Path("data/raw/vcc/controls-final")
SAIDA = Path("runs/vcc/caracterizacao")
CONTEXTOS = ("D", "E", "F")

# Contextos plausiveis para um desafio de Perturb-seq: as linhagens dos grandes
# conjuntos publicos (K562, RPE1, HepG2, Jurkat, H1) e os destinos de
# diferenciacao de celula-tronco. Cada painel tem marcadores que outro painel
# nao compartilha, para que dois concordantes signifiquem alguma coisa.
PAINEIS: dict[str, list[str]] = {
    "pluripotente (H1/iPSC)": ["POU5F1", "NANOG", "LIN28A", "DPPA4", "TDGF1", "L1TD1"],
    "eritroleucemia (K562)": ["HBG1", "HBG2", "GATA1", "TAL1", "GYPA", "HBZ"],
    "hepatocito (HepG2)": ["ALB", "APOA1", "APOA2", "AFP", "TTR", "FGB"],
    "linfocito T (Jurkat)": ["CD3E", "CD3D", "LCK", "CD2", "CD247"],
    "epitelio pigmentar (RPE1)": ["KRT8", "KRT18", "VIM", "TGFBI", "CAV1"],
    "progenitor neural": ["PAX6", "NES", "SOX1", "HES5", "FABP7"],
    "neuronio": ["STMN2", "DCX", "MAP2", "ELAVL4", "SYT1"],
    "endoderma definitivo": ["SOX17", "FOXA2", "GATA6", "CXCR4", "CER1"],
    "mesoderma": ["TBXT", "MIXL1", "MESP1", "HAND1", "PDGFRA"],
    "cardiomiocito": ["TNNT2", "MYL7", "NKX2-5", "ACTC1", "MYH6"],
    "fibroblasto/mesenquima": ["COL1A1", "COL1A2", "COL3A1", "DCN", "LUM"],
    "endotelio": ["PECAM1", "CDH5", "KDR", "TEK", "VWF"],
    "epitelio (carcinoma)": ["EPCAM", "KRT19", "CDH1", "CLDN4", "ELF3"],
    "mieloide": ["LYZ", "CD14", "CSF1R", "TYROBP", "FCER1G"],
    "linfocito B": ["CD79A", "CD79B", "MS4A1", "CD19", "PAX5"],
}

# Para um marcador contar como "concordante" ele precisa aparecer numa fracao
# razoavel das celulas, nao so ter media alta puxada por poucas.
FRACAO_MINIMA = 0.30


@dataclass
class Contexto:
    qc: dict[str, float | int | str]
    paineis: list[dict[str, float | int | str]]
    alvos: list[dict[str, float | bool | str]]
    media_log: np.ndarray
    fracao: np.ndarray


def _vetor(m: object) -> np.ndarray:
    """Soma de matriz esparsa vem como np.matrix 2D; achata para 1D."""
    return np.asarray(m, dtype=np.float64).ravel()


def caracterizar(ctx: str, genes: pd.Index, alvos: list[str]) -> Contexto:
    adata = ad.read_h5ad(BASE / f"context_{ctx}.h5ad")
    n_celulas = int(adata.n_obs)
    x = sp.csr_matrix(adata.X)
    x.eliminate_zeros()
    if x.data.size and not np.allclose(x.data, np.round(x.data)):
        raise ValueError(f"contexto {ctx}: X nao e contagem inteira")

    total = _vetor(x.sum(axis=1))
    detectados = _vetor(x.getnnz(axis=1))
    mt = np.asarray(genes.str.startswith("MT-"), dtype=bool)
    ribo = np.asarray(genes.str.match(r"^RP[SL]\d"), dtype=bool)
    pct_mt = 100 * _vetor(x[:, np.flatnonzero(mt)].sum(axis=1)) / total
    pct_ribo = 100 * _vetor(x[:, np.flatnonzero(ribo)].sum(axis=1)) / total

    # log1p de CP10k: a normalizacao que paineis de marcadores pressupoem.
    # O CSR e reconstruido em vez de mutado no lugar.
    escalada = sp.csr_matrix(sp.diags(1e4 / total) @ x)
    log_norm = sp.csr_matrix((np.log1p(escalada.data), escalada.indices, escalada.indptr), shape=escalada.shape)
    media_log = _vetor(log_norm.mean(axis=0))
    fracao = _vetor(x.getnnz(axis=0)) / n_celulas
    media_crua = _vetor(x.mean(axis=0))

    qc: dict[str, float | int | str] = {
        "contexto": ctx,
        "celulas": n_celulas,
        "contagens_mediana": float(np.median(total)),
        "contagens_p10": float(np.percentile(total, 10)),
        "contagens_p90": float(np.percentile(total, 90)),
        "genes_detectados_mediana": float(np.median(detectados)),
        "pct_mt_mediana": round(float(np.median(pct_mt)), 2),
        "pct_ribo_mediana": round(float(np.median(pct_ribo)), 2),
        "genes_mt": int(mt.sum()),
        "max_contagem": float(x.data.max()) if x.data.size else 0.0,
    }

    indice = {g: i for i, g in enumerate(genes)}
    paineis: list[dict[str, float | int | str]] = []
    for nome, marcadores in PAINEIS.items():
        presentes = [m for m in marcadores if m in indice]
        idx = [indice[m] for m in presentes]
        concordantes = [m for m in presentes if fracao[indice[m]] >= FRACAO_MINIMA]
        paineis.append(
            {
                "contexto": ctx,
                "painel": nome,
                "escore_medio_log": round(float(media_log[idx].mean()), 3) if idx else float("nan"),
                "marcadores_presentes": len(presentes),
                "concordantes": len(concordantes),
                "quais": ",".join(concordantes),
            }
        )

    linhas_alvos: list[dict[str, float | bool | str]] = []
    for alvo in alvos:
        i = indice.get(alvo)
        linhas_alvos.append(
            {
                "contexto": ctx,
                "alvo": alvo,
                "na_matriz": i is not None,
                "media_crua": float(media_crua[i]) if i is not None else float("nan"),
                "media_log": float(media_log[i]) if i is not None else float("nan"),
                "fracao_expressa": float(fracao[i]) if i is not None else float("nan"),
            }
        )
    return Contexto(qc=qc, paineis=paineis, alvos=linhas_alvos, media_log=media_log, fracao=fracao)


def main() -> None:
    SAIDA.mkdir(parents=True, exist_ok=True)
    genes = pd.Index(pd.read_csv(BASE / "gene_names.csv")["gene_name"].astype(str))
    alvos = pd.read_csv(BASE / "pert_counts.csv")["target_gene"].astype(str).tolist()

    resultados = {ctx: caracterizar(ctx, genes, alvos) for ctx in CONTEXTOS}

    qc_df = pd.DataFrame([r.qc for r in resultados.values()])
    qc_df.to_csv(SAIDA / "qc_por_contexto.tsv", sep="\t", index=False)
    paineis = pd.DataFrame([p for r in resultados.values() for p in r.paineis])
    paineis.to_csv(SAIDA / "paineis.tsv", sep="\t", index=False)
    alvos_df = pd.DataFrame([a for r in resultados.values() for a in r.alvos])
    alvos_df.to_csv(SAIDA / "alvos_por_contexto.tsv", sep="\t", index=False)

    # Especificidade de contexto: media log no contexto menos a maior media
    # entre os outros dois, so entre genes expressos em >= 30% das celulas.
    medias = pd.DataFrame({ctx: resultados[ctx].media_log for ctx in CONTEXTOS}, index=genes)
    fracoes = pd.DataFrame({ctx: resultados[ctx].fracao for ctx in CONTEXTOS}, index=genes)
    marcadores: dict[str, list[str]] = {}
    for ctx in CONTEXTOS:
        outros = [c for c in CONTEXTOS if c != ctx]
        especificidade = medias[ctx] - medias[outros].max(axis=1)
        elegiveis = especificidade[fracoes[ctx] >= FRACAO_MINIMA]
        marcadores[ctx] = [str(g) for g in elegiveis.sort_values(ascending=False).head(25).index]
    expressos = fracoes.max(axis=1) >= 0.10
    correlacao = medias.loc[expressos].corr(method="spearman").round(4)

    resumo = {
        "genes": len(genes),
        "alvos": len(alvos),
        "alvos_fora_da_matriz": sorted(set(alvos) - set(genes)),
        "genes_para_correlacao": int(expressos.sum()),
        "spearman_pseudobulk": correlacao.to_dict(),
        "marcadores_especificos": marcadores,
    }
    (SAIDA / "resumo.json").write_text(json.dumps(resumo, indent=2, ensure_ascii=False), encoding="utf-8")

    print("== QC por contexto ==")
    print(qc_df.to_string(index=False))
    print("\n== Paineis com >= 2 marcadores concordantes ==")
    fortes = paineis.loc[paineis["concordantes"] >= 2].sort_values(
        by=["contexto", "escore_medio_log"], ascending=[True, False]
    )
    print(fortes[["contexto", "painel", "escore_medio_log", "concordantes", "quais"]].to_string(index=False))
    print(f"\n== Alvos fora da matriz: {resumo['alvos_fora_da_matriz'] or 'nenhum'}")
    print(f"\n== Spearman entre contextos (pseudobulk, {resumo['genes_para_correlacao']} genes) ==")
    print(correlacao.to_string())
    print("\n== Marcadores mais especificos (top 25) ==")
    for ctx, lista in marcadores.items():
        print(f"  {ctx}: {', '.join(lista)}")


if __name__ == "__main__":
    main()
