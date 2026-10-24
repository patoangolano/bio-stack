"""
Etapas 3 e 4 de docs/VCC-2026-ROTA.md: `transfer_fit` — o modelo proprio de
transferencia de resposta em pseudobulk — e a validacao "deixa um par de fora"
que diz se ele generaliza.

As duas etapas moram no mesmo arquivo porque a segunda mede a primeira: separar
duplicaria o carregamento dos arrays, que e a parte caraa. Subcomandos
`ajustar` e `validar`.

## A dificuldade central, dita em voz alta

Nao existe uma unica celula perturbada em D, E ou F. Entao o coeficiente de
transferencia doador -> alvo para os nossos contextos e **inestimavel por
construcao**. O que se pode estimar e o coeficiente TIPICO entre pares de
contexto que existem (hct116 <-> hek293t no Orion, 6 pares no Replogle, pares
do Jiang), e usar isso como prior, modulado pelo unico dado que temos de D/E/F:
a expressao de controle. Gene que o contexto nao transcreve nao se move nele.

A validacao da etapa 4 e exatamente a pergunta honesta: um coeficiente
estimado num par serve para um par que ficou de fora? Se nao servir, o modelo
nao tem direito de ser aplicado a D/E/F e a submissao deve recuar para a
resposta media.

## Por que o alvo e o centrado, e nao a resposta bruta

O zero da pontuacao do vcc2026 e a RESPOSTA MEDIA REAL do contexto as 300
perturbacoes, nao o controle. Predizer a media do contexto pontua exatamente
zero. Logo o sinal que vale e o desvio especifico da perturbacao em relacao a
essa media, e o modelo e ajustado sobre o centrado:

    y~_c(P,g) = y_c(P,g) - media_P' y_c(P',g)

## Por que a variavel transferida e uma razao, e nao um delta em log

`ctrl_means` do PIE e a media, entre celulas, de log1p(CP10k). O nosso `b_ctrl`
e log1p do CP50k do pseudobulk de soma. Sao estimadores diferentes:
sum(expm1(ctrl_means)) da 8.615, nao 10.000, que e a desigualdade de Jensen
aparecendo. Transferir um delta em log arrastaria essa divergencia para dentro
da predicao.

A saida e transferir a RAZAO LINEAR implicada, que e adimensional, e so depois
aplica-la ao nosso proprio baseline:

    r = expm1(ctrl_means + delta_p) / expm1(ctrl_means)
    delta_nosso = log1p(expm1(b_ctrl) * r * renorm) - b_ctrl

Medido: `delta_p ~ log1p(expm1(ctrl_means)*fold_changes) - ctrl_means` com
r=0,963 e erro absoluto mediano 0,003. Usamos `delta_p` e nao `fold_changes`
porque `delta_p` e denso (0,4% de zeros) enquanto `fold_changes` e zero
exatamente onde `tested` e False (51,2% testado), e densidade aqui e dado.

Rodar da raiz do projeto:
    ./.venv/Scripts/python.exe -I scripts/vcc/transfer_fit.py validar
    ./.venv/Scripts/python.exe -I scripts/vcc/transfer_fit.py ajustar
"""

from __future__ import annotations

import argparse
import json
from dataclasses import dataclass
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

ARRAYS = Path("data/raw/pie/arrays")
METAS = Path("data/raw/pie")


def caminho_meta(repo_dir: str) -> Path:
    """`meta.json` de um repositorio PIE, no layout do `baixar_pie.py`.

    Na execucao de 22/10 havia uma copia em `data/raw/pie/<repo>/meta.json`;
    ela se perdeu no re-clone de 23/10 e o `baixar_pie.py` traz so a de
    `preprocessed/`, que e o mesmo arquivo do repositorio no Hub.
    """
    return ARRAYS / repo_dir / "preprocessed" / "meta.json"
VCC = Path("data/raw/vcc/controls-final")
EIXO = Path("runs/vcc/eixo")
CTRL = Path("runs/vcc/ctrl_bulk")
SAIDA = Path("runs/vcc/transfer")

# log2 da razao fora desta faixa e artefato de gene quase ausente no controle,
# nao biologia: 2^5 = 32 vezes para cada lado ja e generoso.
LIMITE_LOG2R = 5.0
# Numero minimo de perturbacoes com medida nos dois contextos para que o
# coeficiente por gene tenha direito a existir.
MIN_PERTS_POR_GENE = 100
# Forca do encolhimento do coeficiente por gene para o global, em unidades de
# "perturbacoes equivalentes". Com 3.277 perts, 50 e um encolhimento leve.
K_ENCOLHIMENTO = 50.0
# Minimo de perturbacoes em comum para um par entrar na validacao. 60 e baixo
# para coeficiente por gene — e de proposito: com menos de MIN_PERTS_POR_GENE
# perturbacoes nenhum gene ganha coeficiente proprio e tudo encolhe para o
# global, que e o comportamento correto num par escasso.
MINIMO_PERTS_PAR = 60
# Minimo de perturbacoes para um contexto entrar no AGRUPAMENTO DA MEDIA. A
# media de um contexto com poucas perturbacoes e ruido: no limite de uma
# perturbacao ela E aquela perturbacao. Os 30 contextos do Jiang tem de 1 a 39
# perturbacoes e ficam todos fora; os 8 restantes tem de 44 a 8.987.
#
# O valor cai UMA perturbacao acima do maior contexto do Jiang (39), o que e
# perto demais para passar sem teste. Medido em pool_media_min10.json (30
# contextos) e pool_media_min200.json (6): `atual` perde nos tres cortes e
# `concordancia` nao vence em nenhum. Em 200 `profundidade` vence, mas esse
# corte exclui o arc_vcc_25 — a propria patologia que torna profundidade ruim —
# entao nao e evidencia a favor dela.
MIN_PERTS_MEDIA = 40
# Programas de resposta retidos na base de baixa dimensao. Com ~1.000
# perturbacoes por par, 32 componentes ficam bem determinados; o modelo passa
# de 18.533 parametros (um por gene) para 32 (um por programa).
K_PROGRAMAS = 32
# Encolhimento da resposta media pela cobertura do gene, em entradas validas
# equivalentes. Com 20, um gene medido em 20 perturbacoes fica com metade do
# peso, um medido em 1 fica com 5%, e um de cobertura cheia (1.253 no Replogle)
# fica praticamente intacto.
K_MEDIA = 20.0
# Fracao minima de perturbacoes em que o gene precisa ter sido testado para a
# media dele entrar na estimativa da resposta do contexto.
#
# ESTE E O PARAMETRO MAIS IMPORTANTE DESTE ARQUIVO, e ele existe porque o
# encolhimento por CONTAGEM (K_MEDIA) nao resolve o problema que parecia
# resolver. A media de log2(fold_changes) mascarada por `tested` e selecionada
# no proprio desfecho: um gene pouco expresso so e "testado" nas perturbacoes
# em que sua contagem saiu alta, entao a razao sai > 1 por construcao.
#
# Medido, em orion.hct116 (8.987 perturbacoes):
#
#   fracao testada   genes   media bruta   % > 0
#   0,00 - 0,01      1.742      +3,358     100,0%
#   0,01 - 0,05        503      +1,539     100,0%
#   0,05 - 0,20        417      +0,895     100,0%
#   0,20 - 0,50        346      +0,492     100,0%
#   0,90 - 1,00      9.871      -0,029       6,2%
#
# Quatro faixas seguidas a 100% positivo nao e biologia. E a faixa >= 0,90 e
# praticamente nao enviesada, o que fixa o corte aqui e nao por gosto.
#
# A correlacao entre fracao testada e media bruta e -0,878 no orion.hct116 e
# -0,831 no arc_vcc_25; no Replogle, que testa 100% dos genes, e +0,040. O
# Replogle e o controle que fecha o diagnostico.
#
# Por que a contagem nao bastava: com 8.987 perturbacoes, um gene testado em 1%
# tem n=90, peso 90/110 = 0,818, e sobra +2,75 de uma media que deveria ser ~0.
# O vies escala com a FRACAO; o encolhimento olhava a contagem.
FRACAO_MINIMA_COBERTURA = 0.9

# Pares de contexto disponiveis para estimar e validar transferencia. O Replogle
# nao tem nenhum dos 300 alvos, e por isso mesmo e o melhor conjunto para medir
# transferencia generica: nao ha como se enganar com vazamento de alvo.
REPOS = {
    "orion": "PIE_x_atlas_orion",
    "replogle": "PIE_replogle_nadig_essential",
    "arc_vcc_25": "PIE_arc_vcc_25",
    "jiang": "PIE_jiang",
}


@dataclass
class Resposta:
    """Resposta centrada de um contexto, nas perturbacoes e genes comuns."""

    repo: str
    contexto: str
    perts: list[str]
    genes: list[str]
    log2r: np.ndarray  # (n_perts, n_genes) centrado por gene
    valido: np.ndarray  # (n_perts, n_genes) bool: entrada utilizavel
    media_por_gene: np.ndarray  # (n_genes,) a media que foi removida
    fdr: np.ndarray  # (n_perts, n_genes)


def _carregar(repo_dir: str, nome: str, linhas: np.ndarray, colunas: np.ndarray) -> np.ndarray:
    """Le um array do PIE por memoria mapeada e recorta linhas e colunas.

    O recorte de coluna acontece por bloco de linhas: fancy indexing nos dois
    eixos de uma vez materializaria a matriz inteira antes de cortar.
    """
    caminho = ARRAYS / repo_dir / "preprocessed" / f"{nome}.npy"
    fonte = np.load(caminho, mmap_mode="r", allow_pickle=False)
    saida = np.empty((linhas.size, colunas.size), dtype=fonte.dtype)
    passo = 512
    for inicio in range(0, linhas.size, passo):
        fim = min(inicio + passo, linhas.size)
        saida[inicio:fim] = np.asarray(fonte[linhas[inicio:fim]])[:, colunas]
    return saida


def media_do_contexto(repo: str, contexto: str, genes_alvo: list[str]) -> np.ndarray:
    """Media por gene de log2(fold_changes) sobre TODAS as perturbacoes do contexto.

    Centrar a resposta do doador sobre apenas o nosso subconjunto de alvos
    degenera quando o subconjunto e pequeno: com uma perturbacao so, o centrado
    e exatamente zero e o contexto contribui nada enquanto consome peso. A media
    do contexto e uma propriedade do contexto, e tem de ser estimada com tudo
    que ele mediu.

    Varre por blocos de linha e acumula soma e contagem: nao materializa a
    matriz, que no Orion hct116 teria 8.987 x 14.690 entradas.
    """
    repo_dir = REPOS[repo]
    meta = json.loads((caminho_meta(repo_dir)).read_text(encoding="utf-8"))
    id_contexto = int(meta["context_to_id"][contexto])
    eixo_repo = {str(g): i for i, g in enumerate(meta["genes"])}
    colunas = np.array([eixo_repo[g] for g in genes_alvo], dtype=np.int64)

    ctx_ids = np.load(ARRAYS / repo_dir / "preprocessed" / "ctx_ids.npy", allow_pickle=False)
    linhas = np.flatnonzero(ctx_ids == id_contexto)
    fonte_fc = np.load(
        ARRAYS / repo_dir / "preprocessed" / "fold_changes.npy", mmap_mode="r", allow_pickle=False
    )
    fonte_te = np.load(
        ARRAYS / repo_dir / "preprocessed" / "tested.npy", mmap_mode="r", allow_pickle=False
    )
    soma = np.zeros(colunas.size, dtype=np.float64)
    contagem = np.zeros(colunas.size, dtype=np.int64)
    passo = 512
    for inicio in range(0, linhas.size, passo):
        fim = min(inicio + passo, linhas.size)
        bloco = linhas[inicio:fim]
        fc = np.asarray(fonte_fc[bloco])[:, colunas].astype(np.float64)
        te = np.asarray(fonte_te[bloco])[:, colunas]
        with np.errstate(divide="ignore", invalid="ignore"):
            valor = np.log2(fc)
        ok = te & (fc > 0) & np.isfinite(valor) & (np.abs(valor) <= LIMITE_LOG2R)
        soma += np.where(ok, valor, 0.0).sum(axis=0)
        contagem += ok.sum(axis=0)
    # Encolhimento pela cobertura, e nao e detalhe. `tested` cobre 51% dos genes
    # no Orion e 87% no arc_vcc_25, entao um gene testado em UMA perturbacao
    # herdava aquele log2 como se fosse a media do contexto. O efeito media por
    # RMS: 1,33 no orion.hct116 contra mediana de 0,033 — a cauda batia no
    # limite de +-5 e dominava tudo. O Replogle, que testa 100% dos genes, nao
    # tem o problema (RMS 0,18, maximo 0,94), o que confirma a causa.
    #
    # Gene com pouca cobertura encolhe para zero, que e "sem estimativa", nao
    # para um palpite.
    bruta = np.where(contagem > 0, soma / np.maximum(contagem, 1), 0.0)
    # Duas defesas, e elas atacam coisas diferentes.
    #
    # A FRACAO corta o vies de selecao: gene testado em poucas perturbacoes tem
    # media positiva por construcao, nao por biologia (ver FRACAO_MINIMA_COBERTURA).
    # Abaixo do corte a media nao e encolhida, e ZERADA — "sem estimativa" e a
    # resposta honesta, nao "estimativa fraca".
    #
    # A CONTAGEM corta ruido amostral dentro dos genes que passaram no corte.
    fracao = contagem / max(linhas.size, 1)
    coberto = fracao >= FRACAO_MINIMA_COBERTURA
    peso = (contagem / (contagem + K_MEDIA)) * coberto
    encolhida = bruta * peso
    print(
        f"  media de {repo}.{contexto} sobre {linhas.size} perturbacoes | "
        f"cobertura mediana {np.median(contagem):.0f} entradas/gene | "
        f"RMS {np.sqrt(np.square(bruta).mean()):.4f} -> {np.sqrt(np.square(encolhida).mean()):.4f}"
    )
    return encolhida


def construir(
    repo: str,
    contexto: str,
    perts_alvo: list[str],
    genes_alvo: list[str],
    media_externa: np.ndarray | None = None,
) -> Resposta:
    """Monta a resposta centrada de um contexto nas perts e genes pedidos.

    `media_externa` substitui a media calculada sobre `perts_alvo`. Use quando
    `perts_alvo` e um subconjunto pequeno — ver media_do_contexto.
    """
    repo_dir = REPOS[repo]
    meta = json.loads((caminho_meta(repo_dir)).read_text(encoding="utf-8"))
    id_contexto = int(meta["context_to_id"][contexto])
    pert_para_id = {str(p): int(i) for p, i in meta["pert_to_id"].items()}
    eixo_repo = {str(g): i for i, g in enumerate(meta["genes"])}

    ctx_ids = np.load(ARRAYS / repo_dir / "preprocessed" / "ctx_ids.npy", allow_pickle=False)
    pert_ids = np.load(ARRAYS / repo_dir / "preprocessed" / "pert_ids.npy", allow_pickle=False)
    # Linha de cada (contexto, pert): o par e unico, verificado no Orion.
    linha_de = {
        (int(c), int(p)): i for i, (c, p) in enumerate(zip(ctx_ids, pert_ids, strict=True))
    }
    linhas = np.array(
        [linha_de[(id_contexto, pert_para_id[p])] for p in perts_alvo], dtype=np.int64
    )
    colunas = np.array([eixo_repo[g] for g in genes_alvo], dtype=np.int64)

    fold_changes = _carregar(repo_dir, "fold_changes", linhas, colunas).astype(np.float64)
    tested = _carregar(repo_dir, "tested", linhas, colunas)
    fdr = _carregar(repo_dir, "fdr", linhas, colunas).astype(np.float32)

    # `fold_changes` JA e a razao linear: mediana de 0,9974 a 1,0006 entre os
    # quatro repos nas entradas testadas. Entao a variavel transferida e o log2
    # dela, direto — sem passar por `ctrl_means`, o que evita duas coisas: a
    # divergencia de estimador (ctrl_means e media de log1p por celula, nosso
    # b_ctrl e log1p do pseudobulk de soma) e a dependencia de um array que o
    # arc_vcc_25 nao publica.
    #
    # `fold_changes` e zero exatamente onde `tested` e False — verificado por
    # igualdade de mascara no Orion. Zero aqui e ausencia de medida, nao razao
    # zero, e tratar os dois como a mesma coisa enviesaria o centramento.
    with np.errstate(divide="ignore", invalid="ignore"):
        log2r = np.log2(fold_changes)
    valido = tested & (fold_changes > 0) & np.isfinite(log2r)
    valido &= np.abs(log2r) <= LIMITE_LOG2R
    log2r = np.where(valido, log2r, 0.0)

    # Centrar por gene usando so as entradas validas: a media removida e a
    # "resposta media do contexto", que e o zero da pontuacao.
    if media_externa is None:
        n_validos = valido.sum(axis=0)
        soma = (log2r * valido).sum(axis=0)
        bruta = np.where(n_validos > 0, soma / np.maximum(n_validos, 1), 0.0)
        # Mesmas duas defesas de media_do_contexto, e pela mesma razao. Se a
        # medicao usar uma regra e a aplicacao outra, o coeficiente medido nao
        # vale para o que foi aplicado — que e exatamente o achado F2 da
        # revisao: 20 das 24 direcoes validadas nao tinham NENHUM gene de
        # cobertura baixa, enquanto 87% da massa aplicada vinha deles.
        fracao = n_validos / max(len(perts_alvo), 1)
        coberto = fracao >= FRACAO_MINIMA_COBERTURA
        media_por_gene = bruta * (n_validos / (n_validos + K_MEDIA)) * coberto
    else:
        media_por_gene = media_externa
    centrado = np.where(valido, log2r - media_por_gene[None, :], 0.0)

    significativos = ((fdr < 0.05) & valido).sum(axis=1)
    print(
        f"  {repo}.{contexto}: {len(perts_alvo)} perts x {len(genes_alvo)} genes | "
        f"entradas validas {valido.mean():.3f} | "
        f"sig/pert mediana {np.median(significativos):.0f} | "
        f"|log2r| centrado mediano {np.median(np.abs(centrado[valido])):.4f}"
    )
    return Resposta(
        repo=repo,
        contexto=contexto,
        perts=perts_alvo,
        genes=genes_alvo,
        log2r=centrado,
        valido=valido,
        media_por_gene=media_por_gene,
        fdr=fdr,
    )


def coeficiente_por_gene(doador: Resposta, alvo: Resposta) -> tuple[np.ndarray, np.ndarray, float]:
    """Coeficiente de transferencia por gene, encolhido para o global.

    Minimos quadrados sem intercepto sobre o centrado — o intercepto ja foi
    removido pelo proprio centramento. Encolhimento tipo James-Stein: gene com
    poucas perturbacoes medidas herda o coeficiente global em vez de inventar
    um proprio.
    """
    mascara = doador.valido & alvo.valido
    n = mascara.sum(axis=0)
    x = doador.log2r * mascara
    y = alvo.log2r * mascara
    numerador = (x * y).sum(axis=0)
    denominador = np.square(x).sum(axis=0)
    beta_global = float(numerador.sum() / denominador.sum()) if denominador.sum() > 0 else 0.0
    beta_cru = np.where(denominador > 0, numerador / np.maximum(denominador, 1e-12), beta_global)
    suficiente = n >= MIN_PERTS_POR_GENE
    peso = np.where(suficiente, n.astype(np.float64), 0.0)
    beta = (peso * beta_cru + K_ENCOLHIMENTO * beta_global) / (peso + K_ENCOLHIMENTO)
    return beta, n, beta_global


def programas_de_resposta(resposta: Resposta, k: int) -> np.ndarray:
    """Os k programas de resposta dominantes do doador, por SVD randomizada.

    Respostas de Perturb-seq nao ocupam o espaco de 18.533 dimensoes: vivem num
    subespaco de poucas dimensoes. O modelo diagonal (um beta por gene) nao
    captura nada disso — trata cada gene como independente.

    SVD randomizada e nao exata porque a exata nao cabe: economia de SVD de uma
    matriz 3.277 x 14.690 custa O(n_perts^2 * n_genes), da ordem de 1e11
    operacoes. A randomizada custa O(n_perts * n_genes * k), da ordem de 1e9,
    com uma iteracao de potencia para nao perder os componentes de cauda.

    A base vem SO do doador. Usar a do alvo seria vazamento — e, para D/E/F,
    impossivel: nao ha resposta medida la de onde tirar base.
    """
    y = resposta.log2r
    rng = np.random.default_rng(0)
    sonda = rng.standard_normal((y.shape[1], k + 10))
    amostra = y @ sonda
    # Uma iteracao de potencia: reaproxima o subespaco dos componentes fortes.
    amostra = y @ (y.T @ amostra)
    q, _ = np.linalg.qr(amostra)
    b = q.T @ y
    _, _, vt = np.linalg.svd(b, full_matrices=False)
    return np.ascontiguousarray(vt[:k].T)


def coeficiente_por_programa(
    doador: Resposta, alvo: Resposta, base: np.ndarray
) -> tuple[np.ndarray, np.ndarray]:
    """Coeficiente de transferencia por programa, nao por gene.

    Projeta as duas respostas na base do doador e ajusta um coeficiente por
    componente. Sao k parametros em vez de 18.533, entao sobra muito menos
    espaco para ajustar ruido — que foi exatamente o que o beta por gene fez no
    Orion, onde o spearman de magnitude DESPENCOU de 0,146 para 0,044 ao passar
    do beta global para o beta por gene.
    """
    carga_doador = doador.log2r @ base
    carga_alvo = alvo.log2r @ base
    numerador = (carga_doador * carga_alvo).sum(axis=0)
    denominador = np.square(carga_doador).sum(axis=0)
    gama = np.where(denominador > 0, numerador / np.maximum(denominador, 1e-12), 0.0)
    return gama, carga_doador


def avaliar(
    previsto: np.ndarray, real: np.ndarray, mascara: np.ndarray, rotulo: str
) -> dict[str, float]:
    """Metricas por perturbacao, na forma que o vcc2026 pontua.

    `cosseno` espelha pds_cosine sobre o delta centrado. `spearman_magnitude` e
    a comparacao que o memorando usa para situar o baseline contra o PIE
    (Rhaister-O 0,741 contra PIE 0,734). `mse_relativo` abaixo de 1 significa
    bater a predicao trivial de resposta media, que e o zero da escala.
    """
    cossenos: list[float] = []
    spearmans: list[float] = []
    sse_modelo = 0.0
    sse_zero = 0.0
    for i in range(real.shape[0]):
        m = mascara[i]
        if m.sum() < 50:
            continue
        p, r = previsto[i][m], real[i][m]
        np_, nr = float(np.linalg.norm(p)), float(np.linalg.norm(r))
        if np_ > 0 and nr > 0:
            cossenos.append(float(p @ r) / (np_ * nr))
        if np.ptp(p) > 0 and np.ptp(r) > 0:
            # Spearman como Pearson dos postos: `spearmanr` devolve um tipo que
            # o verificador nao estreita, e rankdata ja resolve empate por media.
            posto_p = np.asarray(stats.rankdata(np.abs(p)), dtype=np.float64)
            posto_r = np.asarray(stats.rankdata(np.abs(r)), dtype=np.float64)
            spearmans.append(float(np.corrcoef(posto_p, posto_r)[0, 1]))
        sse_modelo += float(np.square(p - r).sum())
        sse_zero += float(np.square(r).sum())
    resultado = {
        "perts_avaliadas": float(len(cossenos)),
        "cosseno_medio": float(np.mean(cossenos)) if cossenos else float("nan"),
        "cosseno_mediano": float(np.median(cossenos)) if cossenos else float("nan"),
        "spearman_magnitude": float(np.mean(spearmans)) if spearmans else float("nan"),
        "mse_relativo": sse_modelo / sse_zero if sse_zero > 0 else float("nan"),
    }
    print(
        f"    {rotulo:28s} cos={resultado['cosseno_medio']:+.4f} "
        f"spearman|mag|={resultado['spearman_magnitude']:+.4f} "
        f"mse/mse_zero={resultado['mse_relativo']:.4f}"
    )
    return resultado


def perts_por_contexto(repo: str) -> dict[str, set[str]]:
    """Perturbacoes medidas em cada contexto de um repo, por nome."""
    repo_dir = REPOS[repo]
    meta = json.loads((caminho_meta(repo_dir)).read_text(encoding="utf-8"))
    id_para_pert = {int(i): str(p) for p, i in meta["pert_to_id"].items()}
    id_para_ctx = {int(i): str(c) for c, i in meta["context_to_id"].items()}
    ctx_ids = np.load(ARRAYS / repo_dir / "preprocessed" / "ctx_ids.npy", allow_pickle=False)
    pert_ids = np.load(ARRAYS / repo_dir / "preprocessed" / "pert_ids.npy", allow_pickle=False)
    saida: dict[str, set[str]] = {}
    for c, p in zip(ctx_ids.tolist(), pert_ids.tolist(), strict=True):
        saida.setdefault(id_para_ctx[int(c)], set()).add(id_para_pert[int(p)])
    return saida


def genes_do_repo(repo: str) -> list[str]:
    meta = json.loads((caminho_meta(REPOS[repo])).read_text(encoding="utf-8"))
    return [str(g) for g in meta["genes"]]


def pares_disponiveis(minimo_perts: int) -> list[tuple[str, str, str, str, list[str], list[str]]]:
    """Todo par (repo_a, ctx_a, repo_b, ctx_b) com perturbacoes e genes em comum.

    Inclui pares ENTRE repos, nao so dentro. Isso nao estava no plano e o dado
    exigiu: a profundidade por perturbacao varia em uma ordem de grandeza entre
    repos (mediana de 2 genes significativos no Orion contra 2.016 no
    arc_vcc_25), entao o par que mede transferencia util pode nao estar dentro
    de um repo so.
    """
    genes_painel = set(pd.read_csv(VCC / "gene_names.csv")["gene_name"].astype(str))
    catalogo = {repo: perts_por_contexto(repo) for repo in REPOS}
    eixos = {repo: [g for g in genes_do_repo(repo) if g in genes_painel] for repo in REPOS}

    entradas: list[tuple[str, str]] = [
        (repo, ctx) for repo, por_ctx in catalogo.items() for ctx in sorted(por_ctx)
    ]
    pares: list[tuple[str, str, str, str, list[str], list[str]]] = []
    for i, (repo_a, ctx_a) in enumerate(entradas):
        for repo_b, ctx_b in entradas[i + 1 :]:
            perts = sorted(catalogo[repo_a][ctx_a] & catalogo[repo_b][ctx_b])
            if len(perts) < minimo_perts:
                continue
            genes = sorted(set(eixos[repo_a]) & set(eixos[repo_b]))
            if len(genes) < 1000:
                continue
            pares.append((repo_a, ctx_a, repo_b, ctx_b, perts, genes))
    return pares


def validar() -> None:
    """Etapa 4: estima o coeficiente num par e aplica a um par que ficou fora."""
    SAIDA.mkdir(parents=True, exist_ok=True)
    pares = pares_disponiveis(MINIMO_PERTS_PAR)
    print(f"pares com >= {MINIMO_PERTS_PAR} perturbacoes em comum: {len(pares)}")

    ajustes: dict[str, dict[str, object]] = {}
    for repo_a, ctx_a, repo_b, ctx_b, perts, genes in pares:
        entre_repos = repo_a != repo_b
        marca = "ENTRE REPOS" if entre_repos else "dentro do repo"
        print(
            f"\n== {repo_a}.{ctx_a} <-> {repo_b}.{ctx_b} | {len(perts)} perts, "
            f"{len(genes)} genes | {marca} =="
        )
        a = construir(repo_a, ctx_a, perts, genes)
        b = construir(repo_b, ctx_b, perts, genes)
        mascara = a.valido & b.valido

        # As duas direcoes, porque o coeficiente e assimetrico: hct116->hepg2
        # da 0,218 e hepg2->hct116 da 0,088 no mesmo par. Quem e doador e quem
        # e alvo muda a resposta, e o dado ja esta carregado.
        for doador, alvo in ((a, b), (b, a)):
            beta, n_perts_gene, beta_global = coeficiente_por_gene(doador, alvo)
            rotulo = f"{doador.repo}.{doador.contexto}->{alvo.repo}.{alvo.contexto}"
            proprios = int((n_perts_gene >= MIN_PERTS_POR_GENE).sum())
            print(
                f"  -- {rotulo} | beta global {beta_global:+.4f} | "
                f"genes com coef. proprio {proprios} --"
            )
            resultados = {
                "zero_resposta_media": avaliar(
                    np.zeros_like(alvo.log2r), alvo.log2r, mascara, "zero (resposta media)"
                ),
                "copia_doador": avaliar(
                    doador.log2r, alvo.log2r, mascara, "copia do doador (beta=1)"
                ),
                "beta_global": avaliar(
                    beta_global * doador.log2r, alvo.log2r, mascara, "beta global"
                ),
                "beta_por_gene": avaliar(
                    beta[None, :] * doador.log2r, alvo.log2r, mascara, "beta por gene encolhido"
                ),
            }
            base = programas_de_resposta(doador, K_PROGRAMAS)
            gama, carga = coeficiente_por_programa(doador, alvo, base)
            previsao_programas = (carga * gama[None, :]) @ base.T
            resultados["gama_por_programa"] = avaliar(
                previsao_programas, alvo.log2r, mascara, f"gama por programa (k={K_PROGRAMAS})"
            )
            # Diagonal e baixa dimensao capturam coisas diferentes: a primeira
            # acerta gene a gene, a segunda acerta a estrutura. A soma das duas
            # entra na escada para dizer se sao complementares ou redundantes.
            previsao_diagonal = beta[None, :] * doador.log2r
            resultados["soma_diagonal_e_programas"] = avaliar(
                0.5 * previsao_diagonal + 0.5 * previsao_programas,
                alvo.log2r,
                mascara,
                "media de diagonal e programas",
            )
            # O peso otimo tem forma fechada: minimizar ||w*d + (1-w)*p - y||^2
            # em w e um problema de uma variavel. Reportar a distribuicao do
            # otimo entre as direcoes diz se 0,5 e uma escolha defensavel ou um
            # palpite redondo que deu sorte.
            diferenca = (previsao_diagonal - previsao_programas) * mascara
            residuo = (previsao_programas - alvo.log2r) * mascara
            denominador_peso = float(np.square(diferenca).sum())
            peso_otimo = (
                -float((residuo * diferenca).sum()) / denominador_peso
                if denominador_peso > 0
                else 0.5
            )
            peso_limitado = float(np.clip(peso_otimo, 0.0, 1.0))
            resultados["peso_otimo"] = avaliar(
                peso_limitado * previsao_diagonal + (1.0 - peso_limitado) * previsao_programas,
                alvo.log2r,
                mascara,
                f"peso otimo w={peso_limitado:.3f}",
            )
            ajustes[rotulo] = {
                "perts": len(perts),
                "genes": len(genes),
                "entre_repos": entre_repos,
                "beta_global": beta_global,
                "beta_por_gene_mediano": float(np.median(beta)),
                "beta_por_gene_p05": float(np.percentile(beta, 5)),
                "beta_por_gene_p95": float(np.percentile(beta, 95)),
                "genes_com_coeficiente_proprio": proprios,
                "peso_otimo_diagonal": peso_limitado,
                "peso_otimo_sem_limite": peso_otimo,
                "metricas": resultados,
            }
            np.savez_compressed(
                SAIDA
                / (
                    f"beta_{doador.repo}_{doador.contexto}"
                    f"_para_{alvo.repo}_{alvo.contexto}.npz"
                ),
                genes=np.array(genes, dtype=object),
                beta=beta,
                n_perts=n_perts_gene,
                beta_global=np.float64(beta_global),
                media_por_gene_doador=doador.media_por_gene,
                media_por_gene_alvo=alvo.media_por_gene,
                # Base e gama viajam juntos: o gama so faz sentido na base em
                # que foi medido, e a base vale so nos genes deste par.
                base_programas=base,
                gama=gama,
                allow_pickle=True,
            )

    # Generalizacao entre pares: o beta global de um par serve para outro?
    globais = {k: float(v["beta_global"]) for k, v in ajustes.items()}  # pyright: ignore[reportArgumentType]
    print("\n== beta global por par ==")
    for chave, valor in sorted(globais.items()):
        print(f"  {chave:46s} {valor:+.4f}")
    if len(globais) >= 2:
        valores = np.array(list(globais.values()))
        print(
            f"  dispersao entre pares: mediana {np.median(valores):.4f}, "
            f"min {valores.min():.4f}, max {valores.max():.4f}, "
            f"coef. de variacao {valores.std() / abs(valores.mean()):.3f}"
        )

    relatorio = {
        "limite_log2r": LIMITE_LOG2R,
        "variavel_transferida": "log2(fold_changes), mascarado por tested",
        "min_perts_por_gene": MIN_PERTS_POR_GENE,
        "k_encolhimento": K_ENCOLHIMENTO,
        "pares": ajustes,
        "beta_global_por_par": globais,
        "leitura": (
            "mse_relativo < 1 significa bater a predicao de resposta media, que e o "
            "zero da escala do vcc2026. cosseno > 0 significa direcao correta. "
            "A dispersao do beta global entre pares e o limite superior da confianca "
            "que se pode ter ao aplicar qualquer beta a D/E/F, onde nao ha medida."
        ),
    }
    (SAIDA / "validacao_loco.json").write_text(
        json.dumps(relatorio, indent=2, ensure_ascii=False), encoding="utf-8"
    )
    print(f"\nescrito em {SAIDA / 'validacao_loco.json'}")


def politica_de_beta(validacao: dict[str, object]) -> tuple[dict[str, str], dict[str, float]]:
    """Qual direcao medida serve de proxy para cada tipo de doador.

    O coeficiente doador -> D/E/F e inestimavel: nao ha perturbacao medida em
    D, E ou F. Entao escolhemos, para cada doador, a direcao MEDIDA que mais se
    parece com a nossa situacao, e dizemos qual foi.

    O que decide a semelhanca e a profundidade do ALVO, nao a do doador. Com o
    doador fixo em orion.hct116, mse/mse0 vai de 0,954 quando o alvo e Replogle
    (profundo) a 0,9984 quando o alvo e o proprio Orion (raso). Com o alvo fixo
    em orion.hct116, o melhor doador nao passa de 0,962. O alvo e o teto.

    Nosso alvo e o ensaio da propria competicao: mediana de ~2.016 genes
    significativos por perturbacao no arc_vcc_25, contra 46 no Replogle e 2 no
    Orion. Portanto a direcao analoga e "doador raso -> alvo profundo", e a
    melhor medida dela e orion.hct116 -> replogle.rpe1.
    """
    pares = validacao["pares"]
    assert isinstance(pares, dict)
    # Direcao de referencia por repo de doador: a de menor mse/mse0 entre as
    # que tem alvo profundo (Replogle ou arc_vcc_25) e amostra suficiente.
    melhor_por_doador: dict[str, tuple[float, str]] = {}
    for rotulo, dados in pares.items():
        assert isinstance(dados, dict)
        doador, alvo = str(rotulo).split("->")
        repo_doador = doador.split(".")[0]
        repo_alvo = alvo.split(".")[0]
        if repo_alvo not in ("replogle", "arc_vcc_25"):
            continue
        if int(dados["perts"]) < 400:
            continue
        mse = float(dados["metricas"]["beta_por_gene"]["mse_relativo"])
        atual = melhor_por_doador.get(repo_doador)
        if atual is None or mse < atual[0]:
            melhor_por_doador[repo_doador] = (mse, str(rotulo))

    referencia = {repo: rotulo for repo, (_, rotulo) in melhor_por_doador.items()}
    mse_esperado = {repo: mse for repo, (mse, _) in melhor_por_doador.items()}
    return referencia, mse_esperado


def ajustar() -> None:
    """Etapa 3: aplica a transferencia medida aos 300 alvos em D, E e F."""
    SAIDA.mkdir(parents=True, exist_ok=True)
    validacao_path = SAIDA / "validacao_loco.json"
    if not validacao_path.exists():
        raise FileNotFoundError(
            "rode `validar` antes de `ajustar`: sem a validacao nao se sabe se o "
            "coeficiente tem direito de ser aplicado a D/E/F"
        )
    validacao = json.loads(validacao_path.read_text(encoding="utf-8"))
    referencia, mse_esperado = politica_de_beta(validacao)
    print("== politica de beta, por repo de doador ==")
    for repo, rotulo in sorted(referencia.items()):
        print(f"  {repo:12s} usa o beta medido em {rotulo} (mse/mse0 {mse_esperado[repo]:.4f})")
    if not referencia:
        raise ValueError("nenhuma direcao com alvo profundo e amostra suficiente")

    genes_painel = pd.read_csv(VCC / "gene_names.csv")["gene_name"].astype(str).tolist()
    alvos = pd.read_csv(VCC / "pert_counts.csv")["target_gene"].astype(str).tolist()
    indice_gene = {g: i for i, g in enumerate(genes_painel)}
    indice_alvo = {a: i for i, a in enumerate(alvos)}
    n_genes, n_alvos = len(genes_painel), len(alvos)

    por_gene = pd.read_csv(CTRL / "por_gene.tsv", sep="\t")
    gate = pd.read_csv(EIXO / "gate_alvos.tsv", sep="\t")
    eixo_resumo = json.loads((EIXO / "resumo.json").read_text(encoding="utf-8"))
    sem_evidencia = set(eixo_resumo["evidencia_doadora"]["alvos_sem_evidencia"])
    reservados = set(eixo_resumo["evidencia_doadora"]["reservados_para_validacao"])

    # Beta por gene de cada direcao de referencia, remapeado para o nosso eixo.
    beta_por_doador: dict[str, np.ndarray] = {}
    # (genes da base, base de programas, gama) por repo de doador. So os repos
    # com direcao de referencia tem base: um gama so faz sentido na base em que
    # foi medido, e nao existe "base mediana" que signifique alguma coisa.
    programas_por_doador: dict[str, tuple[list[str], np.ndarray, np.ndarray]] = {}
    for repo, rotulo in referencia.items():
        doador, alvo = rotulo.split("->")
        caminho = SAIDA / (
            f"beta_{doador.split('.')[0]}_{doador.split('.')[1]}"
            f"_para_{alvo.split('.')[0]}_{alvo.split('.')[1]}.npz"
        )
        pacote = np.load(caminho, allow_pickle=True)
        genes_beta = [str(g) for g in pacote["genes"]]
        valores = np.asarray(pacote["beta"], dtype=np.float64)
        global_beta = float(pacote["beta_global"])
        programas_por_doador[repo] = (
            genes_beta,
            np.asarray(pacote["base_programas"], dtype=np.float64),
            np.asarray(pacote["gama"], dtype=np.float64),
        )
        # Gene fora da direcao de referencia herda o beta global dela.
        espalhado = np.full(n_genes, global_beta, dtype=np.float64)
        for gene, valor in zip(genes_beta, valores, strict=True):
            posicao = indice_gene.get(gene)
            if posicao is not None:
                espalhado[posicao] = valor
        beta_por_doador[repo] = espalhado
        print(f"  beta de {repo}: mediana {np.median(espalhado):+.4f} no nosso eixo")

    # arc_vcc_25 e jiang nao tem nenhuma direcao de referencia: sao doadores sem
    # par de alvo profundo e amostra suficiente. O fallback e a MEDIANA, gene a
    # gene, dos betas que existem — nao a primeira entrada do dicionario, que
    # dependeria da ordem de iteracao e seria uma escolha invisivel.
    #
    # COM DOIS DOADORES MEDIDOS A MEDIANA E A MEDIA DELES, e a palavra "mediana"
    # promete uma robustez a outlier que n=2 nao tem. A aritmetica continua
    # correta — a media de dois betas e a resposta defensavel aqui — mas o numero
    # de betas que entrou fica registrado, porque e ele que diz quanto peso a
    # palavra carrega.
    betas_medidos = np.stack(list(beta_por_doador.values()))
    beta_mediano = np.median(betas_medidos, axis=0)
    sem_referencia = [r for r in REPOS if r not in beta_por_doador]
    if sem_referencia:
        rotulo = "media" if betas_medidos.shape[0] == 2 else "mediana"
        print(
            f"  sem direcao de referencia ({', '.join(sem_referencia)}): usam a "
            f"{rotulo} gene a gene de {betas_medidos.shape[0]} betas medidos "
            f"({', '.join(sorted(beta_por_doador))}), mediana global "
            f"{np.median(beta_mediano):+.4f}"
        )

    # Doadores disponiveis, com a profundidade medida de cada repo. A
    # profundidade pondera a combinacao quando um alvo foi medido em mais de um.
    profundidade = {"arc_vcc_25": 2016.0, "replogle": 46.0, "jiang": 13.0, "orion": 2.0}
    split = json.loads(
        (METAS / "PIE_splits" / "replogle_xdataset" / "train.json").read_text(encoding="utf-8")
    )
    perts_treino: dict[str, set[str]] = {}
    for chave, lista in split.items():
        perts_treino.setdefault(str(chave).split(".", 1)[0], set()).update(
            str(p) for p in lista
        )

    # A RESPOSTA MEDIA do contexto, que e o zero da escala. Medida em
    # runs/vcc/transfer/resposta_media.json, de onde os numeros sao LIDOS e nao
    # copiados: duas revisoes acharam este comentario com valores de uma corrida
    # anterior, e a segunda trouxe a propria correcao tambem obsoleta.
    #
    # Sem esta parcela, `delta = 0` significa "prediga o controle", que pontua
    # ABAIXO de zero — e nao "prediga a resposta media", que pontua zero. Isso
    # atingiria os 143 alvos sem evidencia por inteiro e deslocaria os 157.
    caminho_media = SAIDA / "resposta_media.json"
    if not caminho_media.exists():
        raise FileNotFoundError(
            "rode `resposta_media.py` antes de `ajustar`: sem ele a predicao cola no "
            "controle, que pontua abaixo de zero"
        )
    medida = json.loads(caminho_media.read_text(encoding="utf-8"))
    usar_media = bool(medida["usar_resposta_media_transferida"])
    coeficiente_media = float(medida["coeficiente_mediano"]) if usar_media else 0.0
    print(
        f"resposta media: {'USAR' if usar_media else 'NAO USAR'} | coeficiente "
        f"{coeficiente_media:.4f} | ganho mediano medido "
        f"{float(medida['ganho_mediano_sobre_o_controle']):+.4f} em "
        f"{int(medida['direcoes_com_ganho_positivo'])}/{int(medida['direcoes_totais'])} direcoes"
    )

    # Resposta do doador por alvo, no nosso eixo, centrada e ponderada.
    numerador = np.zeros((n_alvos, n_genes), dtype=np.float64)
    denominador = np.zeros((n_alvos, n_genes), dtype=np.float64)
    # Acumuladores da resposta media agrupada entre contextos doadores.
    media_numerador = np.zeros(n_genes, dtype=np.float64)
    media_denominador = np.zeros(n_genes, dtype=np.float64)
    procedencia: dict[str, list[str]] = {}

    # ========================================================================
    # AGRUPAMENTO DA RESPOSTA MEDIA — laco PROPRIO, e isso e a correcao.
    #
    # Antes a media era acumulada dentro do laco do desvio, atras do gate
    # `set(alvos) & disponiveis & split; if not usaveis: continue`. O gate e
    # correto para o DESVIO, que precisa do nosso alvo medido naquele contexto.
    # A media nao precisa de alvo nosso nenhum: `media_do_contexto()` percorre
    # TODAS as perturbacoes do doador. Um gate servindo duas parcelas com
    # requisitos diferentes excluia os 4 contextos do Replogle — os de melhor
    # transferencia medida de media do conjunto inteiro — porque zero dos
    # nossos 300 alvos esta no split de treino do Replogle.
    #
    # O PESO CONTINUA SENDO PROFUNDIDADE, e isso e uma decisao de NAO mexer.
    #
    # Peso igual venceu por larga margem no primeiro teste (+0,2176 contra
    # +0,1660 de profundidade). Mas aquele teste deixava um CONTEXTO de fora e
    # mantinha os irmaos do mesmo estudo no agrupamento — e com `ARC_H1`
    # retido, profundidade punha 91% do peso em `ARC_H1_VAL`, que e o mesmo
    # ensaio. Atalho de mesmo-estudo, indisponivel para D/E/F, que nao tem
    # doador do proprio estudo.
    #
    # Refeito com `--sem-mesmo-repo` (pool_media_sem_mesmo_repo.json), o sinal
    # fica MISTO: profundidade ganha a mediana (+0,0983 contra +0,0579) e peso
    # igual ganha a media (+0,0873 contra +0,0859), e as medias diferem em
    # 0,0014. Sinal misto da ordem do ruido nao escolhe vencedor, entao o peso
    # fica como estava e o que muda e so o gate — que ganha em 4 de 4
    # configuracoes testadas (cortes de 10, 40 e 200 perturbacoes, com e sem
    # mesmo estudo).
    #
    # Registrado em aberto, porque nao esta resolvido: com o gate fora,
    # profundidade da 95,5% do agrupamento aos dois contextos do arc_vcc_25
    # (2 x 2016 contra 4 x 46 do Replogle e 2 x 2 do Orion) — concentracao
    # MAIOR que os 93,2% de antes, e justamente no doador cuja media e
    # ortogonal a do Orion. Nas duas linhas de Orion do teste sem mesmo
    # estudo, profundidade mede +0,0031 e +0,0009 contra +0,0329 e +0,0038 de
    # peso igual. Qual se parece com D/E/F e desconhecido, e por isso nenhum
    # dos dois pode ser escolhido com a evidencia que existe.
    # ========================================================================
    medias_de_contexto: dict[tuple[str, str], np.ndarray] = {}
    contextos_na_media: list[str] = []
    for repo in REPOS:
        catalogo = perts_por_contexto(repo)
        genes_repo = [g for g in genes_do_repo(repo) if g in indice_gene]
        colunas_media = np.array([indice_gene[g] for g in genes_repo], dtype=np.int64)
        for contexto, disponiveis in sorted(catalogo.items()):
            if len(disponiveis) < MIN_PERTS_MEDIA:
                continue
            media = media_do_contexto(repo, contexto, genes_repo)
            medias_de_contexto[(repo, contexto)] = media
            contextos_na_media.append(f"{repo}.{contexto}")
            tem_media = media != 0.0
            media_numerador[colunas_media] += profundidade[repo] * media * tem_media
            media_denominador[colunas_media] += profundidade[repo] * tem_media
    print(
        f"  resposta media agrupada por PROFUNDIDADE sobre {len(contextos_na_media)} "
        f"contextos (>= {MIN_PERTS_MEDIA} perts, sem o gate do desvio): "
        f"{', '.join(contextos_na_media)}"
    )

    for repo in REPOS:
        meta = json.loads((caminho_meta(REPOS[repo])).read_text(encoding="utf-8"))
        dataset = str(meta["dataset"])
        catalogo = perts_por_contexto(repo)
        genes_repo = [g for g in genes_do_repo(repo) if g in indice_gene]
        for contexto, disponiveis in sorted(catalogo.items()):
            # So o split de treino: usar linha reservada para validacao e vazamento.
            usaveis = sorted(
                set(alvos) & disponiveis & perts_treino.get(dataset, set())
            )
            if not usaveis:
                continue
            print(f"  {repo}.{contexto}: {len(usaveis)} dos 300 alvos no split de treino")
            # Reusa a media ja calculada no laco de cima quando o contexto entrou
            # la: `orion.hct116` percorre 8.987 perturbacoes e ler duas vezes e
            # desperdicio puro. Contexto abaixo de MIN_PERTS_MEDIA nao esta no
            # cache e ainda precisa da media PROPRIA para centrar o desvio — que
            # e uso diferente de entrar no agrupamento.
            media = medias_de_contexto.get((repo, contexto))
            if media is None:
                media = media_do_contexto(repo, contexto, genes_repo)
            resposta = construir(repo, contexto, usaveis, genes_repo, media_externa=media)
            beta = beta_por_doador.get(repo, beta_mediano)
            peso = profundidade[repo]
            colunas = np.array([indice_gene[g] for g in genes_repo], dtype=np.int64)

            # Parte de baixa dimensao: existe so onde a base foi medida. Os
            # genes da base sao um subconjunto (6.177 na direcao de referencia
            # do Orion), e fora dele vale apenas a parte diagonal.
            pacote_programas = programas_por_doador.get(repo)
            if pacote_programas is None:
                linhas_repo = np.empty(0, dtype=np.int64)
                linhas_base = np.empty(0, dtype=np.int64)
                base_doador = np.empty((0, K_PROGRAMAS), dtype=np.float64)
                gama_doador = np.zeros(K_PROGRAMAS, dtype=np.float64)
            else:
                genes_base, base_doador, gama_doador = pacote_programas
                posicao_repo = {g: i for i, g in enumerate(genes_repo)}
                pares_alinhados = [
                    (posicao_repo[g], j) for j, g in enumerate(genes_base) if g in posicao_repo
                ]
                linhas_repo = np.array([p for p, _ in pares_alinhados], dtype=np.int64)
                linhas_base = np.array([j for _, j in pares_alinhados], dtype=np.int64)
                print(
                    f"    baixa dimensao em {len(pares_alinhados)} dos {len(genes_repo)} "
                    f"genes do eixo deste repo"
                )

            base_alinhada = base_doador[linhas_base] if linhas_base.size else base_doador
            colunas_programas = colunas[linhas_repo] if linhas_repo.size else linhas_repo

            for linha, alvo_nome in enumerate(usaveis):
                destino = indice_alvo[alvo_nome]
                valido = resposta.valido[linha].astype(np.float64)
                vetor = resposta.log2r[linha]
                # Diagonal cheia em todo o eixo do repo. Nos genes onde existe
                # base, metade dela e trocada pela parte de programas logo
                # abaixo; fora da base, a diagonal responde sozinha e inteira.
                diagonal = beta[colunas] * vetor * valido
                numerador[destino, colunas] += peso * diagonal
                denominador[destino, colunas] += peso * valido
                if linhas_repo.size:
                    # Baixa dimensao: projeta, escala por programa, reconstroi.
                    carga = (vetor[linhas_repo] * valido[linhas_repo]) @ base_alinhada
                    reconstruido = base_alinhada @ (carga * gama_doador)
                    numerador[destino, colunas_programas] += (
                        0.5 * peso * (reconstruido * valido[linhas_repo] - diagonal[linhas_repo])
                    )
                procedencia.setdefault(alvo_nome, []).append(f"{repo}.{contexto}")

    log2r_previsto = np.where(denominador > 0, numerador / np.maximum(denominador, 1e-12), 0.0)
    com_sinal = int((np.abs(log2r_previsto).sum(axis=1) > 0).sum())
    print(f"\nalvos com algum sinal de doador: {com_sinal} de {n_alvos}")

    media_agrupada = np.where(
        media_denominador > 0, media_numerador / np.maximum(media_denominador, 1e-12), 0.0
    )
    # O coeficiente 0,41 foi medido doador-unico contra alvo-unico. Agrupar
    # varios doadores reduz o ruido do preditor, entao o coeficiente otimo para
    # a media agrupada seria MAIOR (menos atenuacao). Usar o 0,41 e a escolha
    # conservadora, e deliberada: errar por baixo deixa a predicao entre o
    # controle e a media real, que ainda e melhor que o controle; errar por
    # cima passa da media e pode piorar.
    deslocamento = coeficiente_media * media_agrupada
    print(
        f"resposta media agrupada: {int((media_agrupada != 0).sum())} genes com estimativa, "
        f"|media| tipica {np.linalg.norm(media_agrupada) / np.sqrt(n_genes):.4f}, "
        f"deslocamento aplicado |.| {np.linalg.norm(deslocamento) / np.sqrt(n_genes):.4f}"
    )

    # Conversao para o nosso espaco, por contexto, com renormalizacao
    # composicional: a proporcao perturbada e p * 2^log2r RENORMALIZADA para
    # somar 1. Sem isso, muitos genes se movendo deslocam todo o eixo.
    saida: dict[str, np.ndarray] = {}
    diagnostico: dict[str, dict[str, float]] = {}
    for contexto in ("D", "E", "F"):
        b_ctrl = np.asarray(por_gene[f"b_ctrl_{contexto}"], dtype=np.float64)
        proporcao = np.expm1(b_ctrl)
        gate_ctx = gate.loc[gate["contexto"] == contexto].set_index("alvo")
        delta = np.zeros((n_alvos, n_genes), dtype=np.float32)
        n_gateados = 0
        n_com_desvio = 0
        for alvo_nome, destino in indice_alvo.items():
            # A resposta media vale para TODOS os 300. E o efeito comum do
            # procedimento de perturbacao — entrega de guia, selecao, atividade
            # de Cas —, pelo qual todos os 300 passaram, inclusive os que nao
            # transcrevem o proprio alvo. E e o zero da escala: sem ela,
            # delta = 0 significaria "prediga o controle", que pontua abaixo de
            # zero, nao zero.
            #
            # O desvio especifico da perturbacao e que depende de evidencia e de
            # o alvo estar transcrito.
            sem_doador = alvo_nome in sem_evidencia or alvo_nome in reservados
            transcrito = str(gate_ctx.loc[alvo_nome, "gate"]) == "ativo"
            if not sem_doador and not transcrito:
                # CRISPRi reprime transcricao: sem transcrito nao ha o que
                # derrubar, e o desvio especifico seria invencao. A resposta
                # media continua valendo.
                n_gateados += 1
            usa_desvio = not sem_doador and transcrito
            if usa_desvio:
                log2r_total = deslocamento + log2r_previsto[destino]
                n_com_desvio += 1
            else:
                log2r_total = deslocamento
            perturbada = proporcao * np.exp2(log2r_total)
            soma = perturbada.sum()
            if soma <= 0:
                continue
            perturbada *= proporcao.sum() / soma
            delta[destino] = (np.log1p(perturbada) - b_ctrl).astype(np.float32)
        saida[f"delta_{contexto}"] = delta
        movidos = int((np.abs(delta).sum(axis=1) > 0).sum())
        nao_nulos = np.abs(delta) > 0
        diagnostico[contexto] = {
            "alvos_com_delta": float(movidos),
            "alvos_com_desvio_especifico": float(n_com_desvio),
            "alvos_so_com_resposta_media": float(n_alvos - n_com_desvio),
            "alvos_gateados_por_expressao": float(n_gateados),
            "delta_abs_mediano": float(np.median(np.abs(delta[nao_nulos]))),
            "delta_abs_p99": float(np.percentile(np.abs(delta[nao_nulos]), 99)),
        }
        print(
            f"  {contexto}: {movidos}/{n_alvos} alvos com delta | "
            f"{n_com_desvio} com desvio especifico, {n_alvos - n_com_desvio} so com a "
            f"resposta media | {n_gateados} gateados | "
            f"|delta| mediano {diagnostico[contexto]['delta_abs_mediano']:.5f}"
        )

    np.savez_compressed(
        SAIDA / "delta_previsto.npz",
        genes=np.array(genes_painel, dtype=object),
        alvos=np.array(alvos, dtype=object),
        **saida,
        allow_pickle=True,
    )

    # Faixa do peso otimo entre as 24 direcoes, lida da validacao em vez de
    # transcrita: o literal anterior era de uma corrida que nao existe mais.
    pesos_otimos = np.array(
        [float(par["peso_otimo_diagonal"]) for par in validacao["pares"].values()],
        dtype=np.float64,
    )

    resumo = {
        "formato_saida": {
            "arquivo": str(SAIDA / "delta_previsto.npz"),
            "genes": "array de 18.533 simbolos, na ordem de gene_names.csv",
            "alvos": "array de 300 simbolos, na ordem de pert_counts.csv",
            "delta_D/delta_E/delta_F": "(300, 18533) float32, delta em espaco log1p de CP50k",
            "convencao": (
                "b_previsto = b_ctrl + delta. Delta zero significaria predizer o "
                "CONTROLE, que pontua abaixo de zero — por isso nenhum dos 300 alvos "
                "tem delta zero: todos carregam a resposta media estimada do contexto, "
                "que e o zero da escala."
            ),
        },
        "resposta_media": {
            "usada": usar_media,
            "coeficiente_aplicado": coeficiente_media,
            "ganho_mediano_medido": float(medida["ganho_mediano_sobre_o_controle"]),
            "direcoes_com_ganho_positivo": int(medida["direcoes_com_ganho_positivo"]),
            "direcoes_totais": int(medida["direcoes_totais"]),
            "correlacao_mediana": float(medida["correlacao_mediana"]),
            "por_que": (
                "O zero da escala do vcc2026 e a resposta media real do contexto as 300 "
                "perturbacoes, nao o controle. Sem esta parcela os 143 alvos sem "
                "evidencia pontuariam ABAIXO de zero em vez de zero. A parcela vale para "
                "todos os 300, inclusive os gateados: e o efeito comum do procedimento "
                "de perturbacao, pelo qual todos passaram. Transferiu em 24 de 24 "
                "direcoes medidas."
            ),
            "conservadorismo": (
                "O coeficiente foi medido doador-unico contra alvo-unico. Agrupar "
                "doadores reduz o ruido do preditor, entao o otimo para a media agrupada "
                "seria maior. Usar o valor medido e errar por baixo de proposito: a "
                "predicao fica entre o controle e a media real, que ainda bate o "
                "controle, em vez de passar da media."
            ),
        },
        "politica_de_beta": {
            "direcao_de_referencia_por_doador": referencia,
            "mse_relativo_esperado": mse_esperado,
            "por_que": (
                "O coeficiente doador -> D/E/F e inestimavel: zero perturbacao medida "
                "em D, E ou F. A profundidade do ALVO e o teto — com doador fixo em "
                "orion.hct116, mse/mse0 vai de 0,954 (alvo Replogle) a 0,9984 (alvo "
                "Orion). Nosso alvo e o ensaio da competicao, profundo, entao a "
                "direcao analoga e doador raso -> alvo profundo."
            ),
        },
        "ponderacao_por_profundidade": profundidade,
        "doadores_sem_direcao_de_referencia": sem_referencia,
        "fallback": {
            "regra": "mediana gene a gene dos betas medidos",
            "betas_que_entraram": sorted(beta_por_doador),
            "n_betas": int(betas_medidos.shape[0]),
            "ressalva": (
                "Com dois betas medidos a mediana E a media dos dois: nao ha "
                "robustez a outlier, so ausencia de escolha arbitraria entre eles."
                if betas_medidos.shape[0] == 2
                else "Mediana com mais de dois betas: robusta a um outlier."
            ),
        },
        "alvos_com_sinal_de_doador": com_sinal,
        "alvos_sem_evidencia": len(sem_evidencia),
        "alvos_reservados_para_validacao": sorted(reservados),
        "por_contexto": diagnostico,
        "limite_honesto": (
            "Duas parcelas, de tamanhos muito diferentes. A resposta media do contexto "
            f"transferiu em {int(medida['direcoes_com_ganho_positivo'])} de "
            f"{int(medida['direcoes_totais'])} direcoes, com ganho mediano de "
            f"{float(medida['ganho_mediano_sobre_o_controle']):.3f} do erro "
            "quadratico sobre predizer o controle — e a parcela grande, e vale para os "
            "300 alvos. O desvio especifico da perturbacao tambem bate o zero em 24 de "
            "24, mas por margem pequena: mse/mse0 mediano 0,955 com o modelo combinado, "
            "melhor 0,872, e vale so para os 157 com evidencia. Copiar o doador cru "
            "(beta=1) pontua de 1,2 a 3,4 vezes PIOR que nao predizer nada, entao a "
            "magnitude do encolhimento vale mais que qualquer sofisticacao de modelo."
        ),
        "peso_diagonal_contra_programas": {
            "usado": 0.5,
            "otimo_medido_mediano": round(float(np.median(pesos_otimos)), 4),
            "faixa_medida": [
                round(float(pesos_otimos.min()), 4),
                round(float(pesos_otimos.max()), 4),
            ],
            "por_que_fica_em_meio_a_meio": (
                "O otimo rende 0,0003 de mse/mse0 DENTRO da amostra e a faixa entre "
                f"direcoes vai de {pesos_otimos.min():.2f} a {pesos_otimos.max():.2f} "
                "— nao esta determinado. Perseguir isso custaria uma reemissao inteira "
                "por ganho abaixo do ruido. Mediana e faixa saem de validacao_loco.json "
                "a cada corrida, nao de literal no codigo."
            ),
        },
    }
    (SAIDA / "transfer_fit.json").write_text(
        json.dumps(resumo, indent=2, ensure_ascii=False), encoding="utf-8"
    )
    print(f"\nescrito em {SAIDA / 'delta_previsto.npz'} e transfer_fit.json")


def main() -> None:
    analisador = argparse.ArgumentParser(description=__doc__)
    analisador.add_argument("subcomando", choices=("validar", "ajustar"))
    argumentos = analisador.parse_args()
    if argumentos.subcomando == "validar":
        validar()
    else:
        ajustar()


if __name__ == "__main__":
    main()
