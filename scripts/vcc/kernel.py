"""O kernel de emissao de um grupo (contexto, alvo): 400 celulas de contagem inteira.

A ideia e transportar celulas REAIS de controle em vez de sortear de um modelo.
Um modelo parametrico -- multinomial sobre o perfil medio, Poisson em torno das
contagens, gama-misturado -- produz 400 celulas mais parecidas entre si do que
as reais, e `de_sig_jaccard` e `de_direction_fidelity` leem exatamente essa
diferenca: celulas identicas davam cerca de 9.700 genes "significativos" por
perturbacao contra cerca de 2.700 reais.

O caminho tem seis passos e UMA propriedade que os amarra:

    sob `delta == 0`, a emissao devolve as doadoras bit a bit.

Isso nao e enfeite. Significa que nenhum vies entra pelo caminho, e que o ruido
de amostragem das 400 doadoras fica DENTRO do espalhamento entre perturbacoes
que o `expr_mse` usa como orcamento -- logo `r = min(1, B_pred/claim)` sai 1 e
o aviso de `PRED_CORRECTION_BUDGET_FLAG_RATIO = 0,7` nao dispara. As duas coisas
pelo mesmo mecanismo.

Duas decisoes que parecem detalhe e nao sao:

**Alcancabilidade e por `pi` do POOL, nao por `Q` do sorteio.** `M[g] =
Q[g]*q[g]/pi[g]` devolve `Q[g]` exato sob delta nulo, inclusive quando `Q[g]`
e zero. Testar `Q == 0` em vez de `pi == 0` faz genes presentes no pool mas
ausentes das 400 celulas sorteadas receberem `S0*q` do ar, e a identidade
quebra (medido: `max|M - Q|` passa de 2,9e-11 para mais de 1).

**A inteirizacao e PPS sistematico, nunca `np.round`.** Arredondar desloca a
media do grupo em 0,6% no arranjo medido, e a media do grupo e exatamente o que
`de_lfc_nmae` e `expr_mse` cobram. E nunca `rng.binomial` sobre o `.data` do
h5ad: `Cannot cast array data from dtype('float32') to dtype('int64')`.
"""

from __future__ import annotations

from dataclasses import dataclass
from hashlib import blake2b

import numpy as np
import scipy.sparse as sp
from conversao import PISO_LOG2, TETO_LOG2, AlvoDeComposicao, alvo_de_composicao
from momentos import MomentosControle
from pseudobulk import bulk_lognorm

N_CELULAS: int = 400
TETO_DE_CONTAGEM: int = 1_000_000
TENTATIVAS_DE_LINHA_NULA: int = 8


def _hash63(texto: str) -> int:
    """blake2b de 8 bytes truncado em 63 bits: cabe num int64 sem virar negativo."""
    return int.from_bytes(blake2b(texto.encode("utf-8"), digest_size=8).digest(), "big") % (2**63)


def semente_do_grupo(semente_base: int, contexto: str, alvo: str) -> np.random.SeedSequence:
    """Semente derivada do ROTULO do grupo, nao da ordem de emissao.

    Com `spawn_key` vindo do hash do contexto e do alvo, reemitir um unico grupo
    devolve a mesma matriz sem precisar reproduzir os 899 outros. Um `rng` global
    consumido em sequencia pelos 900 grupos nao tem essa propriedade, e e o que
    torna impossivel consertar um grupo so.
    """
    return np.random.SeedSequence(
        entropy=semente_base, spawn_key=(_hash63(contexto), _hash63(alvo))
    )


def sortear_doadoras(
    rng: np.random.Generator, lote: np.ndarray, n_celulas: int = N_CELULAS
) -> np.ndarray:
    """`n_celulas` indices de controle, COM REPOSICAO, estratificados por `ntc_id`.

    Com reposicao porque e isso que justifica `tr(Sigma)/n` como a correcao de
    jackknife certa, sem fator de populacao finita -- o mesmo regime em que o
    `baseline` do organizador calcula a dele.

    Estratificado porque a biblioteca media por lote varia de 20.313 a 22.462:
    um sorteio simples deixaria a profundidade do grupo flutuar por composicao
    de lote, e profundidade entra direto no CPM que a DE le. Com 46 lotes e 400
    celulas, 32 lotes dao 9 celulas e 14 dao 8; quais dao 9 e sorteado.
    """
    estratos = np.unique(lote)
    n = estratos.size
    base, resto = divmod(n_celulas, n)
    quantos = np.full(n, base, dtype=np.int64)
    quantos[rng.permutation(n)[:resto]] += 1
    partes = []
    for e, k in zip(estratos, quantos, strict=True):
        candidatos = np.flatnonzero(lote == e)
        if candidatos.size == 0:
            raise ValueError(f"estrato {e} vazio")
        partes.append(candidatos[rng.integers(0, candidatos.size, size=int(k))])
    return np.concatenate(partes).astype(np.int64)


def massa_por_gene(
    q: np.ndarray,
    razao: np.ndarray,
    alcancavel: np.ndarray,
    somas_doadoras: np.ndarray,
    soma_total_doadoras: float,
) -> np.ndarray:
    """A soma de contagens que cada gene tem de ter no grupo emitido.

    Gene alcancavel: `Q[g] * razao[g]`. A razao vem da conversao e vale 1,0
    EXATO onde o delta efetivo e zero, logo `M == Q` exato ali -- e dai sai a
    identidade bit a bit. Calcular `Q*q/pi` em vez disso introduz dois
    arredondamentos independentes e a razao sai `1 +- 4,4e-16`, o que basta para
    `floor` decrementar toda contagem nao nula.

    Gene inalcancavel (`pi == 0`): `S0 * q`. E este ramo que recupera os 2.975
    genes que o pool de D nunca viu; sem ele, delta positivo num gene desses
    simplesmente nao sai, e a inclinacao contra o delta cai de 0,998 para 0,893.
    """
    return np.where(alcancavel, somas_doadoras * razao, soma_total_doadoras * q)


def distribuir_por_celula(
    doadoras: np.ndarray,
    bibliotecas: np.ndarray,
    massa: np.ndarray,
    somas_doadoras: np.ndarray,
    alcancavel: np.ndarray,
) -> np.ndarray:
    """Espalha `massa[g]` pelas 400 celulas. Uma regra, dois ramos, continua em `M = Q`.

    Para baixo (`M <= Q`, gene alcancavel) a regra e multiplicativa e o suporte
    fica contido no da doadora -- celula que tinha zero continua com zero, que e
    o comportamento certo de um knockdown. Para cima, o excedente `M - Q` e
    distribuido proporcional a biblioteca da celula, o que alcanca celulas em que
    a doadora tinha zero e faz o `zerofrac` cair na direcao certa sem nenhuma
    curva ajustada. E o mesmo caminho dos genes de `pi == 0`, em que `Q == 0`.

    `sum_c Y[c,g] == massa[g]` nos dois ramos, exatamente (medido: 5,8e-11).

    Nao existe ramo especial para gene alcancavel de doadora nula, e a ausencia e
    deliberada: `massa = Q*razao` ja da zero quando `Q == 0`, e o ramo de subida
    entao devolve `D + (0-0)*L/S0 = 0`. Uma guarda `y[:, alcancavel & (Q==0)] = 0`
    estava aqui e saiu depois do teste de mutacao: ela e inerte sob operacao
    correta e ABSORVIA a mutacao que troca a alcancabilidade do pool pela do
    sorteio, que sem ela a identidade sob delta nulo pega na hora. Codigo morto
    que engole mutacao e pior que codigo nenhum.
    """
    total = float(bibliotecas.sum())
    if total <= 0:
        raise ValueError("doadoras com biblioteca total nula")
    para_baixo = alcancavel & (massa <= somas_doadoras) & (somas_doadoras > 0)

    fator = np.ones_like(massa)
    np.divide(massa, somas_doadoras, out=fator, where=para_baixo)
    y = doadoras * np.where(para_baixo, fator, 1.0)

    # Ramo de subida (e os inalcancaveis): soma, nao escala.
    subida = ~para_baixo
    if np.any(subida):
        excedente = np.where(subida, massa - somas_doadoras, 0.0)
        y[:, subida] = doadoras[:, subida] + np.outer(
            bibliotecas / total, excedente[subida]
        )
    return np.maximum(y, 0.0)


def inteirizar_por_coluna(fracionario: np.ndarray, rng: np.random.Generator) -> np.ndarray:
    """PPS sistematico por coluna: um deslocamento uniforme por gene, sem ordenacao.

    `picks = floor(cumsum(frac) + u) - floor(cumsum_anterior + u)`, que vale 0 ou
    1 por construcao porque `frac < 1`. Propriedades:

    * `E[Z] == Y` exatamente, para `u ~ U(0,1)` -- medido 2,7e-05 de erro
      relativo da media em 500 sementes, contra 0,60% do `np.round`;
    * `|sum_c Z[c,g] - massa[g]| <= 1` sempre, porque a soma da coluna e
      `floor(sum frac + u)`;
    * `frac == 0` implica `picks == 0`, logo delta nulo devolve a doadora.

    Custo O(n_celulas * n_genes), sem ordenar nada -- a alocacao por maior-resto
    custaria uma ordenacao por coluna para comprar 3e-03 de exatidao contra um
    residuo de amostragem de 3e-02 que mantemos de proposito.
    """
    y = np.asarray(fracionario, dtype=np.float64)
    piso = np.floor(y)
    frac = y - piso
    u = rng.random(y.shape[1])
    acumulado = np.cumsum(frac, axis=0) + u[None, :]
    chao = np.floor(acumulado)
    picks = np.empty_like(chao)
    picks[0] = chao[0] - np.floor(u)
    picks[1:] = chao[1:] - chao[:-1]
    return (piso + picks).astype(np.int64)


@dataclass(frozen=True)
class GrupoEmitido:
    contexto: str
    alvo: str
    matriz: sp.csr_matrix
    semente: int
    nnz: int
    soma_maxima_por_celula: int
    soma_minima_por_celula: int
    linhas_nulas: int
    n_genes_criados: int
    zerofrac_no_gate: float
    residuo_bulk_rms: float
    residuo_bulk_rms_no_gate: float
    log2fc_do_alvo: float
    alvo_de_composicao: AlvoDeComposicao


def emitir_grupo(
    momentos: MomentosControle,
    pool: sp.csr_matrix,
    delta_linha: np.ndarray,
    alvo: str,
    semente_base: int,
    piso_log2: float = PISO_LOG2,
    teto_log2: float = TETO_LOG2,
    n_celulas: int = N_CELULAS,
) -> GrupoEmitido:
    """As 400 celulas de um grupo, em CSR canonico de valores inteiros."""
    semente = semente_do_grupo(semente_base, momentos.contexto, alvo)
    rng = np.random.default_rng(semente)
    alvo_q = alvo_de_composicao(momentos, delta_linha, piso_log2, teto_log2)

    indices = sortear_doadoras(rng, momentos.lote, n_celulas)
    for _ in range(TENTATIVAS_DE_LINHA_NULA):
        doadoras = np.asarray(pool[indices].todense(), dtype=np.float64)
        bibliotecas = doadoras.sum(axis=1)
        if np.all(bibliotecas > 0):
            break
        # A doadora de menor nnz do pool tem 787 entradas, entao isto e
        # praticamente impossivel -- mas uma linha de soma zero reprova o item
        # `celula_nao_nula` do contrato, e o preco de re-sortear e nulo.
        nulas = np.flatnonzero(bibliotecas <= 0)
        indices = indices.copy()
        indices[nulas] = rng.integers(0, momentos.n_celulas, size=nulas.size)
    else:
        raise ValueError(f"{momentos.contexto}/{alvo}: linha nula apos as tentativas")

    somas = doadoras.sum(axis=0)
    massa = massa_por_gene(
        alvo_q.q, alvo_q.razao, momentos.alcancavel, somas, float(doadoras.sum())
    )
    y = distribuir_por_celula(doadoras, bibliotecas, massa, somas, momentos.alcancavel)
    z = inteirizar_por_coluna(y, rng)

    por_celula = z.sum(axis=1)
    if int(por_celula.max()) > TETO_DE_CONTAGEM:
        raise ValueError(
            f"{momentos.contexto}/{alvo}: celula com {int(por_celula.max())} contagens, "
            f"acima do teto de {TETO_DE_CONTAGEM}"
        )
    matriz = sp.csr_matrix(z.astype(np.float32))
    matriz.eliminate_zeros()
    matriz.sort_indices()

    soma_emitida = z.sum(axis=0).astype(np.float64)
    emitido_bulk = bulk_lognorm(soma_emitida[None, :])[0]
    pretendido = momentos.bulk + alvo_q.delta_efetivo
    residuo = emitido_bulk - pretendido
    gate = momentos.gate
    j_alvo = momentos.genes.index(alvo) if alvo in momentos.genes else -1
    if j_alvo >= 0 and momentos.w[j_alvo] > 0 and soma_emitida[j_alvo] >= 0:
        log2fc_alvo = float(np.log2(max(np.expm1(emitido_bulk[j_alvo]), 1e-12) / momentos.w[j_alvo]))
    else:
        log2fc_alvo = float("nan")

    return GrupoEmitido(
        contexto=momentos.contexto,
        alvo=alvo,
        matriz=matriz,
        semente=_hash63(f"{semente_base}|{momentos.contexto}|{alvo}"),
        nnz=int(matriz.nnz),
        soma_maxima_por_celula=int(por_celula.max()),
        soma_minima_por_celula=int(por_celula.min()),
        linhas_nulas=int((por_celula == 0).sum()),
        n_genes_criados=int(((somas == 0) & (soma_emitida > 0)).sum()),
        zerofrac_no_gate=float((z[:, gate] == 0).mean()) if np.any(gate) else float("nan"),
        residuo_bulk_rms=float(np.sqrt((residuo**2).mean())),
        residuo_bulk_rms_no_gate=(
            float(np.sqrt((residuo[gate] ** 2).mean())) if np.any(gate) else float("nan")
        ),
        log2fc_do_alvo=log2fc_alvo,
        alvo_de_composicao=alvo_q,
    )
