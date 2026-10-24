"""Do delta em `bulk_lognorm` ao alvo de composicao do grupo.

O arquivo de delta promete `m[g] = b[g] + delta[g]` no espaco em que o placar
compara. A composicao que isso pede e `q[g] = expm1(m[g])/ALVO_SOMA`, porque
`bulk_lognorm` e `log1p(ALVO_SOMA * q)`. Tres coisas podem estar erradas no
pedido, e cada uma tem uma resposta diferente:

**1. `m < 0` e impossivel.** Contagem nao negativa implica `bulk_lognorm >= 0`.
Com `w` mediano de 0,0028 fora do gate, qualquer `delta <= -0,003` ali pede o
impossivel. A resposta e `max(m, 0)` -- a projecao ortogonal no fecho do
conjunto factivel. Nao e um botao, e a fronteira. Medido no `lfc_D` resgatado:
22,95% dos pares sao clampeados e 20,9% da variancia do delta se perde ali, o
que derruba a inclinacao contra o delta CRU para 0,769. Esse numero e
reportado, nao suprimido: o conserto e gravar `delta >= -b[g]` no arquivo de
delta, e ai a inclinacao contra o proprio delta sai +0,9997.

**2. Zerar um gene do gate e patologico.** `m = 0` num gene do gate significa
emitir zero contagem nas 400 celulas, e `de_lfc_nmae` le
`log2((0 + 1e-9)/ref)`, cerca de -33, quando nao predizer nada vale 1,0 naquele
gene. Medido em D sem piso: 626 pares do gate com `q == 0`, 544 genes, 251 dos
300 alvos. Com 400 celulas contra 18.400 o Mann-Whitney satura -- -4 e -8 sao
indistinguiveis em significancia -- e `de_lfc_nmae` e LINEAR no erro, logo
sobretiro custa proporcionalmente. Dai um piso declarado, so no gate, so onde o
LFC e lido. `PISO_LOG2 = 4,0` custa 0,008 de inclinacao; 3,0 custa 0,019.

**3. A composicao nao fecha.** `sum_g expm1(m0) != ALVO_SOMA`. Fecha-se com UM
escalar ADITIVO em espaco `bulk_lognorm`, nao multiplicativo na razao: um
escalar sobre a razao `R` nao muda nada, porque `pi*R/sum(pi*R)` e invariante a
escalar. Um `lambda` em `m` muda o `bulk_lognorm` emitido de forma UNIFORME
entre genes, logo cai inteiro no intercepto da regressao de
`validar_numerico.py` e nao toca a inclinacao nem o `deslocamento_por_gene`.

NAO existe parametro de amplitude aqui. `alpha`, `k`, `k_desvio`: nenhum, e o
motivo e numerico. A conversao e exata, logo o delta emitido seria exatamente
`alpha*delta` e a inclinacao medida por `validar_numerico._ajustar` SERIA
`alpha` -- com residuo zero para disfarcar. `FAIXA = (0,85; 1,15)` reprovaria
qualquer `alpha` fora disso. Amplitude e propriedade do arquivo de delta.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from momentos import MomentosControle
from pseudobulk import ALVO_SOMA

PISO_LOG2: float = 4.0
TETO_LOG2: float = 6.0
ITERACOES_DE_BISSECCAO: int = 200
LIMITE_DE_BISSECCAO: float = 40.0


@dataclass(frozen=True)
class AlvoDeComposicao:
    q: np.ndarray  # [n_genes] float64, >= 0, soma 1
    razao: np.ndarray  # [n_genes] float64, q/pi -- EXATAMENTE 1 onde delta_efetivo == 0
    delta_efetivo: np.ndarray  # [n_genes] float64
    lambda_de_fechamento: float
    residuo_de_fechamento: float
    n_no_piso_factivel: int  # m0 + lambda <= 0
    n_no_piso_do_gate: int
    n_no_teto_do_gate: int
    log2fc_min_no_gate: float
    log2fc_max_no_gate: float


def fechar_por_lambda(
    m0: np.ndarray,
    alvo_soma: float = ALVO_SOMA,
    iteracoes: int = ITERACOES_DE_BISSECCAO,
    limite: float = LIMITE_DE_BISSECCAO,
) -> tuple[float, float]:
    """`(lambda, residuo)` de `f(l) = sum_g expm1(max(m0+l, 0)) - alvo_soma`.

    `f` e continua e nao decrescente, vale `-alvo_soma` em `-inf` e cresce sem
    limite, logo a raiz e unica e a bissecao e o metodo certo. Newton nao e: a
    derivada e exatamente zero nos genes clampeados em `max(.,0)`, e com 23% dos
    pares clampeados ele nao converge.
    """

    def f(lam: float) -> float:
        return float(np.expm1(np.maximum(m0 + lam, 0.0)).sum() - alvo_soma)

    f0 = f(0.0)
    if f0 == 0.0:
        return 0.0, 0.0
    baixo, alto = -limite, limite
    for _ in range(iteracoes):
        meio = 0.5 * (baixo + alto)
        if f(meio) < 0.0:
            baixo = meio
        else:
            alto = meio
    lam = 0.5 * (baixo + alto)
    return lam, f(lam)


def alvo_de_composicao(
    momentos: MomentosControle,
    delta_linha: np.ndarray,
    piso_log2: float = PISO_LOG2,
    teto_log2: float = TETO_LOG2,
    alvo_soma: float = ALVO_SOMA,
) -> AlvoDeComposicao:
    """O alvo de composicao de UM grupo: clamp na fronteira factivel e fechamento."""
    b = momentos.bulk
    delta = np.asarray(delta_linha, dtype=np.float64)
    if delta.shape != b.shape:
        raise ValueError(f"delta tem forma {delta.shape}, esperado {b.shape}")
    gate = momentos.gate
    w = momentos.w

    # Piso e teto SO no gate, e em espaco `bulk_lognorm`, nao em log2 da razao:
    # `log1p(w * 2^+-k)` e o nivel que corresponde a um log2FC de +-k.
    lo = np.full_like(b, -np.inf)
    hi = np.full_like(b, np.inf)
    if np.isfinite(piso_log2):
        lo[gate] = np.log1p(w[gate] * 2.0**-piso_log2)
    if np.isfinite(teto_log2):
        hi[gate] = np.log1p(w[gate] * 2.0**teto_log2)
    pedido = b + delta
    m1 = np.clip(pedido, lo, hi)
    n_piso_gate = int((gate & (pedido < lo)).sum())
    n_teto_gate = int((gate & (pedido > hi)).sum())

    m0 = np.maximum(m1, 0.0)
    lam, residuo = fechar_por_lambda(m0, alvo_soma)
    m = np.maximum(m0 + lam, 0.0)
    q = np.expm1(m) / alvo_soma
    delta_efetivo = m - b

    # A razao de composicao contra o pool, que e o que o kernel multiplica pelas
    # contagens da doadora. E `1,0` EXATO onde o delta efetivo e exatamente zero,
    # e isso nao e um snap de tolerancia: um gene cujo nivel pretendido nao muda
    # tem de emitir as contagens da doadora sem mudanca. Calcular
    # `expm1(m)/w` ali daria `1 +- 4,44e-16`, porque `expm1(log1p(w))` nao
    # devolve `w` exato (medido: erro de 2,27e-13), e `floor(D * (1 - 1e-16))` e
    # `D - 1` -- um ulp na razao decrementa toda contagem nao nula do grupo.
    #
    # Sob delta nulo o lambda e 0 e o delta efetivo e 0 em todos os genes, logo a
    # razao e identicamente 1 e a emissao devolve as doadoras bit a bit. Com
    # delta real o lambda e da ordem de 1e-3 e nenhum gene fica em 1,0 -- o que
    # esta certo: se a composicao muda, a fatia de todo gene muda de verdade.
    inalterado = delta_efetivo == 0.0
    razao = np.ones_like(q)
    calcular = ~inalterado & (w > 0)
    razao[calcular] = np.expm1(m[calcular]) / w[calcular]

    # log2FC emitido, so onde o placar le: razao de composicao contra o pool.
    no_gate = gate & (w > 0)
    if np.any(no_gate):
        with np.errstate(divide="ignore"):
            log2fc = np.log2(np.expm1(m[no_gate]) / w[no_gate])
        log2fc = log2fc[np.isfinite(log2fc)]
        minimo = float(log2fc.min()) if log2fc.size else float("nan")
        maximo = float(log2fc.max()) if log2fc.size else float("nan")
    else:
        minimo = maximo = float("nan")

    return AlvoDeComposicao(
        q=q,
        razao=razao,
        delta_efetivo=delta_efetivo,
        lambda_de_fechamento=lam,
        residuo_de_fechamento=residuo,
        # `< 0` e nao `<= 0`: o piso so ATUA onde mudaria o valor. Com `<=`, os
        # 2.975 genes de `b == 0` contariam sob delta nulo e a identidade do T02
        # pareceria quebrada por um contador.
        n_no_piso_factivel=int((m0 + lam < 0.0).sum()),
        n_no_piso_do_gate=n_piso_gate,
        n_no_teto_do_gate=n_teto_gate,
        log2fc_min_no_gate=minimo,
        log2fc_max_no_gate=maximo,
    )


@dataclass(frozen=True)
class Viabilidade:
    contexto: str
    fracao_de_pares_infactiveis: float  # sobre TODOS os pares
    fracao_dos_nao_nulos_infactiveis: float  # sobre os pares de delta != 0
    genes_atingidos: int
    variancia_perdida: float
    inclinacao_contra_delta_cru: float
    massa_em_genes_inalcancaveis: float


def viabilidade_do_delta(momentos: MomentosControle, delta_matriz: np.ndarray) -> Viabilidade:
    """Mede, em segundos e sem emitir uma celula, se o arquivo de delta e factivel.

    A quantidade que importa e a VARIANCIA do delta perdida no clamp, nao a massa
    L1: com 98,9% de zeros exatos, a massa em genes inalcancaveis parece
    desprezivel (1,3e-03 no `lfc_D`) enquanto a variancia ali vale 7%, e e a
    variancia que entra na inclinacao da regressao do validador numerico.

    As tres fracoes sao centradas, como `validar_numerico._ajustar` centra, para
    que `inclinacao_contra_delta_cru` possa ser comparada com o que o validador
    vai medir no arquivo emitido. Medido no `lfc_D`: a previsao analitica da
    0,7692 e o validador mede 0,7623 sobre o arquivo emitido -- 0,007 de
    diferenca, que e o custo do piso do gate e do lambda.

    Duas fracoes, e as duas sao necessarias porque a primeira engana. Com 98,9%
    de zeros exatos no `lfc_D`, a fracao sobre TODOS os pares vale 0,0026 e se
    le como desprezivel; sobre os pares que o preditor realmente prediz ela vale
    0,2393. E a segunda que diz que quase um quarto das predicoes pede contagem
    negativa. Nenhum portao deve usar a primeira.

    `variancia_perdida` e `1 - Var(clipado)/Var(cru)`, e NAO e o complemento da
    inclinacao: clipar nao e projetar na direcao do delta. Medido, as duas
    diferem por 0,022 (0,7473 contra 0,7692). Quem decide se a submissao passa e
    a inclinacao.
    """
    d = np.asarray(delta_matriz, dtype=np.float64)
    b = momentos.bulk[None, :]
    infactivel = b + d < 0.0
    clipado = np.maximum(b + d, 0.0) - b

    x = d.ravel()
    y = clipado.ravel()
    xm, ym = x.mean(), y.mean()
    vx = float(((x - xm) ** 2).mean())
    vy = float(((y - ym) ** 2).mean())
    cxy = float(((x - xm) * (y - ym)).mean())

    inalcancavel = np.broadcast_to(~momentos.alcancavel, d.shape) & (d > 0.0)
    parcela = float((((d - xm) * inalcancavel) ** 2).mean() / vx) if vx > 0 else float("nan")
    nao_nulos = int((d != 0.0).sum())

    return Viabilidade(
        contexto=momentos.contexto,
        fracao_de_pares_infactiveis=float(infactivel.mean()),
        fracao_dos_nao_nulos_infactiveis=(
            float(infactivel.sum() / nao_nulos) if nao_nulos else 0.0
        ),
        genes_atingidos=int(infactivel.any(axis=0).sum()),
        variancia_perdida=float(1.0 - vy / vx) if vx > 0 else float("nan"),
        inclinacao_contra_delta_cru=float(cxy / vx) if vx > 0 else float("nan"),
        massa_em_genes_inalcancaveis=parcela,
    )
