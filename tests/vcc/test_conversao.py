"""Do delta em `bulk_lognorm` ao alvo de composicao: a conversao e seus clamps.

Esta e a peca em que os tres projetos de emissao divergiam, e cada teste aqui
nomeia a mudanca de producao que o faz falhar. Quatro armadilhas que os
projetos descartados tinham e que T01/T02/T06/T08 pegam:

* `psi = expm1(b+delta)/expm1(b)` da 0/0 nos genes de `b == 0` (2.975 em D) e
  NaN sobrevive a qualquer truncagem, porque `|NaN| < tau` e `False`;
* testar alcancabilidade por `Q == 0` (o sorteio) em vez de `pi == 0` (o pool)
  quebra a identidade sob delta nulo;
* um teto global em `|log2 R|` derruba a inclinacao de 0,99 para 0,52;
* sem piso no gate, 626 pares de D emitem zero contagem num gene que o
  `de_lfc_nmae` le, e `log2((0+1e-9)/ref)` vale cerca de -33.
"""

from __future__ import annotations

import warnings

import numpy as np
import pytest
from conversao import (
    PISO_LOG2,
    TETO_LOG2,
    alvo_de_composicao,
    fechar_por_lambda,
    viabilidade_do_delta,
)
from pseudobulk import ALVO_SOMA


def inclinacao(y: np.ndarray, x: np.ndarray) -> float:
    """Inclinacao da regressao de y em x, centrada -- a conta de `validar_numerico._ajustar`."""
    x = np.asarray(x, dtype=np.float64).ravel()
    y = np.asarray(y, dtype=np.float64).ravel()
    xm, ym = x.mean(), y.mean()
    vx = float(((x - xm) ** 2).mean())
    return float(((x - xm) * (y - ym)).mean() / vx) if vx > 0 else float("nan")


def fora_do_alvo(momentos, alvos: list[str]) -> np.ndarray:
    """Mascara [n_alvos, n_genes] com o gene-alvo de cada linha em False."""
    coluna = {g: j for j, g in enumerate(momentos.genes)}
    mascara = np.ones((len(alvos), len(momentos.genes)), dtype=bool)
    for i, a in enumerate(alvos):
        if a in coluna:
            mascara[i, coluna[a]] = False
    return mascara


def converter_todos(momentos, delta: np.ndarray, **kw) -> list:
    return [alvo_de_composicao(momentos, delta[i], **kw) for i in range(delta.shape[0])]


# --------------------------------------------------------------------------
# T01 / T02 / T03: finitude, identidade e fechamento
# --------------------------------------------------------------------------


def test_conversao_nunca_produz_nan_nem_negativo(momentos, delta_lfcD):
    """T01 -- em dois deltas, 5.559.900 pares cada: nada de NaN, nada negativo."""
    m = momentos("D")
    sintetico = np.random.default_rng(3).normal(0.0, 0.3, size=delta_lfcD.shape)
    for nome, d in (("lfc_D", delta_lfcD), ("N(0;0,3)", sintetico)):
        with warnings.catch_warnings():
            warnings.simplefilter("error", RuntimeWarning)
            alvos = converter_todos(m, d[:12])
        for a in alvos:
            assert np.all(np.isfinite(a.q)), nome
            assert np.all(a.q >= 0.0), nome
            assert np.all(np.isfinite(a.delta_efetivo)), nome
            assert np.isfinite(a.lambda_de_fechamento), nome


def test_conversao_de_todos_os_300_alvos_e_finita(momentos, delta_lfcD):
    """T01, varredura inteira: a propriedade tem de valer nos 300, nao numa amostra."""
    m = momentos("D")
    for a in converter_todos(m, np.asarray(delta_lfcD, dtype=np.float64)):
        assert np.all(np.isfinite(a.q)) and np.all(a.q >= 0.0)
        assert np.all(np.isfinite(a.delta_efetivo))


def test_delta_nulo_e_identidade_na_conversao(momentos):
    """T02 -- a propriedade mais carregada do desenho, no lado da conversao."""
    m = momentos("D")
    a = alvo_de_composicao(m, np.zeros(len(m.genes)))
    assert abs(a.lambda_de_fechamento) <= 1e-15
    assert float(np.abs(a.q - m.pi).max()) <= 1e-17
    assert float(np.abs(a.delta_efetivo).max()) <= 1e-15
    assert a.n_no_piso_factivel == 0
    assert a.n_no_piso_do_gate == 0
    assert a.n_no_teto_do_gate == 0


def test_fechamento_resolve_a_composicao(momentos, delta_lfcD):
    """T03 -- `sum expm1(m) == 50.000` nos 300 alvos, com lambda pequeno."""
    m = momentos("D")
    alvos = converter_todos(m, np.asarray(delta_lfcD, dtype=np.float64))
    residuos = np.array([a.residuo_de_fechamento for a in alvos])
    lambdas = np.array([a.lambda_de_fechamento for a in alvos])
    assert float(np.abs(residuos).max()) <= 1e-8
    assert float(np.sqrt((lambdas**2).mean())) == pytest.approx(1.012e-03, rel=0.25)
    assert float(np.abs(lambdas).max()) <= 0.02
    for a in alvos[:20]:
        assert float(abs(np.expm1(np.log1p(ALVO_SOMA * a.q)).sum() - ALVO_SOMA)) <= 1e-6
        assert abs(float(a.q.sum()) - 1.0) <= 1e-12


def test_funcao_de_fechamento_e_monotona(momentos, delta_lfcD):
    """T03 -- a raiz e unica porque `f` nao decresce; e o que justifica a bissecao."""
    m = momentos("D")
    m0 = np.maximum(m.bulk + np.asarray(delta_lfcD[0], dtype=np.float64), 0.0)

    def f(lam: float) -> float:
        return float(np.expm1(np.maximum(m0 + lam, 0.0)).sum() - ALVO_SOMA)

    pontos = np.linspace(-5.0, 5.0, 50)
    valores = np.array([f(p) for p in pontos])
    assert np.all(np.diff(valores) >= -1e-9)
    assert f(-40.0) == pytest.approx(-ALVO_SOMA, rel=1e-9)
    assert f(40.0) > 0
    lam, residuo = fechar_por_lambda(m0)
    assert abs(residuo) <= 1e-8
    assert abs(f(lam)) <= 1e-8


# --------------------------------------------------------------------------
# T05: o gene-alvo nao e variavel de folga
# --------------------------------------------------------------------------


def test_gene_alvo_nao_absorve_folga(momentos, ref, delta_lfcD):
    """T05 -- `pi` do gene-alvo e minusculo, logo 1/pi e uma alavanca de 16.000x."""
    import pandas as pd

    m = momentos("D")
    alvos = pd.read_csv(ref / "pert_counts.csv")["target_gene"].astype(str).tolist()
    coluna = {g: j for j, g in enumerate(m.genes)}
    pi_alvo = np.array([m.pi[coluna[a]] for a in alvos if a in coluna])
    assert pi_alvo.size == 300
    assert float(np.median(pi_alvo)) == pytest.approx(6.231e-05, rel=0.05)
    assert float(pi_alvo.max()) < 4e-03
    assert float(m.pi.max()) < 0.03
    # Um knockdown total desloca a composicao global em log2(1/(1-pi)), nao em 0,17.
    assert float(np.log2(1.0 / (1.0 - np.median(pi_alvo)))) < 2e-04

    # E a emissao trata o alvo como os outros 18.532: nenhum ramo especial.
    a = alvo_de_composicao(m, np.asarray(delta_lfcD[0], dtype=np.float64))
    j = coluna[alvos[0]]
    reconstruido = float(np.expm1(m.bulk[j] + a.delta_efetivo[j]) / ALVO_SOMA)
    assert abs(float(a.q[j]) - reconstruido) <= 1e-15


# --------------------------------------------------------------------------
# T06 / T08: os dois clamps, e o teto que NAO existe
# --------------------------------------------------------------------------


def test_sem_teto_em_log2_da_razao(momentos, ref, delta_lfcD):
    """T06 -- um teto global em `|log2 R|` leva a inclinacao de 0,99 para 0,52."""
    import pandas as pd

    m = momentos("D")
    alvos = pd.read_csv(ref / "pert_counts.csv")["target_gene"].astype(str).tolist()
    d = np.asarray(delta_lfcD, dtype=np.float64)
    usar = fora_do_alvo(m, alvos)

    # Referencia: o delta clipado na fronteira FACTIVEL -- o unico clamp legitimo.
    clipado = np.maximum(m.bulk[None, :] + d, 0.0) - m.bulk[None, :]
    efetivo = np.stack([a.delta_efetivo for a in converter_todos(m, d)])
    assert inclinacao(efetivo[usar], clipado[usar]) == pytest.approx(1.0, abs=0.03)

    # Controle: teto global de 3,0 em log2 da razao de composicao.
    with np.errstate(divide="ignore", invalid="ignore"):
        r = np.where(m.pi[None, :] > 0, np.expm1(m.bulk[None, :] + d) / m.w[None, :], 1.0)
        log2r = np.where(np.isfinite(r) & (r > 0), np.log2(np.maximum(r, 1e-300)), 0.0)
    com_teto = np.clip(log2r, -3.0, 3.0)
    delta_com_teto = np.log1p(m.w[None, :] * 2.0**com_teto) - m.bulk[None, :]
    com = inclinacao(delta_com_teto[usar], clipado[usar])
    sem = inclinacao(efetivo[usar], clipado[usar])
    # Medido: 0,7067 com teto contra 0,9978 sem. So 0,165% dos pares alcancaveis
    # batem num teto de 3,0, e eles custam 0,29 de inclinacao -- sao os pares de
    # maior |delta|, que e exatamente o que a regressao pesa.
    assert com == pytest.approx(0.707, abs=0.02)
    assert sem - com > 0.25
    assert float((np.abs(log2r) > 3.0)[np.broadcast_to(m.alcancavel, d.shape)].mean()) < 0.01
    # O maximo emitido no gate passa longe de 3,0 e nao e ele que o teto corta.
    assert float(log2r[:, m.gate].max()) == pytest.approx(2.9299, abs=1e-3)


def test_piso_do_gate_impede_log2fc_patologico(momentos, delta_lfcD):
    """T08 -- sem piso, 626 pares do gate em D emitem contagem zero num gene lido."""
    m = momentos("D")
    d = np.asarray(delta_lfcD, dtype=np.float64)

    sem_piso = converter_todos(m, d, piso_log2=np.inf, teto_log2=np.inf)
    zerados = np.stack([(a.q == 0.0) & m.gate for a in sem_piso])
    assert int(zerados.sum()) == 626
    assert int(zerados.any(axis=0).sum()) == 544
    assert int(zerados.any(axis=1).sum()) == 251

    com_piso = converter_todos(m, d)
    for a in com_piso:
        assert int(((a.q == 0.0) & m.gate).sum()) == 0
        # O piso age em `m1`, ANTES do fechamento; o lambda desloca `m` depois,
        # e num gene de `w` pequeno isso vale varios decimos em log2. Medido:
        # -4,711 no pior dos 300. Sem piso o mesmo numero vale cerca de -33.
        assert a.log2fc_min_no_gate >= -PISO_LOG2 - 0.8
        assert a.log2fc_max_no_gate <= TETO_LOG2 + 1e-9
    # O teto e inerte na medicao: o log2FC emitido maximo no gate fica em ~2,93.
    assert max(a.log2fc_max_no_gate for a in com_piso) == pytest.approx(2.9279, abs=1e-3)
    assert min(a.log2fc_min_no_gate for a in com_piso) == pytest.approx(-4.711, abs=0.01)
    assert sum(a.n_no_teto_do_gate for a in com_piso) == 0
    # O piso atua em 869 pares dos 5.559.900, e so no gate.
    assert sum(a.n_no_piso_do_gate for a in com_piso) == 869


def test_piso_mais_apertado_custa_inclinacao(momentos, ref, delta_lfcD):
    """T08 -- o piso e encolhimento declarado, e o custo dele e medido, nao suposto."""
    import pandas as pd

    m = momentos("D")
    alvos = pd.read_csv(ref / "pert_counts.csv")["target_gene"].astype(str).tolist()
    d = np.asarray(delta_lfcD[:60], dtype=np.float64)
    usar = fora_do_alvo(m, alvos[:60])
    clipado = np.maximum(m.bulk[None, :] + d, 0.0) - m.bulk[None, :]
    medido = {}
    for piso in (3.0, 4.0, 6.0):
        efetivo = np.stack([a.delta_efetivo for a in converter_todos(m, d, piso_log2=piso)])
        medido[piso] = inclinacao(efetivo[usar], clipado[usar])
    assert medido[3.0] < medido[4.0] < medido[6.0]
    assert medido[4.0] == pytest.approx(0.992, abs=0.015)


# --------------------------------------------------------------------------
# T18: a viabilidade do arquivo de delta, medida antes de emitir uma celula
# --------------------------------------------------------------------------


def test_viabilidade_do_delta_e_medida_antes_de_emitir(momentos, delta_lfcD):
    """T18 -- a quantidade que move a inclinacao e a VARIANCIA do delta, nao a massa L1."""
    m = momentos("D")
    v = viabilidade_do_delta(m, np.asarray(delta_lfcD, dtype=np.float64))
    assert v.contexto == "D"
    # `lfc_D` tem 98,9% de zeros exatos (60.298 nao nulos em 5.559.900), logo a
    # fracao sobre TODOS os pares NAO PODE passar de 0,0108. Sobre os pares
    # preditos, quase um quarto pede contagem negativa.
    assert float((np.asarray(delta_lfcD) != 0).mean()) == pytest.approx(0.01085, abs=1e-4)
    assert v.fracao_de_pares_infactiveis == pytest.approx(0.002595, abs=2e-5)
    assert v.fracao_dos_nao_nulos_infactiveis == pytest.approx(0.2393, abs=0.002)
    assert v.genes_atingidos == 7571
    assert v.variancia_perdida == pytest.approx(0.2527, abs=0.004)
    assert v.inclinacao_contra_delta_cru == pytest.approx(0.7692, abs=0.002)
    # A parcela em genes inalcancaveis com delta positivo: 7,04% da variancia.
    assert v.massa_em_genes_inalcancaveis == pytest.approx(0.0704, abs=0.004)
    # `variancia_perdida` NAO e o complemento da inclinacao -- clipar nao e
    # projetar na direcao do delta. Medido, as duas diferem por 0,022.
    assert abs((1.0 - v.variancia_perdida) - v.inclinacao_contra_delta_cru) < 0.03
    # O numero decisivo: 0,769 esta ABAIXO do piso 0,85 de `validar_numerico.FAIXA`.
    assert v.inclinacao_contra_delta_cru < 0.85


def test_viabilidade_de_um_delta_factivel_e_quase_perfeita(momentos, delta_lfcD):
    """T17/T18 -- o conserto e do preditor: gravar `delta >= -b_ctrl` da +0,9997."""
    m = momentos("D")
    cru = np.asarray(delta_lfcD, dtype=np.float64)
    factivel = np.maximum(cru, -m.bulk[None, :])
    v = viabilidade_do_delta(m, factivel)
    assert v.fracao_de_pares_infactiveis == pytest.approx(0.0, abs=1e-9)
    assert v.variancia_perdida == pytest.approx(0.0, abs=1e-9)
    assert v.inclinacao_contra_delta_cru == pytest.approx(1.0, abs=1e-9)


def test_l1_e_variancia_quase_coincidem_neste_delta(momentos, delta_lfcD):
    """T18, controle que FALSIFICA a premissa de que L1 esconderia a parcela.

    Um dos projetos descartados media a massa inalcancavel em L1, achava
    desprezivel e cortava o ramo aditivo. Medi as duas fracoes e elas quase
    coincidem: 0,0772 em L1 contra 0,0704 em variancia. Nenhuma das duas esconde
    nada, e o argumento contra o ramo aditivo nao vem de escolher a norma.

    O que decide o ramo aditivo e o efeito na INCLINACAO, e isso se mede
    emitindo -- e o `test_genes_inalcancaveis_recebem_massa_pelo_ramo_aditivo`
    de `test_kernel.py`, nao aqui.
    """
    m = momentos("D")
    d = np.asarray(delta_lfcD, dtype=np.float64)
    inalcancavel = np.broadcast_to(~m.alcancavel, d.shape) & (d > 0)
    em_l1 = float(np.abs(d[inalcancavel]).sum() / np.abs(d).sum())
    v = viabilidade_do_delta(m, d)
    assert em_l1 == pytest.approx(0.0772, abs=0.002)
    assert v.massa_em_genes_inalcancaveis == pytest.approx(0.0704, abs=0.004)
    assert abs(em_l1 - v.massa_em_genes_inalcancaveis) < 0.02
    assert int(inalcancavel.sum()) == 4758
