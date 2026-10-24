"""O kernel de emissao: do alvo de composicao as contagens inteiras por celula.

A propriedade que carrega o desenho e a de T02: sob delta nulo a emissao devolve
as doadoras BIT A BIT. Ela nao e enfeite -- e o que garante ao mesmo tempo que
(i) nenhum vies entra pelo caminho e (ii) o ruido de amostragem das 400
doadoras fica DENTRO de `B_pred`, que e o numerador da razao que decide o fator
`r` do `expr_mse`. Qualquer peca que "melhore" a emissao e quebre essa
identidade esta mexendo nos dois ao mesmo tempo.

Os controles negativos aqui nao sao hipoteticos: cada um e um caminho que um dos
tres projetos de emissao propunha.
"""

from __future__ import annotations

from typing import Any, cast

import numpy as np
import pytest
import scipy.sparse as sp
from conversao import alvo_de_composicao
from kernel import (
    N_CELULAS,
    distribuir_por_celula,
    emitir_grupo,
    inteirizar_por_coluna,
    massa_por_gene,
    semente_do_grupo,
    sortear_doadoras,
)


def densas(pool, indices: np.ndarray) -> np.ndarray:
    return np.asarray(pool[indices].todense(), dtype=np.float64)


def inclinacao(y: np.ndarray, x: np.ndarray) -> float:
    x = np.asarray(x, dtype=np.float64).ravel()
    y = np.asarray(y, dtype=np.float64).ravel()
    xm, ym = x.mean(), y.mean()
    vx = float(((x - xm) ** 2).mean())
    return float(((x - xm) * (y - ym)).mean() / vx) if vx > 0 else float("nan")


# --------------------------------------------------------------------------
# T02: identidade sob delta nulo
# --------------------------------------------------------------------------


def test_delta_nulo_e_identidade_bit_a_bit(momentos, pool, ref):
    """T02 -- com `delta == 0` a emissao devolve exatamente as celulas sorteadas."""
    import pandas as pd

    m = momentos("D")
    p = pool("D")
    alvo = pd.read_csv(ref / "pert_counts.csv")["target_gene"].astype(str).tolist()[0]
    g = emitir_grupo(m, p, np.zeros(len(m.genes)), alvo, semente_base=0)

    rng = np.random.default_rng(semente_do_grupo(0, "D", alvo))
    esperado = sp.csr_matrix(densas(p, sortear_doadoras(rng, m.lote)).astype(np.float32))
    emitido = g.matriz
    assert np.array_equal(emitido.indptr, esperado.indptr)
    assert np.array_equal(emitido.indices, esperado.indices)
    assert np.array_equal(emitido.data, esperado.data)
    assert g.n_genes_criados == 0
    assert np.all(g.alvo_de_composicao.razao == 1.0)

    # E o `residuo_bulk_rms` NAO e zero, e nao deve ser: a identidade e nas
    # CONTAGENS, nao no pseudobulk. 400 celulas nao sao as 18.400 do pool, entao
    # `bulk_lognorm` das 400 fica a 0,028 rms do pool. Esse residuo e ruido de
    # amostragem puro e e exatamente o que mantem o termo de amostragem dentro de
    # `B_pred` -- zera-lo seria pinar o agregado do grupo e derrubar a razao do
    # `expr_mse` para a faixa que o aviso de 0,7 existe para pegar.
    assert g.residuo_bulk_rms == pytest.approx(0.0277, abs=0.004)
    assert g.residuo_bulk_rms_no_gate == pytest.approx(0.0, abs=0.08)


def test_massa_por_gene_e_identidade_sob_delta_nulo(momentos, pool):
    """T02 -- `M == Q` exige alcancabilidade por `pi` do POOL, nao por `Q` do sorteio."""
    m = momentos("D")
    p = pool("D")
    rng = np.random.default_rng(semente_do_grupo(0, "D", "ABCC1"))
    s = sortear_doadoras(rng, m.lote)
    d = densas(p, s)
    somas, total = d.sum(axis=0), float(d.sum())
    a = alvo_de_composicao(m, np.zeros(len(m.genes)))

    massa = massa_por_gene(a.q, a.razao, m.alcancavel, somas, total)
    # Exato, nao aproximado: a razao e 1,0 em todos os genes sob delta nulo.
    assert np.array_equal(massa, somas)
    assert np.all(a.razao == 1.0)
    assert float(np.abs(massa[~m.alcancavel]).max()) == 0.0

    # Controle: alcancabilidade pelo SORTEIO. Genes presentes no pool mas ausentes
    # das 400 celulas sorteadas recebem `S0*q` do ar e a identidade quebra.
    por_sorteio = massa_por_gene(a.q, a.razao, somas > 0, somas, total)
    assert float(np.abs(por_sorteio - somas).max()) > 1.0

    # Controle: a forma `Q*q/pi`, com dois arredondamentos independentes.
    ingenuo = somas * a.q / np.where(m.alcancavel, m.pi, 1.0)
    assert not np.array_equal(ingenuo, somas)
    assert float(np.abs(ingenuo / np.maximum(somas, 1e-30) - 1.0).max()) > 0.0


# --------------------------------------------------------------------------
# T07: o ramo aditivo e o que recupera os genes inalcancaveis
# --------------------------------------------------------------------------


def test_genes_inalcancaveis_recebem_massa_pelo_ramo_aditivo(momentos, pool):
    """T07 -- sem o ramo aditivo, delta positivo em gene de `pi == 0` nao sai."""
    m = momentos("D")
    p = pool("D")
    inalcancaveis = np.flatnonzero(~m.alcancavel)[:50]
    delta = np.zeros(len(m.genes))
    delta[inalcancaveis] = 1.5

    g = emitir_grupo(m, p, delta, "ABCC1", semente_base=0)
    soma_emitida = np.asarray(g.matriz.sum(axis=0), dtype=np.float64).ravel()
    assert np.all(soma_emitida[inalcancaveis] >= 1.0)
    assert g.n_genes_criados >= 50

    # Controle: multiplicacao pura. `M = Q*q/pi` e 0/0 nos inalcancaveis, e quem
    # trata isso como zero nunca cria contagem onde a doadora nao tinha nenhuma.
    rng = np.random.default_rng(semente_do_grupo(0, "D", "ABCC1"))
    s = sortear_doadoras(rng, m.lote)
    d = densas(p, s)
    somas = d.sum(axis=0)
    a = alvo_de_composicao(m, delta)
    multiplicativo = np.where(m.alcancavel, somas * a.q / np.where(m.pi > 0, m.pi, 1.0), 0.0)
    assert float(multiplicativo[inalcancaveis].max()) == 0.0


def test_ramo_aditivo_vale_sete_centesimos_de_inclinacao(momentos, ref, delta_lfcD):
    """T07 -- quanto custa descartar os genes de `pi == 0`, medido na inclinacao.

    A conta e analitica, no nivel de `q`: o kernel apenas realiza `q` com ruido
    de amostragem, entao e `q` que fixa a inclinacao. Medido: 0,9978 com o ramo
    aditivo e 0,8928 sem -- 0,105 de inclinacao, mais do que a parcela de 7,04%
    da variancia que vive nesses genes, porque eles concentram |delta| grande e
    a regressao pesa justamente por isso.
    """
    import pandas as pd

    m = momentos("D")
    alvos = pd.read_csv(ref / "pert_counts.csv")["target_gene"].astype(str).tolist()
    d = np.asarray(delta_lfcD, dtype=np.float64)
    coluna = {g: j for j, g in enumerate(m.genes)}
    usar = np.ones(d.shape, dtype=bool)
    for i, a in enumerate(alvos):
        if a in coluna:
            usar[i, coluna[a]] = False

    clipado = np.maximum(m.bulk[None, :] + d, 0.0) - m.bulk[None, :]
    efetivo = np.stack([alvo_de_composicao(m, d[i]).delta_efetivo for i in range(d.shape[0])])
    sem_ramo = np.where(m.alcancavel[None, :], efetivo, 0.0)

    com = inclinacao(efetivo[usar], clipado[usar])
    sem = inclinacao(sem_ramo[usar], clipado[usar])
    assert com == pytest.approx(0.9978, abs=0.003)
    assert sem == pytest.approx(0.8928, abs=0.004)
    assert com - sem > 0.08


# --------------------------------------------------------------------------
# T09 / T10: a massa fecha, e a inteirizacao nao enviesa
# --------------------------------------------------------------------------


def test_massa_por_coluna_dentro_de_uma_contagem(momentos, pool, delta_lfcD):
    """T09 -- o fracionario fecha exato e o inteiro fica a menos de uma contagem."""
    m = momentos("D")
    p = pool("D")
    delta = np.asarray(delta_lfcD[0], dtype=np.float64)
    rng = np.random.default_rng(semente_do_grupo(0, "D", "ABCC1"))
    s = sortear_doadoras(rng, m.lote)
    d = densas(p, s)
    somas, total = d.sum(axis=0), float(d.sum())
    a = alvo_de_composicao(m, delta)
    massa = massa_por_gene(a.q, a.razao, m.alcancavel, somas, total)

    y = distribuir_por_celula(d, d.sum(axis=1), massa, somas, m.alcancavel)
    assert float(np.abs(y.sum(axis=0) - massa).max()) <= 1e-9
    assert np.all(y >= 0.0)

    z = inteirizar_por_coluna(y, rng)
    erro = z.sum(axis=0) - massa
    assert float(np.abs(erro).max()) <= 1.0
    assert float(np.sqrt((erro**2).mean())) <= 0.35
    assert z.dtype == np.int64
    assert np.all(z >= 0)


def test_distribuicao_e_continua_na_fronteira(momentos, pool):
    """T09 -- em `M == Q` o incremento e zero e os dois ramos coincidem em `Y = D`."""
    m = momentos("D")
    p = pool("D")
    rng = np.random.default_rng(semente_do_grupo(0, "D", "ACLY"))
    d = densas(p, sortear_doadoras(rng, m.lote))
    somas = d.sum(axis=0)
    y = distribuir_por_celula(d, d.sum(axis=1), somas.copy(), somas, m.alcancavel)
    assert np.array_equal(y, d)


def test_inteirizacao_sem_vies(momentos):
    """T10 -- 500 sementes: o PPS sistematico erra 2,7e-05 relativo; `np.round` erra 0,6%."""
    rng_base = np.random.default_rng(5)
    doadora = rng_base.poisson(3.0, size=N_CELULAS).astype(np.float64)
    y = (1.3 * doadora)[:, None]
    massa = float(y.sum())

    somas = np.array(
        [float(inteirizar_por_coluna(y, np.random.default_rng(s)).sum()) for s in range(500)]
    )
    assert abs(float(somas.mean()) - massa) / massa <= 1e-3
    assert float(np.abs(somas - massa).max()) <= 1.0

    # Controle: arredondamento determinstico desloca a media em 0,6%, e e
    # justamente a media do grupo que `de_lfc_nmae` e `expr_mse` medem.
    arredondado = float(np.round(y).sum())
    assert abs(arredondado - massa) / massa == pytest.approx(0.0060, abs=0.002)

    # `frac == 0` tem de ser identidade, qualquer que seja o deslocamento.
    inteiro = doadora[:, None]
    for s in (0, 1, 7):
        assert np.array_equal(inteirizar_por_coluna(inteiro, np.random.default_rng(s)), inteiro)


def test_inteirizacao_nunca_passa_de_um_por_celula(momentos):
    """T10 -- `picks` e 0 ou 1 por construcao: diff de piso de uma soma crescente."""
    rng = np.random.default_rng(13)
    y = rng.random((N_CELULAS, 64)) * 4.0
    z = inteirizar_por_coluna(y, rng)
    assert np.all((z - np.floor(y)) >= 0)
    assert np.all((z - np.floor(y)) <= 1)


# --------------------------------------------------------------------------
# T11 / T12 / T13: dtype, guardas e CSR canonico
# --------------------------------------------------------------------------


def test_kernel_aceita_o_float32_do_h5ad(momentos, pool, delta_lfcD):
    """T11 -- nada no caminho pode pedir inteiro de um array float."""
    m = momentos("D")
    p = pool("D")
    delta = np.asarray(delta_lfcD[1], dtype=np.float64)
    como_esta = emitir_grupo(m, p, delta, "ACLY", semente_base=0)
    como_float32 = emitir_grupo(
        m, sp.csr_matrix((p.data.astype(np.float32), p.indices, p.indptr), shape=p.shape),
        delta, "ACLY", semente_base=0,
    )
    assert np.array_equal(como_esta.matriz.data, como_float32.matriz.data)
    assert np.array_equal(como_esta.matriz.indices, como_float32.matriz.indices)

    # Controle explicito: e assim que o caminho ingenuo morre. O `cast` existe
    # porque o basedpyright TAMBEM recusa isto estaticamente -- o que e evidencia
    # a mais, nao a menos; o teste afirma que falha em tempo de execucao.
    em_float32 = cast("Any", np.array([3, 7, 0], dtype=np.float32))
    with pytest.raises(TypeError, match="safe"):
        np.random.default_rng(0).binomial(em_float32, 0.5)


@pytest.mark.parametrize("i", [0, 1, 2, 3, 4])
def test_guardas_do_grupo(momentos, pool, ref, delta_lfcD, i):
    """T12 -- as guardas sao medidas sobre o `Z` real, nunca reconstruidas."""
    import pandas as pd

    m = momentos("D")
    p = pool("D")
    alvos = pd.read_csv(ref / "pert_counts.csv")["target_gene"].astype(str).tolist()
    g = emitir_grupo(m, p, np.asarray(delta_lfcD[i], dtype=np.float64), alvos[i], semente_base=0)
    denso = np.asarray(g.matriz.todense())
    assert g.linhas_nulas == 0
    assert g.soma_minima_por_celula > 0
    assert g.soma_maxima_por_celula <= 100_000  # o teto do contrato e 1e6
    assert np.all(denso >= 0)
    assert np.array_equal(denso, np.floor(denso))
    assert np.all(np.isfinite(denso))
    assert g.matriz.shape == (N_CELULAS, len(m.genes))
    # Todo gene do gate com `q > 0` sai com pelo menos uma contagem.
    soma = denso.sum(axis=0)
    do_gate = m.gate & (g.alvo_de_composicao.q > 0)
    assert int((soma[do_gate] < 1).sum()) == 0
    # O residuo contra o pretendido e ruido de amostragem das 400, nao vies.
    assert 0.01 < g.residuo_bulk_rms < 0.08


def test_csr_canonico_do_grupo(momentos, pool, delta_lfcD):
    """T13 -- indices estritamente crescentes, sem repetido e sem zero armazenado."""
    m = momentos("D")
    p = pool("D")
    g = emitir_grupo(m, p, np.asarray(delta_lfcD[2], dtype=np.float64), "ADNP", semente_base=0)
    x = g.matriz
    assert x.has_sorted_indices
    assert int(x.indptr[-1]) == x.data.size == x.indices.size
    assert not np.any(x.data == 0)
    for linha in range(0, N_CELULAS, 37):
        faixa = x.indices[x.indptr[linha] : x.indptr[linha + 1]]
        assert np.all(np.diff(faixa) > 0)
        assert np.unique(faixa).size == faixa.size
    assert g.nnz == x.nnz


# --------------------------------------------------------------------------
# T21: determinismo e regeneracao por grupo
# --------------------------------------------------------------------------


def test_semente_do_grupo_e_estavel_e_distinta():
    """T21 -- derivada por hash do rotulo, nao por consumo sequencial de um rng."""
    a = semente_do_grupo(0, "D", "ABCC1")
    b = semente_do_grupo(0, "D", "ABCC1")
    assert a.spawn_key == b.spawn_key and a.entropy == b.entropy
    chaves = {
        semente_do_grupo(0, c, alvo).spawn_key
        for c in ("D", "E", "F")
        for alvo in ("ABCC1", "ACLY", "ADNP")
    }
    assert len(chaves) == 9
    assert semente_do_grupo(0, "D", "ABCC1").spawn_key != semente_do_grupo(1, "D", "ABCC1").entropy
    # blake2b truncado em 63 bits: cabe em int64 sem sinal negativo.
    for parte in semente_do_grupo(0, "D", "ABCC1").spawn_key:
        assert 0 <= parte < 2**63


def test_regeneracao_de_um_grupo_nao_depende_dos_outros(momentos, pool, delta_lfcD):
    """T21 -- reemitir so um grupo da o mesmo resultado, em qualquer ordem."""
    m = momentos("D")
    p = pool("D")
    d = np.asarray(delta_lfcD, dtype=np.float64)
    direto = emitir_grupo(m, p, d[3], "AKT1", semente_base=0)
    # Emite outros dois antes, consumindo aleatoriedade, e repete.
    emitir_grupo(m, p, d[0], "ABCC1", semente_base=0)
    emitir_grupo(m, p, d[1], "ACLY", semente_base=0)
    depois = emitir_grupo(m, p, d[3], "AKT1", semente_base=0)
    assert np.array_equal(direto.matriz.data, depois.matriz.data)
    assert np.array_equal(direto.matriz.indices, depois.matriz.indices)
    assert direto.semente == depois.semente


def test_doadoras_sao_estratificadas_e_com_reposicao(momentos):
    """As 400 saem dos 46 `ntc_id` em 8 ou 9 por estrato, com reposicao."""
    m = momentos("D")
    rng = np.random.default_rng(semente_do_grupo(0, "D", "ABCC1"))
    s = sortear_doadoras(rng, m.lote)
    assert s.size == N_CELULAS
    assert s.dtype == np.int64
    assert np.all((s >= 0) & (s < m.n_celulas))
    contagem = np.bincount(m.lote[s], minlength=46)
    assert contagem.size == 46
    assert set(np.unique(contagem).tolist()) <= {8, 9}
    assert int(contagem.sum()) == N_CELULAS
    assert int((contagem == 9).sum()) == 32 and int((contagem == 8).sum()) == 14
    # Com reposicao: em 400 sorteios de 18.400 com estratos de 400, repeticao
    # e esperada em alguma das sementes. Nao e defeito, e o que justifica
    # `tr(Sigma)/n` como correcao sem fator de populacao finita.
    repetidos = sum(
        int(N_CELULAS - np.unique(sortear_doadoras(np.random.default_rng(k), m.lote)).size)
        for k in range(20)
    )
    assert repetidos > 0


def test_gene_alcancavel_sem_doadora_emite_zero(momentos, pool):
    """A decisao sobre `Q == 0` num gene que o pool TEM: emite zero, e por que isso serve.

    `M[g] = Q[g] * razao[g]` da zero quando as 400 doadoras sorteadas nao tem
    nenhuma contagem daquele gene, qualquer que seja o delta. Isso e deliberado e
    e o que faz a identidade sob delta nulo fechar -- mas e uma decisao, nao um
    acidente, e o que a torna defensavel e uma medida: num sorteio real de D,
    3.128 dos 15.558 genes alcancaveis tem `Q == 0`, e NENHUM deles esta no gate
    (`min Q` no gate medido 23). Nenhum gene que o placar le e afetado.

    Sem este comportamento, a massa `S0*q` desses genes criaria 692 contagens
    fantasma em 406 genes num unico grupo sob delta NULO.
    """
    m = momentos("D")
    p = pool("D")
    rng = np.random.default_rng(semente_do_grupo(0, "D", "ABCC1"))
    d = densas(p, sortear_doadoras(rng, m.lote))
    somas = d.sum(axis=0)
    sem_doadora = m.alcancavel & (somas == 0)
    assert int(sem_doadora.sum()) > 1000
    assert int((m.gate & (somas == 0)).sum()) == 0
    assert float(somas[m.gate].min()) >= 20.0

    # Delta positivo grande nesses genes, e a emissao continua zero.
    delta = np.zeros(len(m.genes))
    delta[sem_doadora] = 2.0
    g = emitir_grupo(m, p, delta, "ABCC1", semente_base=0)
    emitido = np.asarray(g.matriz.sum(axis=0), dtype=np.float64).ravel()
    assert float(emitido[sem_doadora].max()) == 0.0
    # E os genes de `pi == 0` com o mesmo delta SAEM, pelo ramo aditivo: a
    # diferenca entre os dois casos e o pool, nao o sorteio.
    delta2 = np.zeros(len(m.genes))
    inalcancaveis = np.flatnonzero(~m.alcancavel)[:30]
    delta2[inalcancaveis] = 2.0
    g2 = emitir_grupo(m, p, delta2, "ABCC1", semente_base=0)
    emitido2 = np.asarray(g2.matriz.sum(axis=0), dtype=np.float64).ravel()
    assert np.all(emitido2[inalcancaveis] >= 1.0)


def test_massa_fantasma_de_um_ulp_na_razao_seria_catastrofica(momentos, pool):
    """Por que a razao tem de ser 1,0 EXATO, e nao 1 +- um ulp.

    `floor(D * (1 - 2.22e-16))` e `D - 1` para todo `D` nao nulo. Um ulp abaixo
    de 1 na razao decrementa o piso de cada entrada do grupo e joga 2,27e6
    contagens na inteirizacao estocastica, que as devolve em media mas nao bit a
    bit. E assim que a identidade sob delta nulo morre sem nenhum aviso.
    """
    m = momentos("D")
    p = pool("D")
    rng = np.random.default_rng(semente_do_grupo(0, "D", "ABCC1"))
    d = densas(p, sortear_doadoras(rng, m.lote))
    um_ulp_abaixo = np.nextafter(1.0, 0.0)
    assert um_ulp_abaixo < 1.0
    piso_certo = np.floor(d * 1.0)
    piso_errado = np.floor(d * um_ulp_abaixo)
    assert np.array_equal(piso_certo, d)
    assert not np.array_equal(piso_errado, d)
    assert int((piso_errado != d).sum()) > 2_000_000
