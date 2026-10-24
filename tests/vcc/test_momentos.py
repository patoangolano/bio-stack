"""Momentos do pool de controle: o denominador de todo o resto.

Os valores esperados aqui foram medidos contra os `context_<c>.h5ad` reais e
contra o `cell-eval2 0.19.0` instalado. Dois sao identidades, nao medidas, e e
por isso que entram como igualdade e nao como tolerancia larga:

* `b == bulk_lognorm` do pool inteiro -- o mesmo caminho de contas, bit a bit;
* `sum_g A[g] == 1e6` -- cada celula contribui exatamente 1e6 de CPM, logo a
  media aritmetica entre celulas tambem soma 1e6. O PASSO 4 de um dos projetos
  descartados corrigia um "fator global" que esta identidade mostra ser 1.
"""

from __future__ import annotations

import numpy as np
import pytest
from momentos import PISO_CPM_DO_GATE, carregar_pool, momentos_do_controle
from pseudobulk import ALVO_SOMA, bulk_lognorm

# Medidos nos tres contextos. O gate pelo CPM PONDERADO (do pseudobulk) e a
# comparacao que mostra que o gate do cell-eval2 nao e o do pseudobulk.
MEDIDO = {
    "D": {"gate": 9633, "gate_ponderado": 9565, "pi_zero": 2975, "nnz": 104_487_981},
    "E": {"gate": 9743, "gate_ponderado": 9744, "pi_zero": 2709, "nnz": 101_608_841},
    "F": {"gate": 9537, "gate_ponderado": 9499, "pi_zero": 2162, "nnz": 100_855_632},
}
CONTEXTOS = ["D", "E", "F"]


@pytest.mark.parametrize("c", CONTEXTOS)
def test_bulk_e_identico_ao_pseudobulk_do_pool_inteiro(momentos, c):
    """`b = log1p(50_000 * pi)` tem de ser o `bulk_lognorm` do pool, bit a bit."""
    m = momentos(c)
    esperado = bulk_lognorm(m.soma_por_gene[None, :])[0]
    assert np.array_equal(m.bulk, esperado)
    assert m.w.shape == (len(m.genes),)
    assert abs(float(np.expm1(m.bulk).sum()) - ALVO_SOMA) <= 1e-8


@pytest.mark.parametrize("c", CONTEXTOS)
def test_cpm_medio_por_celula_soma_um_milhao(momentos, c):
    """Identidade: nao existe fator global a corrigir em `de_lfc_nmae`."""
    m = momentos(c)
    assert abs(float(m.cpm_medio_por_celula.sum()) - 1e6) <= 1e-5
    assert np.all(np.isfinite(m.cpm_medio_por_celula))


@pytest.mark.parametrize("c", CONTEXTOS)
def test_gate_reproduz_as_contagens_medidas(momentos, c):
    m = momentos(c)
    esperado = MEDIDO[c]
    assert int(m.gate.sum()) == esperado["gate"]
    assert int((~m.alcancavel).sum()) == esperado["pi_zero"]
    # Nenhum gene inalcancavel esta no gate: o ramo aditivo nunca cria gene lido.
    assert int((m.gate & ~m.alcancavel).sum()) == 0
    assert m.nnz == esperado["nnz"]


@pytest.mark.parametrize("c", CONTEXTOS)
def test_gate_e_estritamente_maior_e_nao_e_o_do_pseudobulk(momentos, c):
    """`de_compute.py:1257` usa `>` sobre a media ARITMETICA de CPM por celula."""
    m = momentos(c)
    a = m.cpm_medio_por_celula
    assert not np.any(a == PISO_CPM_DO_GATE), "gene exatamente em 5,0: > e >= divergem"
    assert np.array_equal(m.gate, a > PISO_CPM_DO_GATE)
    ponderado = 1e6 * m.pi > PISO_CPM_DO_GATE
    assert int(ponderado.sum()) == MEDIDO[c]["gate_ponderado"]
    diferenca = int((m.gate != ponderado).sum())
    assert diferenca == {"D": 82, "E": 29, "F": 48}[c]


def _cpm_medio_em_float32(caminho) -> np.ndarray:
    """O MESMO calculo de `A`, acumulado em float32: o controle negativo.

    Usa `momentos._blocos` de proposito -- o ponto do controle e variar so o
    dtype do acumulador, mantendo tudo o mais identico.
    """
    import scipy.sparse as sp
    from momentos import _abrir, _blocos

    f, x, obs, genes = _abrir(caminho)
    acumulado = np.zeros(len(genes), dtype=np.float32)
    with f:
        for _, bloco in _blocos(x, len(obs), 2048, len(genes)):
            lib = np.asarray(bloco.sum(axis=1), dtype=np.float64).ravel()
            escala = np.zeros_like(lib)
            np.divide(1e6, lib, out=escala, where=lib > 0)
            em32 = sp.csr_matrix(sp.diags(escala.astype(np.float32)) @ bloco.astype(np.float32))
            acumulado += np.asarray(em32.sum(axis=0), dtype=np.float32).ravel()
    return acumulado / np.float32(len(obs))


def test_recast_para_float64_e_obrigatorio(momentos, ref):
    """Controle: acumular `A` em float32 erra 4e8 vezes mais do que em float64.

    Medido em D: `sum_g A - 1e6` vale -1,16e-10 acumulando em float64 e +4,87e-02
    em float32; o erro relativo por gene chega a 8,79e-07 e o rms absoluto a
    5,54e-05. `A` e o denominador do gate E do log2FC que `de_lfc_nmae` le, e
    5,5e-05 de rms num CPM de 5,0 e 1,1e-05 relativo no limiar.

    O que este teste NAO afirma, porque medi e e falso: que o gate vire. Nos
    tres contextos o conjunto `A > 5,0` e o mesmo nos dois dtypes (9.633 em D,
    zero genes de diferenca). Nenhum gene real passa perto o bastante de 5,0.
    O recast e obrigatorio pelo LFC, nao pela filtragem.
    """
    m = momentos("D")
    em_float32 = _cpm_medio_em_float32(ref / "context_D.h5ad")
    erro64 = abs(float(m.cpm_medio_por_celula.sum()) - 1e6)
    erro32 = abs(float(np.float64(em_float32).sum()) - 1e6)
    assert erro64 < 1e-8
    assert erro32 > 1e-2
    assert erro32 / erro64 > 1e6
    rms = float(np.sqrt(np.mean((np.float64(em_float32) - m.cpm_medio_por_celula) ** 2)))
    assert rms == pytest.approx(5.54e-05, rel=0.1)
    assert np.array_equal(em_float32 > np.float32(5.0), m.gate)


@pytest.mark.parametrize("c", CONTEXTOS)
def test_bibliotecas_e_lotes(momentos, c):
    m = momentos(c)
    assert m.n_celulas == 18_400
    assert m.bibliotecas.shape == (18_400,)
    assert m.bibliotecas.dtype == np.float64
    assert np.all(m.bibliotecas > 0)
    assert int(np.unique(m.lote).size) == 46
    assert float(np.median(m.bibliotecas)) == pytest.approx(
        {"D": 20_235.5, "E": 19_681.0, "F": 19_845.0}[c]
    )
    # O teto de contagem por celula do contrato e 1e6; a maior biblioteca e 5e4.
    assert float(m.bibliotecas.max()) == pytest.approx(
        {"D": 53_865.0, "E": 51_808.0, "F": 52_365.0}[c]
    )


def test_zerofrac_e_contagem_maxima(momentos):
    m = momentos("D")
    assert m.zerofrac.shape == (len(m.genes),)
    assert np.all((m.zerofrac >= 0.0) & (m.zerofrac <= 1.0))
    # Gene inalcancavel esta zerado em todas as celulas, por definicao.
    assert np.all(m.zerofrac[~m.alcancavel] == 1.0)
    assert float(m.zerofrac[m.gate].mean()) == pytest.approx(0.4171, abs=5e-4)
    assert float(m.nnz) / m.n_celulas == pytest.approx(5678.69, abs=0.01)
    # A contagem maxima decide o dtype do pool: uint16 cabe ate 65.535.
    assert 0 < m.contagem_maxima < 65_535


def test_bloco_nao_muda_os_momentos(controle_pequeno):
    """Tamanho de bloco e so memoria: dois tamanhos primos entre si dao o mesmo."""
    a = momentos_do_controle(controle_pequeno, linhas_por_bloco=2048)
    b = momentos_do_controle(controle_pequeno, linhas_por_bloco=7)
    assert np.array_equal(a.soma_por_gene, b.soma_por_gene)
    assert np.array_equal(a.bibliotecas, b.bibliotecas)
    assert np.array_equal(a.gate, b.gate)
    assert np.array_equal(a.zerofrac, b.zerofrac)
    assert a.cpm_medio_por_celula == pytest.approx(b.cpm_medio_por_celula, rel=0, abs=1e-9)
    assert a.nnz == b.nnz


def test_pool_e_uint16_e_bate_com_as_somas(controle_pequeno):
    p = carregar_pool(controle_pequeno)
    m = momentos_do_controle(controle_pequeno)
    assert p.shape == (m.n_celulas, len(m.genes))
    assert p.dtype == np.uint16
    assert p.indices.dtype == np.int32
    assert p.nnz == m.nnz
    assert np.array_equal(np.asarray(p.sum(axis=0), dtype=np.float64).ravel(), m.soma_por_gene)
    assert np.array_equal(np.asarray(p.sum(axis=1), dtype=np.float64).ravel(), m.bibliotecas)
    assert p.has_sorted_indices
