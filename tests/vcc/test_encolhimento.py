"""O agrupamento do fator por gene: forma sem nivel, e nada imputado.

`encolhimento.py` mede `k_g = sd_real / sd_previsto` em 16 direcoes de doador e
agrupa. O que estes testes fixam e o que o modulo PROMETE no docstring e o que a
aplicacao em `varredura_k.py` depende:

* o peso e FORMA — media geometrica 1 sobre os genes medidos, entao agrupar nao
  pode mover o nivel, que continua sendo `k_desvio`;
* gene sem medida sai em 1,0 exato, nao imputado: 37% do painel esta nesse caso
  e tratar ausencia como estimativa inventaria sinal em 12.243 genes;
* o nivel de cada direcao e dividido fora ANTES do agrupamento, porque o global
  medido varia de 3,40 a 8,69 entre direcoes e misturar niveis contaminaria a
  forma com qual direcao calhou de medir o gene.
"""

from __future__ import annotations

import encolhimento as enc
import numpy as np
import pytest

GENES_PAINEL = [f"G{i}" for i in range(8)]


def medida(
    rotulo: str,
    genes: list[str],
    k_relativo: dict[str, float],
    k_global: float = 5.0,
    entradas: float = 1000.0,
) -> enc.MedidaDeDirecao:
    """Uma direcao sintetica em que o gene `g` tem `k_g = k_global * k_relativo[g]`.

    `sq_real / sq_previsto = k_g**2`, entao fixar `sq_previsto = 1` e
    `sq_real = k_g**2` produz exatamente a razao pedida.
    """
    k = np.array([k_global * k_relativo.get(g, 1.0) for g in genes], dtype=np.float64)
    n = np.array([entradas if g in k_relativo else entradas for g in genes], dtype=np.float64)
    return enc.MedidaDeDirecao(
        rotulo=rotulo,
        genes=genes,
        entradas=n,
        sq_previsto=np.ones(len(genes), dtype=np.float64),
        sq_real=np.square(k),
        sq_previsto_diagonal=np.ones(len(genes), dtype=np.float64),
        k_global=k_global,
    )


def test_peso_tem_media_geometrica_um_sobre_os_medidos():
    """Agrupar nao move o nivel: o nivel e `k_desvio` e a forma e neutra."""
    m = medida("a->b", GENES_PAINEL[:6], {"G0": 3.0, "G1": 2.0, "G2": 0.5, "G3": 0.25})
    peso, entradas, _ = enc.agrupar([m], GENES_PAINEL)
    medidos = entradas > 0
    assert np.exp(np.log(peso[medidos]).mean()) == pytest.approx(1.0, abs=1e-12)


def test_gene_sem_medida_sai_em_um_exato():
    """G6 e G7 nao estao em nenhuma direcao: herdam o escalar, nao uma imputacao."""
    m = medida("a->b", GENES_PAINEL[:6], {"G0": 4.0})
    peso, entradas, diagnostico = enc.agrupar([m], GENES_PAINEL)
    assert entradas[6] == 0.0
    assert entradas[7] == 0.0
    assert peso[6] == 1.0
    assert peso[7] == 1.0
    assert diagnostico["genes_medidos"] == 6.0


def test_o_nivel_da_direcao_e_dividido_fora():
    """Duas direcoes com a MESMA forma e niveis diferentes dao o mesmo peso.

    Se o nivel nao fosse dividido fora, a direcao de `k_global` maior dominaria
    e a forma viria contaminada por qual direcao calhou de medir o gene.
    """
    forma = {"G0": 3.0, "G1": 0.4, "G2": 1.5}
    baixo = medida("a->b", GENES_PAINEL[:6], forma, k_global=3.4)
    alto = medida("c->d", GENES_PAINEL[:6], forma, k_global=8.7)
    peso_baixo, _, _ = enc.agrupar([baixo], GENES_PAINEL)
    peso_alto, _, _ = enc.agrupar([alto], GENES_PAINEL)
    peso_juntos, _, _ = enc.agrupar([baixo, alto], GENES_PAINEL)
    assert np.allclose(peso_baixo, peso_alto, atol=1e-12)
    assert np.allclose(peso_juntos, peso_baixo, atol=1e-12)


def test_gene_consistentemente_sub_movido_sai_acima_de_um():
    """A ordem e preservada: mais sub-movido que a media da peso maior."""
    m = medida("a->b", GENES_PAINEL[:6], {"G0": 3.0, "G1": 1.0, "G2": 0.33})
    peso, _, _ = enc.agrupar([m], GENES_PAINEL)
    assert peso[0] > peso[1] > peso[2]
    assert peso[0] > 1.0
    assert peso[2] < 1.0


def test_recorte_relativo_barra_o_estouro():
    """Gene cujo preditor quase nao se move teria `k_g` estourado por divisao.

    O teto aplicado NAO e `TETO_RELATIVO`: o recorte vem antes do encolhimento e
    os dois compoem em `TETO_RELATIVO ** (n/(n+K))`. Com n = 1.000 isso da 7,25.
    Este teste fixa o valor composto, porque foi ele que a primeira versao errou
    — o comentario do modulo prometia 8x e o codigo aplicava 7,25x.
    """
    entradas = 1000.0
    m = medida("a->b", GENES_PAINEL[:6], {"G0": 500.0}, entradas=entradas)
    peso, _, _ = enc.agrupar([m], GENES_PAINEL)
    teto = enc.TETO_RELATIVO ** (entradas / (entradas + enc.K_ENCOLHIMENTO))
    # A normalizacao desloca tudo, entao o que se confere e a razao entre o
    # recortado e um gene neutro.
    assert peso[0] / peso[1] == pytest.approx(teto, rel=1e-9)
    assert teto < enc.TETO_RELATIVO


def test_encolhimento_morde_o_gene_pouco_medido():
    """Com poucas entradas o gene vai para o global da direcao, nao para o k_g dele."""
    muitas = medida("a->b", GENES_PAINEL[:6], {"G0": 4.0}, entradas=100_000.0)
    poucas = medida("a->b", GENES_PAINEL[:6], {"G0": 4.0}, entradas=5.0)
    peso_muitas, _, _ = enc.agrupar([muitas], GENES_PAINEL)
    peso_poucas, _, _ = enc.agrupar([poucas], GENES_PAINEL)
    # Mesmo k_g bruto, mas o de 5 entradas fica muito mais perto do neutro.
    assert peso_muitas[0] / peso_muitas[1] > 3.0
    assert peso_poucas[0] / peso_poucas[1] < 1.3


def test_direcao_sem_alvo_profundo_nao_entra():
    """O gate do alvo e o do `politica_de_beta`: so repo profundo serve de denominador."""
    direcoes = enc.direcoes_com_alvo_profundo()
    if not direcoes:
        pytest.skip("arrays do PIE ausentes")
    for _, _, repo_alvo, _, perts, _ in direcoes:
        assert repo_alvo in enc.REPOS_ALVO_PROFUNDO
        assert len(perts) >= enc.MIN_PERTS_ALVO


def test_spearman_concorda_com_o_posto():
    """O idioma de posto do `transfer_fit`, e ele tem de ser invariante a monotona."""
    rng = np.random.default_rng(3)
    a = rng.normal(size=400)
    b = 2.0 * a + rng.normal(scale=0.3, size=400)
    assert enc.spearman(a, b) == pytest.approx(enc.spearman(a, np.exp(b)), abs=1e-12)
    assert enc.spearman(a, a) == pytest.approx(1.0, abs=1e-12)
    assert np.isnan(enc.spearman(a, np.zeros_like(a)))
